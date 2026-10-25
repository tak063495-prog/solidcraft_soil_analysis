"""Analytic polynomial/volume checks for two actual nonmatching CAD meshes."""
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import time
from tet10_overlap import Overlay, np
import freeze_sources


EXPONENTS = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [2, 0, 0],
                      [0, 2, 0], [0, 0, 2], [1, 1, 0], [0, 1, 1], [1, 0, 1]])


def moment(exponents, retained):
    a, b, c = map(int, exponents)

    def integral(lo, hi, ylo, yhi, zlo):
        top = math.fsum(math.comb(c+1, k)*10**(c+1-k)*(-1/3)**k*
                        (hi**(a+k+1)-lo**(a+k+1))/(a+k+1) for k in range(c+2))
        top -= zlo**(c+1)*(hi**(a+1)-lo**(a+1))/(a+1)
        return (yhi**(b+1)-ylo**(b+1))/(b+1)*top/(c+1)
    full = integral(0., 12., 0., 8., 0.)
    return full-integral(4., 8., 2., 6., 7.) if retained else full


def analytic_squared(coefficients, retained):
    return math.fsum(float(coefficients[i]@coefficients[j])*moment(EXPONENTS[i]+EXPONENTS[j], retained)
                     for i in range(10) for j in range(10))


def field(xyz):
    x, y, z = np.asarray(xyz).T
    return np.column_stack((x*x+.3*y-.4*z*z, .2*x*y, 1-.1*z+.03*x*z))


def delta(xyz):
    x, y, z = np.asarray(xyz).T
    return np.column_stack((.01*x*y, .02*y*z, .03*x*x-.01*z*z))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(mesh_a, mesh_b, output):
    output.mkdir(parents=True, exist_ok=False)
    freeze_sources.freeze(output)
    bindings = {str(p): sha(p) for p in (mesh_a, mesh_b)}
    a, b = [json.loads(p.read_text(encoding='utf-8-sig')) for p in (mesh_a, mesh_b)]
    coefficient_a, coefficient_delta = np.zeros((10, 3)), np.zeros((10, 3))
    coefficient_a[4, 0], coefficient_a[2, 0], coefficient_a[6, 0] = 1., .3, -.4
    coefficient_a[7, 1] = .2
    coefficient_a[0, 2], coefficient_a[3, 2], coefficient_a[9, 2] = 1., -.1, .03
    coefficient_delta[7, 0], coefficient_delta[8, 1] = .01, .02
    coefficient_delta[4, 2], coefficient_delta[6, 2] = .03, -.01
    records = []
    for order, reversed_meshes in ((3, False), (4, True)):
        start = time.monotonic()
        left, right = (b, a) if reversed_meshes else (a, b)
        overlay = Overlay(left, right, order,
                          progress=lambda e, n, pairs: print(f'Overlay order {order}: {e}/{n}, intersections {pairs}', flush=True))
        u_left = field(left['nodes'])+(delta(left['nodes']) if reversed_meshes else 0)
        u_right = field(right['nodes'])+(0 if reversed_meshes else delta(right['nodes']))
        for retained in (False, True):
            active_left, active_right = [np.ones(len(m['tetrahedra']), dtype=bool) for m in (left, right)]
            if retained:
                active_left[left['excavation_sets']['Pit']] = False
                active_right[right['excavation_sets']['Pit']] = False
            result = overlay.compare(u_left, u_right, active_left, active_right)
            expected = analytic_squared(coefficient_delta, retained)
            expected_reference = analytic_squared(coefficient_a if reversed_meshes else coefficient_a+coefficient_delta, retained)
            error = abs(result['difference_squared_integral_m5']-expected)/expected
            reference_error = abs(result['reference_squared_integral_m5']-expected_reference)/expected_reference
            assert error < 1e-10 and reference_error < 1e-10, (error, reference_error)
            assert abs(result['active_volume_m3']-(752 if retained else 768)) < 1e-8
            zero = overlay.compare(field(left['nodes']), field(right['nodes']), active_left, active_right)
            assert zero['rms_difference_m'] < 5e-11, zero
            constant = overlay.compare(np.ones_like(u_left)*[1, 2, 3], np.zeros_like(u_right), active_left, active_right)
            assert abs(constant['rms_difference_m']-np.sqrt(14)) < 1e-11
            records.append(dict(quadrature_order=order, mesh_order_reversed=reversed_meshes,
                                retained_domain=retained, result=result,
                                analytic_squared_difference_m5=expected,
                                analytic_relative_error=error, analytic_reference_relative_error=reference_error,
                                identical_quadratic_field_rms_error_m=zero['rms_difference_m'],
                                constant_difference_rms_error_m=abs(constant['rms_difference_m']-np.sqrt(14)),
                                geometry=overlay.metadata, elapsed_seconds=time.monotonic()-start))
        refused = False
        bad = np.ones(len(left['tetrahedra']), dtype=bool)
        bad[left['excavation_sets']['Pit']] = False
        try:
            overlay.compare(u_left, u_right, bad, np.ones(len(right['tetrahedra']), dtype=bool))
        except ValueError as error:
            refused = 'same active physical domains' in str(error)
        assert refused
    curved = copy.deepcopy(a)
    curved['nodes'][curved['tetrahedra10'][0][4]][2] += .01
    curved_refused = False
    try:
        Overlay(curved, b)
    except ValueError as error:
        curved_refused = 'curved cells are unsupported' in str(error)
    assert curved_refused
    coordinates = copy.deepcopy(b)
    coordinates['origin_m'][0] += 1
    coordinate_refused = False
    try:
        Overlay(a, coordinates)
    except ValueError as error:
        coordinate_refused = 'coordinate metadata differ' in str(error)
    assert coordinate_refused
    assert bindings == {str(p): sha(p) for p in (mesh_a, mesh_b)}
    assert json.loads((output/'source-snapshot.json').read_text()) == freeze_sources.manifest()
    report = dict(schema='soil-fem-affine-tet10-overlap-verification/1', passed=True, analytic_cases=records,
                  changed_active_domain_refused=True, curved_geometry_refused=True, coordinate_change_refused=True,
                  input_sha256=bindings, goal_complete=False, p6_complete=False, spatial_convergence_verified=False,
                  limitations='Affine geometric displacement integration only; not raw constitutive-field interpolation or nonlinear FEM acceptance')
    with (output/'verification.json').open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
    print(json.dumps(dict(passed=True, analytic_cases=len(records),
                          maximum_analytic_relative_error=max(r['analytic_relative_error'] for r in records),
                          maximum_identical_field_rms_m=max(r['identical_quadratic_field_rms_error_m'] for r in records))), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ('mesh_a', 'mesh_b', 'output'):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    verify(*(getattr(args, k).resolve() for k in ('mesh_a', 'mesh_b', 'output')))
