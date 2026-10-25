"""Locate spatial error as the exact MC elastic/plastic front crosses a Tet10.

This is a diagnostic for the existing column, not an alternative acceptance
criterion. All original final-state stress/plastic tolerances remain in force.
"""
import argparse
import copy
import csv
import json
from pathlib import Path
import numpy as np
from verify_2d_column_reference import read, sha, maximum, input_subset, serialized_identity


def exact(z, pressure):
    z = np.asarray(z)
    young, nu, gamma, height, cohesion = 1000., .25, 20., .2, 2.
    sine, cosine = np.sin(np.deg2rad(10.)), np.cos(np.deg2rad(10.))
    ratio = nu/(1-nu)
    yield_vertical = 2*cohesion*cosine/(1-ratio-(1+ratio)*sine)
    front = np.clip(height+(pressure-yield_vertical)/gamma, 0., height)
    bulk = young/(3*(1-2*nu))
    modulus = young*(1-nu)/((1+nu)*(1-2*nu))
    a, b = 3*(1+sine)/(3-sine), 4*cohesion*cosine/(3-sine)
    sv = pressure+gamma*(height-z)
    plastic = sv > yield_vertical
    sh = np.where(plastic, (3*(sv-b)/a-sv)/2, ratio*sv)
    epv = np.where(plastic, (sv-b)/(a*bulk)-(sv-2*nu*sh)/young, 0.)
    history = np.column_stack((-.5*epv, -.5*epv, epv, np.zeros((len(z), 3))))
    stress = np.column_stack((sh, sh, sv, np.zeros((len(z), 3))))
    low = np.minimum(z, front)
    up = np.maximum(z, front)
    # Integrate strain from the fixed base, split at the moving yield front.
    u = -((pressure+gamma*height-b)*low-gamma*low*low/2)/(a*bulk)
    u -= ((pressure+gamma*height)*(up-front)-gamma*(up*up-front*front)/2)/modulus
    return u, stress, history, float(front), float(yield_vertical)


def prepare(base, out):
    out.mkdir(exist_ok=False)
    mesh = read(base/'mesh-3d.json')
    case = copy.deepcopy(read(base/'3d-n2-fine/stages.json'))
    case['stages'] = [dict(kind='boundary', name='partial_yield_front', target_factor=.5)]
    for name, value in (('mesh.json', mesh), ('stages.json', case)):
        with (out/name).open('x', encoding='utf-8') as stream:
            json.dump(value, stream, indent=2)
    print('Prepared unchanged MC column with 6 kPa final load and an internal yield front.', flush=True)


def assess(base, out):
    result = read(out/'result.json')
    mesh = read(out/'mesh.json')
    assert result['complete'] and len(result['frames']) == 2
    input_subset(read(out/'stages.json'), result['case'])
    input_subset(read(base/'ground-3d.json'), result['ground'])
    for key in ('ground', 'case'):
        assert serialized_identity(out/'result.json', key) == result['provenance'][key+'_sha256']
    assert result['provenance']['solver_binary_sha256'] == sha(Path(__file__).resolve().parents[1]/'target/release/soil-fem-3d.exe')
    frame = result['frames'][-1]
    assert frame['relative_force_residual'] < 1e-8
    assert abs(frame['support_reaction_kn'][2]-2.) < 1e-7
    xyz = np.asarray(mesh['nodes'])
    points = np.asarray([p['position_m'] for p in result['integration_points']])
    u, _, _, front, threshold = exact(xyz[:, 2], 6.)
    ug, _, _, _, _ = exact(xyz[:, 2], 0.)
    _, stress, history, _, _ = exact(points[:, 2], 6.)
    actual_u = np.asarray(frame['state']['displacement_m'])
    actual_s = np.asarray([s['effective_stress_compression_positive_kpa'] for s in frame['state']['material_states']])
    actual_ep = np.asarray([s['plastic_strain_compression_positive'] for s in frame['state']['material_states']])
    # Independently integrate the analytic strain on either side of the front
    # and check the derivative of its displacement at non-switching points.
    test = np.linspace(.001, .199, 103)
    test = test[abs(test-front) > 1e-4]
    h = 1e-7
    plus = exact(test+h, 6.)[0]
    minus = exact(test-h, 6.)[0]
    _, sig, ep, _, _ = exact(test, 6.)
    strain = (sig[:, 2]-.25*(sig[:, 0]+sig[:, 1]))/1000.+ep[:, 2]
    derivative_error = maximum(-(plus-minus)/(2*h), strain)
    assert derivative_error < 1e-10
    errors = np.abs(actual_s-stress)
    worst = np.unravel_index(np.argmax(errors), errors.shape)
    report = dict(schema='soil-fem-mc-moving-yield-front-diagnostic/1', diagnostic_completed=True,
        goal_complete=False, p6_complete=False, final_column_acceptance=False,
        scope='One partially yielding confined column; original final-load acceptance bounds unchanged',
        pressure_kpa=6., analytic_yield_vertical_stress_kpa=threshold,
        analytic_front_z_m=front, front_inside_first_coarse_cell=0. < front < .1,
        analytic_displacement_derivative_check_maximum=derivative_error,
        maximum_displacement_error_m=maximum(actual_u[:, 2], u-ug),
        maximum_stress_error_kpa=maximum(actual_s, stress),
        maximum_plastic_strain_error=maximum(actual_ep, history),
        maximum_shear_stress_kpa=float(abs(actual_s[:, 3:]).max()),
        worst_stress_point_index=int(worst[0]), worst_stress_component=int(worst[1]),
        worst_stress_point_position_m=points[worst[0]].tolist(),
        relative_force_residual=frame['relative_force_residual'],
        interpretation='A front lies inside the coarse element. Endpoint errors alone cannot prove their cause or spatial convergence.',
        source_sha256=sha(Path(__file__)),
        evidence_sha256={p.name: sha(p) for p in out.iterdir() if p.is_file()}, ground_sha256=sha(base/'ground-3d.json'))
    with (out/'diagnostic.json').open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2)
    print(json.dumps(report), flush=True)


def prepare_excel(base, out):
    for ny in (2, 4):
        fixture = read(base/f'ny{ny}/fixture-2d.json')
        for load in fixture['loads']:
            load[4] *= .5
        directory = out/f'ny{ny}'
        directory.mkdir(exist_ok=False)
        with (directory/'fixture-2d.json').open('x', encoding='utf-8') as stream:
            json.dump(fixture, stream, indent=2)
    print('Prepared 2D Q8 2/4-element partial-front fixtures with 6 kPa load.', flush=True)


def assess_excel(base, out):
    cases = []
    for ny in (2, 4):
        directory = out/f'ny{ny}'
        result = read(directory/'excel/result.json')
        fixture = read(directory/'fixture-2d.json')
        assert result['status'] == 'PASS' and result['application_code_unchanged']
        assert sha(directory/'fixture-2d.json') == result['fixture_sha256'].lower()
        assert sha(base/'source/2DSoilFEM_ui_20261024.xlsm') == result['source_workbook_sha256'].lower()
        assert not read(directory/'excel/code-comparison.json')['program_changes']
        for table, key in (('節点データ', 'nodes'), ('要素データ', 'elements'), ('材料データ', 'materials'), ('載荷', 'loads')):
            assert result['tables'][table][1:] == fixture[key]
        folder = directory/'excel/case_out/stg02_APPLY'
        def csv_rows(name):
            with (folder/(name+'.csv')).open(encoding='utf-8-sig', newline='') as stream:
                return list(csv.DictReader(stream))
        nodes, gauss, cumulative = [csv_rows(n) for n in ('node', 'gauss', 'plastic_cumulative')]
        assert len(nodes) == len(fixture['nodes']) and len(gauss) == len(cumulative) == 4*ny
        z = np.array([float(r['y']) for r in nodes])
        gp_z = np.array([((int(r['elem'])-1)+(1+(-1 if int(r['gp']) < 2 else 1)*.5773502692)/2)*.2/ny for r in gauss])
        u, _, _, front, _ = exact(z, 6.)
        _, stress, history, _, _ = exact(gp_z, 6.)
        actual_u = np.array([float(r['uy']) for r in nodes])
        actual_s = np.array([[float(r[k]) for k in ('sx', 'sz', 'sy')]+[0., 0., float(r['txy'])] for r in gauss])
        lookup = {(int(r['elem']), int(r['gp'])): r for r in cumulative}
        actual_ep = np.array([[float(lookup[(int(r['elem']), int(r['gp']))][k]) for k in ('ex_pl_cumulative', 'ez_pl_cumulative', 'ey_pl_cumulative')]+[0., 0., float(lookup[(int(r['elem']), int(r['gp']))]['gxy_pl_cumulative'])] for r in gauss])
        assert abs(result['tables']['ステージ結果'][2][7]-2.) < 1e-8
        # Preserve observations with no invented looser field-accuracy gate.
        cases.append(dict(vertical_elements=ny, analytic_front_z_m=front,
            maximum_displacement_error_m=maximum(actual_u, u),
            maximum_stress_error_kpa=maximum(actual_s, stress),
            maximum_plastic_strain_error=maximum(actual_ep, history)))
    report = dict(schema='soil-fem-actual-2d-moving-front-diagnostic/1', diagnostic_completed=True,
        goal_complete=False, p6_complete=False, final_column_acceptance=False,
        cases=cases, displacement_error_refinement_ratio=cases[1]['maximum_displacement_error_m']/cases[0]['maximum_displacement_error_m'],
        source_sha256=sha(Path(__file__)),
        evidence_sha256={p.relative_to(out).as_posix(): sha(p) for p in out.rglob('*') if p.is_file()})
    with (out/'excel-diagnostic.json').open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2)
    print(json.dumps({k: report[k] for k in ('diagnostic_completed', 'cases', 'displacement_error_refinement_ratio')}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=['prepare', 'assess', 'prepare-excel', 'assess-excel'])
    parser.add_argument('base', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    {'prepare': prepare, 'assess': assess, 'prepare-excel': prepare_excel,
     'assess-excel': assess_excel}[args.command](args.base.resolve(), args.output.resolve())
