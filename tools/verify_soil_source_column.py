"""Analytic control of the existing SOIL column with existing volume-source law.

The lower rectangular soil block has bottom head H0, an impermeable seam and
no-flow sides. div(-k grad H)=s gives H=H0+s*(L*t-t*t/2)/k,
qz=s*(t-L), t=z-z0. All values come from inputs, not native fields.
This only checks saturated Darcy and SOIL pressure transfer in this fixture.
"""


def _require(condition, message):
    if not condition:
        raise AssertionError(message)

import argparse
import json
from pathlib import Path

import numpy as np
import fem_reference as fem
from joint_water_reference import Water
import verify_wet_structures as wet


def finite(value, shape, name):
    try:
        out = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as error:
        raise ValueError('invalid GP '+name) from error
    if out.shape != shape or not np.isfinite(out).all():
        raise ValueError('finite GP '+name+' with shape '+str(shape)+' required')
    return out


def audit(mesh, ground, case, native):
    _require(native['complete'], "external validation failed: native['complete']")
    _require(mesh == wet.read(wet.ROOT/'examples/joint-water-column-tet10-mesh.json'), 'existing small geometry only')
    _require(case['joint_system']['structure_bodies'] == ['UpperClay'], "external validation failed: case['joint_system']['structure_bodies'] == ['UpperClay']")
    _require(len(case['joint_system']['joints']) == 1, "external validation failed: len(case['joint_system']['joints']) == 1")
    policy = case['joint_system']['joints'][0]
    _require(policy['water'] == dict(kind='soil', soil_side='a'), "external validation failed: policy['water'] == dict(kind='soil', soil_side='a')")
    water = Water(mesh, ground, case)
    xyz = np.asarray(mesh['nodes']); used = np.unique(fem.cells(mesh)[water.soil])
    z0 = float(xyz[used, 2].min()); z1 = float(xyz[used, 2].max()); length = z1-z0
    _require(length > 0, 'external validation failed: length > 0')
    seam = xyz[np.asarray(mesh['face_sets'][policy['side_a']]).ravel(), 2]
    _require(np.all(seam == z1), 'external validation failed: np.all(seam == z1)')
    rows = []
    stages = [(i, stage) for i, stage in enumerate(case['stages']) if stage['kind'] == 'steady_seepage']
    _require(len(stages) == len(native['flow_audits']) == 1, "external validation failed: len(stages) == len(native['flow_audits']) == 1")
    for (index, stage), flow in zip(stages, native['flow_audits'], strict=True):
        _require(flow['stage_index'] == index, "external validation failed: flow['stage_index'] == index")
        cfg = stage['flow']; _require(len(cfg['heads']) == 1, "external validation failed: len(cfg['heads']) == 1")
        boundary = cfg['heads'][0]
        _require(boundary['face_set'] == 'bottom' and boundary['head_m'] == 2., "external validation failed: boundary['face_set'] == 'bottom' and boundary['head_m'] == 2.")
        _require(not np.asarray(boundary.get('gradient', [0., 0., 0.])).any(), "external validation failed: not np.asarray(boundary.get('gradient', [0., 0., 0.])).any()")
        _require(not cfg.get('fluxes') and len(cfg['materials']) == 1, "external validation failed: not cfg.get('fluxes') and len(cfg['materials']) == 1")
        material = cfg['materials'][0]; k = material['hydraulic_conductivity_m_s'][2][2]
        _require(np.array_equal(material['hydraulic_conductivity_m_s'], k*np.eye(3)) and k > 0, "external validation failed: np.array_equal(material['hydraulic_conductivity_m_s'], k*np.eye(3)) and k > 0")
        source = material['volume_source_per_s']; _require(source == 2*k > 0, 'external validation failed: source == 2*k > 0')
        head0 = cfg['heads'][0]['head_m']; gamma = ground['initialization']['water_unit_weight_kn_m3']
        def head(z):
            t = z-z0
            return head0+source/k*(length*t-t*t/2)
        solution = flow['solution']; _require(solution['complete'], "external validation failed: solution['complete']")
        np.testing.assert_array_equal(solution['active_elements'], water.soil)
        field = flow['joint_flow']['original_field']
        np.testing.assert_array_equal(field['defined_nodes'], np.isin(np.arange(len(xyz)), used))
        head_error = float(np.max(abs(np.asarray(field['head_m'])[used]-head(xyz[used, 2]))))
        points = solution['points']
        if not isinstance(points, list) or not points:
            raise ValueError('nonempty GP inventory required')
        # Validate every point before reducing any errors: max(0, NaN) is zero.
        checked = [(point, finite(point['position_m'], (3,), 'position'),
                    finite(point['head_m'], (), 'head'),
                    finite(point['pressure_head_m'], (), 'pressure head'),
                    finite(point['darcy_flux_m_s'], (3,), 'flux')) for point in points]
        max_h = max_p = max_q = 0.
        for point, position, gp_head, gp_pressure, gp_flux in checked:
            _require(water.soil[point['element_index']], "external validation failed: water.soil[point['element_index']]")
            z = position[2]
            flux = np.array([0., 0., source*(z-z0-length)])
            max_h = max(max_h, float(abs(gp_head-head(z))))
            max_p = max(max_p, float(abs(gp_pressure-(head(z)-z))))
            max_q = max(max_q, float(np.max(abs(gp_flux-flux))))
        finite([head_error, max_h, max_p, max_q], (4,), 'maximum errors')
        # Reuse the established Darcy head/flux/pressure tolerances.
        _require(head_error < 1e-8 and max_h < 1e-8 and max_p < 1e-8 and max_q < 1e-12, 'external validation failed: head_error < 1e-8 and max_h < 1e-8 and max_p < 1e-8 and max_q < 1e-12')
        state = native['frames'][index+1]['state']; pressure = gamma*(head(z1)-z1)
        np.testing.assert_allclose(state['joint_water']['physical_pressure_kpa'], pressure, rtol=1e-10, atol=1e-9)
        ip = native['integration_points']; pore = np.asarray(state['pore_pressure_kpa'])
        expected = np.array([gamma*(head(p['position_m'][2])-p['position_m'][2]) if water.soil[p['element_index']] else 0. for p in ip])
        np.testing.assert_allclose(pore, expected, rtol=1e-10, atol=1e-8)
        volume = sum(p['weight_m3'] for p in solution['points']); total = source*volume
        _require(total > 0 and len(solution['balances']) == 1, "external validation failed: total > 0 and len(solution['balances']) == 1")
        balance = solution['balances'][0]
        for key in ('source_m3_s', 'absolute_source_m3_s', 'dirichlet_reaction_outflow_m3_s', 'gradient_boundary_outflow_m3_s'):
            np.testing.assert_allclose(balance[key], total, rtol=1e-7, atol=1e-14)
        _require(balance['neumann_outflow_m3_s'] == 0 and balance['relative_nodal_balance_error'] < 1e-9, "external validation failed: balance['neumann_outflow_m3_s'] == 0 and balance['relative_nodal_balance_error'] < 1e-9")
        bottom_out = sum(p['gradient_outflow_m3_s'] for p in solution['boundary_fluxes'] if p['boundary_kind'] == 'head')
        other = max(abs(p['gradient_outflow_m3_s']) for p in solution['boundary_fluxes'] if p['boundary_kind'] == 'no_flow')
        np.testing.assert_allclose(bottom_out, total, rtol=1e-7, atol=1e-14); _require(other < 1e-12, 'external validation failed: other < 1e-12')
        rows.append(dict(stage=stage['name'], source_per_s=source, soil_volume_m3=volume,
                         source_and_bottom_outflow_m3_s=total, head_at_seam_m=float(head(z1)),
                         physical_joint_pressure_kpa=float(pressure), bottom_flux_m_s=-source*length,
                         seam_flux_m_s=0., maximum_nodal_head_error_m=head_error,
                         maximum_gp_head_error_m=max_h, maximum_pressure_head_error_m=max_p,
                         maximum_flux_error_m_s=max_q, maximum_no_flow_face_outflow_m3_s=other))
    return dict(schema='soil-fem-soil-source-column-analytic/1', passed=True, frames=rows,
                scope='existing 12 Tet10 column; saturated SOIL source flow only; elastic closed sticking with 2 plates and 1 beam',
                p3b_p4_rows_complete=False, p6_complete=False)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--evidence', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(); folder = args.evidence
    report = audit(wet.read(folder/'mesh.json'), fem.normalize_ground(wet.read(folder/'ground.json')),
                   fem.normalize_stage_case(wet.read(folder/'stages.json')), wet.read(folder/'stage-results.json'))
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2, ensure_ascii=False)
