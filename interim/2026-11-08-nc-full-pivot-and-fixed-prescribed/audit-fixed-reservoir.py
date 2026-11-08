"""Input-derived 1D head/flux control for the existing axis-aligned small column."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tools'))
import verify_wet_structures as wet


def audit(folder):
    mesh = wet.read(folder/'mesh.json'); ground = wet.read(folder/'ground.json')
    case = wet.read(folder/'stages.json'); native = wet.read(folder/'stage-results.json')
    xyz = np.array(mesh['nodes']); z = xyz[:, 2]
    policy = case['joint_system']['joints'][0]['water']
    assert policy['kind'] == 'prescribed'
    reservoir = policy['level_m']
    assert reservoir == ground['initialization']['water_level_m']
    flow = next(s for s in case['stages'] if s['kind'] == 'steady_seepage')['flow']
    heads = {h['face_set']: h['head_m'] for h in flow['heads']}
    assert set(heads) == {'bottom', 'top'}
    materials = {m['id']: np.array(m['hydraulic_conductivity_m_s']) for m in flow['materials']}
    seam = np.unique(xyz[np.array(mesh['face_sets'][case['joint_system']['joints'][0]['side_a']]).ravel(), 2])
    assert len(seam) == 1
    zs = float(seam[0]); lo = float(z.min()); hi = float(z.max())
    assert lo < zs < hi
    max_head = max_flux = max_pressure = 0.; gradients = {}; fluxes = {}
    for point in native['flow_audits'][0]['solution']['points']:
        ei = point['element_index']; zz = point['position_m'][2]
        lower = max(z[np.array(mesh['tetrahedra'][ei])]) <= zs
        a, b, ha, hb = (lo, zs, heads['bottom'], reservoir) if lower else (zs, hi, reservoir, heads['top'])
        assert a <= zz <= b
        gradient = (hb-ha)/(b-a)
        tensor = materials[mesh['material_ids'][ei]]
        assert np.array_equal(tensor, np.eye(3)*tensor[2, 2])
        head = ha+(zz-a)*gradient
        flux = np.array([0., 0., -tensor[2, 2]*gradient])
        max_head = max(max_head, abs(point['head_m']-head))
        max_pressure = max(max_pressure, abs(point['pressure_head_m']-(head-zz)))
        max_flux = max(max_flux, float(np.max(abs(np.array(point['darcy_flux_m_s'])-flux))))
        gradients['lower' if lower else 'upper'] = gradient
        fluxes['lower' if lower else 'upper'] = flux.tolist()
    assert max_head < 1e-8 and max_pressure < 1e-8 and max_flux < 1e-12
    pressure = ground['initialization']['water_unit_weight_kn_m3']*(reservoir-zs)
    for frame in native['frames']:
        assert np.max(abs(np.array(frame['state']['joint_water']['physical_pressure_kpa'])-pressure)) < 1e-9
    return dict(schema='soil-fem-fixed-reservoir-analytic-control/1', passed=True,
                scope='existing axis-aligned column, two isotropic layers, fixed PRESCRIBED reservoir only',
                reservoir_head_m=reservoir, seam_z_m=zs, input_derived_head_gradients=gradients,
                input_derived_darcy_flux_m_s=fluxes, input_derived_joint_pressure_kpa=pressure,
                maximum_head_error_m=max_head, maximum_pressure_head_error_m=max_pressure,
                maximum_flux_error_m_s=max_flux, p6_complete=False, nonzero_soil_policy_validated=False)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--evidence', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    with a.output.open('x', encoding='utf-8') as f:
        json.dump(audit(a.evidence), f, indent=2)
