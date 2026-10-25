"""Independent rate-form MCC equilibrium at a prescribed initialization fraction.

The holding/release force is reconstructed from the original seed and geometry.
No Rust solved state, displacement, pressure unknown or holding vector is read.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import mcc_rate_fem as rate
import fem_reference
import mcc_rate_midpoint
import verify_high_order
import high_order_io


def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def solve(mesh_file,ground_file,case_file,identity_file,output,fraction,steps):
    if not 0<fraction<=1 or steps<2:raise ValueError('fraction in (0,1] and at least two steps required')
    mesh,ground,case,identity=map(read,(mesh_file,ground_file,case_file,identity_file))
    cfg=ground['initialization'].get('geostatic')
    if cfg is None or cfg['initial_boundary_factor']!=0 or ground['initialization']['water_level_m'] is not None:
        raise ValueError('this diagnostic requires dry initialization with zero initial construction-boundary factor')
    if any(case.get(k) for k in ('joint_system','beam_system','anchor_system','initially_inactive_bodies')):
        raise ValueError('this independent diagnostic does not include members, contacts or inactive initial domains')
    sources=[mesh_file,ground_file,case_file,identity_file,Path(__file__),Path(rate.__file__),
             Path(fem_reference.__file__),Path(mcc_rate_midpoint.__file__),Path(verify_high_order.__file__),Path(high_order_io.__file__)]
    fingerprints={str(p):sha(p) for p in sources}
    active=np.ones(len(mesh['tetrahedra']),dtype=bool)
    system=rate.System(mesh,ground,active,case);initial=rate.seed(mesh,ground,identity)
    assert np.max(rate.terms(initial.reshape(-1,14),system.parameters)[0])<=1e-8
    start_force=system.force(initial.reshape(-1,14))
    physical=rate.physical(mesh,ground,active,ground['initialization']['surface_surcharge_kpa'],case.get('pressure_on_original_surface_only',False))[0].reshape(-1)
    change=fraction*(physical-start_force);target=start_force+change
    final,u,diagnostics=system.advance(initial.reshape(-1,14),np.zeros(np.asarray(mesh['nodes']).size),change,steps,
        progress=lambda t,n,e:print(f'Independent fractional initialization {fraction*t:.9f}, steps {n}, events {e}',flush=True))
    residual=float(np.linalg.norm(system.reduce(target-system.force(final)))/max(np.linalg.norm(target),1.))
    strain=-np.einsum('eqai,ei->eqa',system.b,u[system.dofs]);maximum=float(abs(strain).max())
    report=dict(schema='soil-fem-independent-mcc-geostatic-fraction/1',diagnostic_only=True,goal_complete=False,p6_complete=False,
        load_fraction=fraction,reference_steps=steps,
        state_origin='Original input stress/pc/e/plastic history and zero seed displacement; no solved Rust state/pressure/force initializes this calculation',
        method='Unchanged independent analytic continuum tangent and event-aware midpoint; independently reconstructed proportional holding-force release',
        displacement_m=u.reshape(-1,3).tolist(),material_state=final.reshape(initial.shape).tolist(),
        seed_material_state=initial.tolist(),active_elements=active.tolist(),target_nodal_force_kn=target.reshape(-1,3).tolist(),
        relative_force_residual=residual,maximum_engineering_strain_from_seed=maximum,
        maximum_normalized_yield=float(np.max(rate.terms(final,system.parameters)[0])),diagnostics=diagnostics,
        input_and_source_sha256=fingerprints,
        limitations='A completed prescribed fraction is not a completed geostatic initialization, excavation cycle, mixed formulation or P6 acceptance')
    assert residual<1e-8 and maximum<cfg['maximum_initial_engineering_strain']
    assert fingerprints=={str(p):sha(p) for p in sources},'reference source/input changed during computation'
    with output.open('x',encoding='utf-8') as f:json.dump(report,f,allow_nan=False)
    print(json.dumps({k:report[k] for k in ('load_fraction','reference_steps','relative_force_residual','maximum_engineering_strain_from_seed','maximum_normalized_yield')}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mesh',type=Path);p.add_argument('ground',type=Path);p.add_argument('case',type=Path)
    p.add_argument('identity',type=Path);p.add_argument('output',type=Path);p.add_argument('--fraction',type=float,required=True);p.add_argument('--steps',type=int,default=50)
    a=p.parse_args();solve(*(getattr(a,k).resolve() for k in ('mesh','ground','case','identity','output')),a.fraction,a.steps)
