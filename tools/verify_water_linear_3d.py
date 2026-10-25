"""Input-only direct 3D elastic solution of the actual SolidCraft water case."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from water_reference import linear_3d
from verify_high_order import prepare
import fem_reference as fem


def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,d):
    with p.open('x') as f:json.dump(d,f,indent=2)


def verify(solver,run,output):
    output.mkdir(parents=True,exist_ok=False);mesh=read(run/'ground/mesh.json');ground=read(run/'ground/ground.json');case=read(run/'stages.json');native=read(run/'stage-results.json')
    assert native['complete'] and native['provenance']['solver_binary_sha256']==sha(solver)
    assert not case.get('boundary') and not case.get('beam_system') and not case.get('anchor_system') and not case.get('joint_system')
    q,offset=prepare(mesh,native['integration_points']);active=np.ones(len(mesh['tetrahedra']),bool);pressure=ground['initialization']['surface_surcharge_kpa']
    pore_fields=None
    if native['schema']=='soil-fem-stage-result/12':
        from darcy_reference import stage_pressures
        pore_fields,_=stage_pressures(mesh,ground,native,q,offset)
    h=ground['initialization']['water_level_m'];h=min(p[2] for p in mesh['nodes']) if h is None else h;rows=[]
    for index,frame in enumerate(native['frames']):
        if index:
            stage=case['stages'][index-1]
            if stage['kind']=='hydrostatic_water':h=stage['target_level_m']
            elif stage['kind']=='steady_seepage':h=stage['exterior_water_level_m']
            elif stage['kind']=='pressure':pressure=stage['target_pressure_kpa']
            elif stage['kind']=='excavate':
                pressure=stage['remaining_surface_pressure_kpa']
                for body in stage['bodies']:
                    if body in mesh.get('excavation_sets',{}):ids=mesh['excavation_sets'][body]
                    else:
                        rid=next(r['id'] for r in mesh['regions'] if r['body']==body);ids=np.flatnonzero(np.asarray(mesh['region_ids'])==rid)
                    active[ids]=False
            else:raise AssertionError('unsupported independent linear path')
        state=frame['state'];assert state['active_elements']==active.tolist() and abs(state['hydrostatic_water_level_m']-h)<1e-12
        reference=linear_3d(mesh,ground,case,active,pressure,h,None if pore_fields is None else pore_fields[index]);u=np.asarray(reference['displacement_m']);actual=np.asarray(state['displacement_m']);used=np.unique(fem.cells(mesh)[active])
        absolute=float(np.max(abs(actual[used]-u[used])));scale=float(np.linalg.norm(u[used]));relative=float(np.linalg.norm(actual[used]-u[used])/max(scale,1e-12))
        assert absolute<1e-8,(frame['name'],absolute)
        if scale>1e-10:assert relative<2e-6,(frame['name'],relative)
        stress=[];target=[];e=[];ee=[]
        for ei in np.flatnonzero(active):
            points=state['material_states'][offset[ei]:offset[ei+1]]
            assert all(p['equivalent_plastic_strain']==0. and p['plastic_work_kpa']==0. for p in points),'linear benchmark yielded'
            stress.extend(p['effective_stress_compression_positive_kpa'] for p in points);target.extend(reference['stress'][ei]);e.extend(p['void_ratio'] for p in points);ee.extend(reference['void_ratio'][ei])
        stress_error=float(np.linalg.norm(np.asarray(stress)-target)/np.linalg.norm(target));void_error=float(np.max(abs(np.asarray(e)-ee)))
        assert stress_error<2e-8 and void_error<1e-9,(frame['name'],stress_error,void_error)
        reference['stress']={str(k):v.tolist() for k,v in reference['stress'].items()};reference['void_ratio']={str(k):v.tolist() for k,v in reference['void_ratio'].items()}
        write(output/f'reference-{index:03}.json',reference)
        rows.append({'stage':frame['name'],'maximum_absolute_active_node_displacement_error_m':absolute,'relative_all_active_node_3d_displacement_error':relative if scale>1e-10 else None,
            'relative_all_active_point_stress_error':stress_error,'maximum_absolute_void_ratio_error':void_error,'reference_linear_force_residual':reference['true_linear_residual']})
        print('Independent staged-water 3D solve matched:',frame['name'],relative if scale>1e-10 else f'near-zero reference; absolute error {absolute:.3e} m',flush=True)
    report={'schema':'soil-fem-water-independent-linear-3d/1','passed':True,'solver_sha256':sha(solver),'frames':rows,
        'scope':'Actual SolidCraft Tet10 mesh and side rollers, independent input-only K0 fields, analytic constant elastic tensor, NumPy full 3D direct solve and independent quadratic-ray partial-wet surface integration. Same declared material volume quadrature; spatial/quadrature convergence is checked separately. MC benchmark remains strictly elastic and is not a general 3D plastic oracle.',
        'inputs_sha256':{str(p):sha(p) for p in (run/'ground/mesh.json',run/'ground/ground.json',run/'stages.json',run/'stage-results.json')}}
    write(output/'verification.json',report)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--run',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();verify(a.solver.resolve(),a.run.resolve(),a.output_dir.resolve())
