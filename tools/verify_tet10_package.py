"""Basis, independent mechanics and curved-geometry/integration sensitivity."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import numpy as np
import fem_reference as ref
from verify_mechanics import verify as verify_mechanics
from export_mechanics import export as export_mechanics
from verify_high_order import physical,verify_stages
from export_stages import export


def geometry(run):
    mesh=json.loads((run/'ground/mesh.json').read_text());cells=ref.cells(mesh);xyz=np.asarray(mesh['nodes'])
    active=np.ones(len(cells),dtype=bool);active[mesh['excavation_sets']['Pit']]=False
    q=[ref.volume(mesh,i,independent=True) for i in range(len(cells))]
    volumes=np.asarray([v['weight'].sum() for v in q]);ground=json.loads((run/'ground/ground.json').read_text())
    _,surface=physical(mesh,ground,active,30.,True)
    curved=max(np.linalg.norm(xyz[c[4+k]]-(xyz[c[a]]+xyz[c[b]])/2) for c in cells for k,(a,b) in enumerate(ref.EDGES))
    return {'nodes':len(xyz),'elements':len(cells),'removed_volume_m3':float(volumes[~active].sum()),
        'removed_volume_relative_error':float(abs(volumes[~active].sum()/(8*np.pi)-1)),
        'cut_area_relative_error':float(abs(surface['cut_surface_area_m2']/(12*np.pi)-1)),
        'maximum_curved_midnode_offset_m':float(curved),**surface}


def verify(curved,finer,solver,output):
    output.mkdir(parents=True,exist_ok=False);root=Path(__file__).resolve().parents[1]
    report={'schema':'soil-fem-tet10-integration-geometry/1','basis':ref.verify_gmsh_basis()}
    mesh=json.loads((curved/'ground/mesh.json').read_text());ground_path=curved/'ground/ground.json'
    report['curved_geometry']=[geometry(curved),geometry(finer)]
    assert all(v['maximum_curved_midnode_offset_m']>1e-5 for v in report['curved_geometry'])
    assert report['curved_geometry'][1]['removed_volume_relative_error']<report['curved_geometry'][0]['removed_volume_relative_error']
    # Bending with two stiffnesses; a constant-strain patch alone cannot assess
    # quadrature error on non-affine mappings.
    case=json.loads((root/'examples/mechanics-pressure.json').read_text())
    case['materials'][1]['young_kpa']=2000.
    case['boundary']={'constraints':[{'node_set':'bottom','direction':d,'value_m':0.} for d in np.eye(3).tolist()],
        'surface_loads':[{'kind':'traction','face_set':'top','traction_kpa':[1.,.2,0.]}]}
    case_path=output/'mechanics-case.json';case_path.write_text(json.dumps(case,indent=2));solutions=[]
    for order in (3,4,5):
        mesh['quadrature_order']=order;mp=output/f'mesh-q{order}.json';mp.write_text(json.dumps(mesh))
        rp=output/f'mechanics-q{order}.json'
        subprocess.run([str(solver),'solve-mechanics',str(mp),str(case_path),str(rp)],check=True)
        solutions.append(json.loads(rp.read_text()))
        export_mechanics(mp,rp,output/f'mechanics-q{order}-vtu')
    reference8=verify_mechanics(output/'mesh-q5.json',case_path,output/'mechanics-q5.json',output/'mechanics-independent-order8.json',reference_order=8)
    reference9=verify_mechanics(output/'mesh-q5.json',case_path,output/'mechanics-q5.json',output/'mechanics-independent.json',reference_order=9)
    report['independent_order8_to_order9_displacement_difference']=float(np.linalg.norm(reference8[0]-reference9[0])/np.linalg.norm(reference9[0]))
    report['independent_order8_to_order9_stress_difference']=float(np.linalg.norm(reference8[1]-reference9[1])/np.linalg.norm(reference9[1]))
    assert report['independent_order8_to_order9_displacement_difference']<1e-8
    assert report['independent_order8_to_order9_stress_difference']<1e-7
    a=np.asarray(solutions[0]['displacement_m']);b=np.asarray(solutions[1]['displacement_m'])
    report['elastic_q3_to_q4_relative_displacement_difference']=float(np.linalg.norm(a-b)/np.linalg.norm(b))
    c=np.asarray(solutions[2]['displacement_m'])
    report['elastic_q4_to_q5_relative_displacement_difference']=float(np.linalg.norm(b-c)/np.linalg.norm(c))
    stage_path=output/'stages.json';stage_path.write_bytes((curved/'stages.json').read_bytes());rp=output/'stage-results.json'
    subprocess.run([str(solver),'solve-stages',str(output/'mesh-q4.json'),str(ground_path),str(stage_path),str(rp)],check=True)
    export(output/'mesh-q4.json',rp,output/'results');verify_stages(output/'mesh-q4.json',ground_path,rp,output/'results')
    q3=json.loads((curved/'stage-results.json').read_text());q4=json.loads(rp.read_text());a=np.asarray(q3['last_committed']['displacement_m']);b=np.asarray(q4['last_committed']['displacement_m'])
    used=np.unique(ref.cells(mesh)[np.asarray(q4['last_committed']['active_elements'])])
    report['mcc_q3_to_q4_relative_displacement_difference']=float(np.linalg.norm(a[used]-b[used])/np.linalg.norm(b[used]))
    report['mcc_q3_q4_monitors']=[r['frames'][-1]['monitors'] for r in (q3,q4)]
    report['interpretation']='Measured sensitivity only; two geometry sizes and finite quadrature comparisons do not prove spatial convergence, absence of locking, or field accuracy.'
    report['solver_sha256']=hashlib.sha256(solver.read_bytes()).hexdigest()
    report['inputs_sha256']={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in
        (curved/'ground/mesh.json',finer/'ground/mesh.json',curved/'stage-results.json',case_path,root/'tools/fem_reference.py',Path(__file__))}
    (output/'verification.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--curved-run',type=Path,required=True);p.add_argument('--fine-curved-run',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--solver',type=Path,default=Path('target/release/soil-fem-3d.exe'))
    a=p.parse_args();verify(a.curved_run.resolve(),a.fine_curved_run.resolve(),a.solver.resolve(),a.output_dir.resolve())
