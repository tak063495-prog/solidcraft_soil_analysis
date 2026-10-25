"""Actual SolidCraft two-material curved Tet10 mesh under transverse traction.

All bodies remain active. This is finite-bulk mixed elasticity, not a mixed
plastic excavation or a general inf-sup certification. Geometry is reused from
the accepted CAD export; only quadrature order and the explicit upper-layer
Poisson ratio differ from the original general-elastic test.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import numpy as np
import fem_reference as ref
from verify_mixed_mechanics import solve_reference, verify


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text(encoding='utf-8'))
def relative(a,b):return float(np.linalg.norm(np.asarray(a)-b)/max(np.linalg.norm(b),1e-30))


def run(solver,out):
    root=Path(__file__).resolve().parents[1];solver=solver.resolve()
    parent=root/'runs/tet10-boundary-test-20261023-224400-66dd339e'
    source=parent/'curved/ground';mp=source/'mesh.json';mesh=read(mp)
    accepted=parent/'curved/accepted-stage-run.json';parent_record=read(accepted)
    assert sha(mp)==parent_record['inputs_and_results_sha256'][str(mp)]
    bodies=read(source/'bodies.json');provenance=read(source/'mesh.provenance.json')
    assert provenance['cad_manifest_sha256']==sha(source/'bodies.json')
    for body in bodies['regions']+bodies.get('excavation_tools',[]):
        assert Path(body['step']).name==body['step'] and sha(source/body['step'])==body['source_sha256']
    by_body={b['body']:b for b in bodies['regions']}
    for region in mesh['regions']:assert region['source_sha256']==by_body[region['body']]['source_sha256']
    assert len(mesh['nodes'])==595 and len(mesh['tetrahedra10'])==312
    xyz=np.asarray(mesh['nodes']);curvature=0.
    for cell in mesh['tetrahedra10']:
        for k,(a,b) in enumerate(ref.EDGES):curvature=max(curvature,float(np.linalg.norm(xyz[cell[4+k]]-(xyz[cell[a]]+xyz[cell[b]])/2)))
    assert curvature>1e-6,'this test must actually contain curved edges'
    inputs=[solver,mp,accepted,source/'bodies.json',source/'mesh.provenance.json',parent/'high-order/mechanics-case.json',
            Path(__file__),root/'tools/verify_mixed_mechanics.py',Path(ref.__file__)]
    fingerprints={str(p):sha(p) for p in inputs}
    out.mkdir(exist_ok=False);geometry=out/'cad-inputs';geometry.mkdir()
    for p in (mp,source/'bodies.json',source/'mesh.provenance.json',accepted):shutil.copyfile(p,geometry/p.name)
    for body in bodies['regions']+bodies.get('excavation_tools',[]):shutil.copyfile(source/body['step'],geometry/body['step'])
    case=read(parent/'high-order/mechanics-case.json');case['materials'][1]['poisson']=.4999
    cp=out/'case.json'
    with cp.open('x') as f:json.dump(case,f,indent=2)
    results=[]
    for order in (3,4,5):
        m=json.loads(json.dumps(mesh));m['quadrature_order']=order
        path=out/f'mesh-q{order}.json';result=out/f'result-q{order}.json'
        with path.open('x') as f:json.dump(m,f)
        subprocess.run([sys.executable,str(root/'tools/run_logged.py'),str(out/f'Rust-q{order}.log'),str(solver),
            'solve-mixed-mechanics',str(path),str(cp),str(result)],cwd=root,check=True)
        results.append(read(result))
    m=read(out/'mesh-q5.json');references=[]
    for order in (8,9):
        expected=solve_reference(m,case,order)
        # Preserve the independent physical solution, not only a pass flag.
        raw={k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in expected.items()}
        with (out/f'reference-gl{order}.json').open('x') as f:json.dump(raw,f)
        verify(out/'mesh-q5.json',cp,out/'result-q5.json',solver,out/f'independent-gl{order}.json',order)
        references.append(expected)
    refinements=[relative(a['displacement_m'],b['displacement_m']) for a,b in zip(results,results[1:])]
    assert refinements[1]<.6*refinements[0] and refinements[1]<1e-5,refinements
    changes={key:relative(references[0][key],references[1][key]) for key in ('displacement','pressure','reaction')}
    assert changes['displacement']<1e-7 and changes['pressure']<1e-6 and changes['reaction']<1e-6,changes
    assert references[-1]['pressure_nullity']==0,references[-1]['pressure_nullity']
    # Audit the native integration inventory against independent geometry.
    count=0;layout=results[-1]['integration_points']
    for element in range(len(m['tetrahedra'])):
        q=ref.volume(m,element)
        for local,weight in enumerate(q['weight']):
            point=layout[count];count+=1
            assert point['element_index']==element and point['local_index']==local
            np.testing.assert_allclose(point['position_m'],q['position'][local],rtol=1e-12,atol=1e-12)
            assert abs(point['weight_m3']-weight)<1e-12
    assert count==len(layout)
    assert fingerprints=={str(p):sha(p) for p in inputs},'source/input/binary changed during test'
    record=dict(passed=True,goal_complete=False,p6_complete=False,
        scope='Actual curved SolidCraft two-layer mesh, finite-bulk mixed elasticity, upper nu=0.4999, transverse traction; all bodies active',
        nodes=595,elements=312,maximum_edge_midpoint_curvature_m=curvature,
        rust_quadrature_displacement_changes=refinements,reference_gl8_gl9_changes=changes,
        independent_all_field_comparisons=[read(out/f'independent-gl{n}.json') for n in (8,9)],
        pressure_coupling_nullity=references[-1]['pressure_nullity'],pressure_inf_sup_positive_beta=references[-1]['pressure_inf_sup_positive_beta'],
        all_native_integration_owners_positions_weights_verified=True,input_and_source_sha256=fingerprints,
        limitations='No plasticity, excavation, native UI, spatial convergence or general inf-sup claim',
        evidence_sha256={p.relative_to(out).as_posix():sha(p) for p in out.rglob('*') if p.is_file()})
    with (out/'verification.json').open('x') as f:json.dump(record,f,indent=2)
    print(json.dumps({k:record[k] for k in ('passed','scope','rust_quadrature_displacement_changes','reference_gl8_gl9_changes','pressure_coupling_nullity')}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('solver',type=Path);p.add_argument('output',type=Path);a=p.parse_args();run(a.solver,a.output.resolve())
