"""Independent 3D elastic control of a finite, perfectly bonded solid pile."""
import hashlib
import json
import os
from pathlib import Path
import sys
for name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu

import argparse
import re
parser=argparse.ArgumentParser()
parser.add_argument('run',type=Path)
parser.add_argument('--output',type=Path)
parser.add_argument('--mesh-labels',nargs='+',default=['coarse','fine'])
args=parser.parse_args()
labels=args.mesh_labels
if len(labels)>8 or len(set(labels))!=len(labels) or any(not re.fullmatch(r'[a-z][a-z0-9_-]{0,40}',name) for name in labels):
    raise ValueError('require 1 to8 distinct safe mesh labels')
folder=args.run.resolve()
root=Path(__file__).resolve().parents[1]
output=args.output.resolve() if args.output else folder/'independent-verification.json'
if output.exists():raise ValueError('independent output exists')
sys.path.insert(0,str(root/'tools'))
import fem_reference as ref
from verify_high_order import physical

def read(path):return json.loads(path.read_text(encoding='utf-8-sig'))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
solver_digest=read(folder/'solver-binding.json')['solver_sha256']
if sha(folder/'solver/soil-fem-3d.exe')!=solver_digest:
    raise ValueError('solver binary differs from frozen binding')
for label in labels:
    binding=read(folder/label/'input-binding.json')
    if binding['cad_manifest_sha256']!=sha(folder/'cad/bodies.json'):
        raise ValueError('CAD manifest differs from frozen binding')
    for name,digest in binding['files_sha256'].items():
        source=(folder/name).resolve()
        if not source.is_relative_to(folder) or sha(source)!=digest:
            raise ValueError('pile input differs from frozen binding: '+name)
    if read(folder/label/'stage-results.json')['provenance']['solver_binary_sha256']!=solver_digest:
        raise ValueError('result solver identity differs from frozen binding')
watched=[Path(__file__),root/'tools/fem_reference.py',root/'tools/verify_high_order.py',
    folder/'ground.json',folder/'stages.json',folder/'design.json',folder/'cad/bodies.json',folder/'solver/soil-fem-3d.exe']
watched += [folder/label/name for label in labels for name in ('mesh.json','stage-results.json','input-binding.json')]
fingerprints={str(p.resolve()):sha(p) for p in watched}
ground=ref.normalize_ground(read(folder/'ground.json'))
case=read(folder/'stages.json')
head_name=case['monitor_sets']['pile_head'];tip_name=case['monitor_sets']['pile_tip']
assert case['boundary']['surface_loads'][0]['face_set']==head_name
materials={m['id']:m for m in ground['materials']}
seeds=ground['initialization']['geostatic']['seed']['states']
assert set(materials)=={1,2} and all('mohr_coulomb' in m for m in materials.values())
assert ground['initialization']['water_level_m'] is None
assert all(s['kind']=='boundary' for s in case['stages'])
reports=[]
for label in labels:
    path=folder/label
    mesh=read(path/'mesh.json');result=read(path/'stage-results.json')
    assert result['complete'] and result['schema']=='soil-fem-stage-result/4' and len(result['frames'])==3
    serialized=result['case']['boundary']
    assert serialized['constraints']==case['boundary']['constraints']==[]
    assert serialized['nodal_loads']==case['boundary']['nodal_loads']==[]
    assert len(serialized['surface_loads'])==len(case['boundary']['surface_loads'])==1
    for key,value in case['boundary']['surface_loads'][0].items():
        assert serialized['surface_loads'][0][key]==value
    np.testing.assert_array_equal(serialized['surface_loads'][0]['gradient_kpa_per_m'],0.)
    np.testing.assert_array_equal(serialized['surface_loads'][0]['reference_point_m'],0.)
    xyz=np.asarray(mesh['nodes']);ndof=xyz.size;active=np.ones(len(mesh['tetrahedra']),bool)
    regions={r['id']:r for r in mesh['regions']}
    assert {r['body'] for r in regions.values()}=={'Soil','Pile'}
    assert all(all(f['state']['active_elements']) for f in result['frames'])
    rr=[];cc=[];vv=[];operators=[];seed_force=np.zeros(ndof);pile_force=np.zeros(ndof)
    independent_volume={1:0.,2:0.}
    for e,material_id in enumerate(mesh['material_ids']):
        law=materials[material_id]['mohr_coulomb'];d=ref.elastic_tensor(law['young_kpa'],law['poisson'])
        q=ref.volume(mesh,e,independent=True);indices=(q['nodes'][:,None]*3+np.arange(3)).reshape(-1)
        ke=np.einsum('qai,ab,qbj,q->ij',q['b'],d,q['b'],q['weight'],optimize=True)
        s0=seeds[regions[mesh['region_ids'][e]]['body']]
        force=-np.einsum('qai,a,q->i',q['b'],np.asarray(s0['effective_stress_compression_positive_kpa']),q['weight'])
        rr.extend(np.repeat(indices,len(indices)));cc.extend(np.tile(indices,len(indices)));vv.extend(ke.reshape(-1))
        seed_force[indices]+=force
        if material_id==2:
            pile_force[indices.reshape(-1,3)[:,2]]-=q['shape'].T@q['weight']*materials[2]['unit_weight_kn_m3']
        independent_volume[material_id]+=q['weight'].sum()
        operators.append((indices,ke,force,d,s0,ref.volume(mesh,e)))
    np.testing.assert_allclose([independent_volume[1],independent_volume[2]],[794.,6.],rtol=1e-10,atol=1e-9)
    k=coo_matrix((vv,(rr,cc)),shape=(ndof,ndof)).tocsc()
    basis,error=ref.staged_constraint_space(mesh,case,dict(active_elements=active,boundary_factor=0.,displacement_m=np.zeros_like(xyz)))
    assert error==0.
    tr=[];tc=[];tv=[];columns=0
    for node,block in basis.items():
        a,b=block.shape
        tr.extend(np.repeat(3*node+np.arange(a),b));tc.extend(np.tile(columns+np.arange(b),a));tv.extend(block.reshape(-1));columns+=b
    t=coo_matrix((tv,(tr,tc)),shape=(ndof,columns)).tocsc()
    reduced=(t.T@k@t).tocsc();lu=splu(reduced)
    external,_=physical(mesh,ground,active,0.,True)
    boundary=ref.conditions_load(mesh,case['boundary'],active)
    head_area=0.;head_weights=np.zeros(len(xyz));tip_weights=np.zeros(len(xyz))
    topology=ref.topology(mesh)
    for name,weights in ((head_name,head_weights),(tip_name,tip_weights)):
        for face in mesh['face_sets'][name]:
            owner=next(e for e in topology[tuple(sorted(face))] if mesh['material_ids'][e]==2)
            q=ref.face(mesh,owner,face);area=np.linalg.norm(q['area_vector'],axis=1)
            weights[q['nodes']]+=q['shape'].T@area
        np.testing.assert_allclose(weights.sum(),1.,atol=1e-11)
    np.testing.assert_allclose(boundary.sum(axis=0),[10.,0.,-100.],atol=1e-10)
    head_area=float(head_weights.sum())
    pile_interface=set()
    for face,owners in topology.items():
        if len(owners)==2 and {mesh['material_ids'][e] for e in owners}=={1,2}:
            owner=next(e for e in owners if mesh['material_ids'][e]==2)
            pile_interface.update(map(int,ref.face(mesh,owner,face)['nodes']))
    pile_interface=sorted(pile_interface)
    frame_reports=[];reference_solutions=[];interface_forces=[]
    for i,frame in enumerate(result['frames']):
        factor=0. if i==0 else case['stages'][i-1]['target_factor']
        force=(external+factor*boundary).reshape(-1)
        rhs=t.T@(force-seed_force);reduced_u=lu.solve(rhs);u=np.asarray(t@reduced_u)
        true=float(np.linalg.norm(reduced@reduced_u-rhs)/max(np.linalg.norm(rhs),1.))
        assert true<1e-10,true
        actual=frame['state'];assert actual['boundary_factor']==factor;np.testing.assert_array_equal(actual['active_elements'],active)
        relative=float(np.linalg.norm(np.asarray(actual['displacement_m']).reshape(-1)-u)/max(np.linalg.norm(u),1e-12))
        assert relative<1e-7,relative
        points=[];void=[];pile_residual=-pile_force-factor*boundary.reshape(-1)
        for e,(indices,ke,fe,d,s0,q) in enumerate(operators):
            eps=np.einsum('qai,i->qa',q['b'],u[indices])
            points.extend(np.asarray(s0['effective_stress_compression_positive_kpa'])-eps@d.T)
            void.extend((1+s0['void_ratio'])*np.exp(eps[:,:3].sum(axis=1))-1.)
            if mesh['material_ids'][e]==2:pile_residual[indices]+=ke@u[indices]+fe
        stored=actual['material_states']
        assert all(p['equivalent_plastic_strain']==0. and p['plastic_work_kpa']==0. for p in stored)
        np.testing.assert_array_equal([p['plastic_strain_compression_positive'] for p in stored],0.)
        points=np.asarray(points)
        stress_error=float(np.linalg.norm(np.asarray([p['effective_stress_compression_positive_kpa'] for p in stored])-points)/max(np.linalg.norm(points),1.))
        void_error=float(np.max(np.abs(np.asarray([p['void_ratio'] for p in stored])-void)))
        assert stress_error<1e-7 and void_error<1e-8,(stress_error,void_error)
        reaction=(k@u+seed_force-force).reshape(-1,3)
        np.testing.assert_allclose(reaction.sum(axis=0),frame['support_reaction_kn'],atol=1e-5,rtol=1e-9)
        transferred=pile_residual.reshape(-1,3)[pile_interface].sum(axis=0)
        expected_transfer=-(pile_force.reshape(-1,3).sum(axis=0)+factor*boundary.sum(axis=0))
        np.testing.assert_allclose(transferred,expected_transfer,atol=1e-5,rtol=1e-9)
        reference_solutions.append(u.reshape(-1,3));interface_forces.append(transferred)
        frame_reports.append(dict(frame=frame['name'],factor=factor,relative_original_matrix_residual=true,
            relative_all_node_displacement_error=relative,relative_all_gp_stress_error=stress_error,maximum_void_ratio_error=void_error,
            head_area_mean_displacement_m=(head_weights@u.reshape(-1,3)/head_area).tolist(),
            tip_area_mean_displacement_m=(tip_weights@u.reshape(-1,3)).tolist(),
            total_soil_force_on_pile_kn=transferred.tolist()))
    np.testing.assert_allclose(reference_solutions[0],reference_solutions[2],atol=1e-12,rtol=1e-12)
    delta=(head_weights@(reference_solutions[1]-reference_solutions[0])).tolist()
    rejected=[a for a in result['attempts'] if not a['accepted']]
    assert len(rejected)==2 and all(a['committed_state_unchanged_on_rejection'] for a in rejected)
    reports.append(dict(mesh=label,nodes=len(xyz),tet10=len(mesh['tetrahedra']),free_dofs=columns,
        head_area_mean_loading_increment_m=delta,frames=frame_reports,all_elastic_gp_fields_checked=True,
        shared_finite_interface_nodes=len(pile_interface),independently_integrated_volumes_m3=independent_volume,
        accepted_increments=sum(a['accepted'] for a in result['attempts']),rejected_increments=len(rejected)))
report=dict(schema='soil-fem-independent-solid-pile-control/1',passed=True,
    reference='Original CAD mesh, input seed and loads; independent GL125 quadrature, sparse LU, nodal SVD constraints; no Rust solved fields generate the solution',
    finite_section_single_pile=True,perfect_bond=True,pile_row_equivalence=False,
    reference_soil_force_convention='Weak consistent nodal interface force on pile, summed over all shaft/toe interface nodes. No arbitrary shaft/toe edge-force split.',
    meshes=reports,mesh_labels=labels,spatial_convergence_accepted=False,capacity_or_interface_slip_validated=False,
    source_sha256={str(p.resolve()):sha(p) for p in (Path(__file__),root/'tools/fem_reference.py',root/'tools/verify_high_order.py',folder/'ground.json',folder/'stages.json',folder/'design.json')},
    artifacts_sha256={str(p.relative_to(folder)):sha(p) for label in labels for p in (folder/label/'mesh.json',folder/label/'stage-results.json',folder/label/'input-binding.json')},
    input_and_source_sha256=fingerprints,solver_sha256=solver_digest,goal_complete=False,p6_complete=False)
if fingerprints!={str(p.resolve()):sha(p) for p in watched}:
    raise ValueError('pile inputs or verification sources changed during independent calculation')
with output.open('x',encoding='utf-8') as stream:json.dump(report,stream,indent=2,allow_nan=False)
print(json.dumps({k:v for k,v in report.items() if k not in ('source_sha256','artifacts_sha256')}),flush=True)
