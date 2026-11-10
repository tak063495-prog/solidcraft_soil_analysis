"""Input-only elastic reference for submerged, closed-sticking crack tips.

The independent water adapter constructs heads, pressure traces and exterior
loads from inputs. No native displacement, stress or pressure generates a solve.
"""
import copy
import json
from pathlib import Path
import numpy as np
import fem_reference as ref
from joint_water_reference import Water
from joint_tip_fixture import fixture
from verify_high_order import prepare, physical
from verify_joint_crack_tip import ROOT, write, solve, restart_and_rejection
from verify_joint_water import audit


def compare(mesh, ground, run):
    case=run['case'];water=Water(mesh,ground,case,run['joint_integration_points'])
    system=water.contact;xyz=np.asarray(mesh['nodes']);nd=xyz.size
    q,offset=prepare(mesh);positions=np.concatenate([v['position'] for v in q]);mat={m['id']:m['mohr_coulomb'] for m in ground['materials']}
    regions={r['id']:r['body'] for r in mesh['regions']}
    active=np.ones(len(mesh['tetrahedra']),bool);level=ground['initialization']['water_level_m'];pressure=ground['initialization']['surface_surcharge_kpa'];factor=0.
    soil_points=np.repeat(water.soil,np.diff(offset))
    pore=np.maximum(level-positions[:,2],0)*water.gamma*soil_points;wet=water.hydrostatic(active,level)
    old_u=np.zeros(nd);old_s=np.zeros((len(positions),6));old_q=np.zeros((len(system.points),3));rows=[];heads=[];hydraulic_field=None
    for index,frame in enumerate(run['frames']):
        if index:
            stage=case['stages'][index-1];kind=stage['kind']
            if kind=='boundary':factor=stage['target_factor']
            elif kind=='pressure':pressure=stage['target_pressure_kpa']
            elif kind=='hydrostatic_water':
                hydraulic_field=None
                level=stage['target_level_m'];pore=np.maximum(level-positions[:,2],0)*water.gamma*soil_points;wet=water.hydrostatic(active,level)
            elif kind=='steady_seepage':
                level=stage['exterior_water_level_m'];flow=water.steady(stage['flow'],active,mesh['material_ids'],level);field=flow['field']
                hydraulic_field=field
                for e,op in enumerate(q):
                    if not water.soil[e]:
                        pore[offset[e]:offset[e+1]]=0.
                        continue
                    ns=op['nodes'];h=field['relative_head'][ns]
                    pore[offset[e]:offset[e+1]]=np.maximum(field['datum'][ns[0]]-op['position'][:,2]+h[0]+op['shape']@(h-h[0]),0)*water.gamma
                wet=water.from_heads(active,field,level)
                native=next(a for a in run['flow_audits'] if a['stage_index']==index-1)['joint_flow']['original_field']['head_m']
                used=np.unique(ref.cells(mesh)[active & water.soil])
                error=float(np.max(np.abs(np.asarray(native)[used]-field['head'][used])))
                assert error<1e-8;heads.append(error)
            elif kind=='excavate':
                active &= np.array([regions[r] not in stage['bodies'] for r in mesh['region_ids']]);pressure=stage['remaining_surface_pressure_kpa'];wet=water.retain(wet,active)
            else:raise AssertionError(kind)
        frame_ground={**ground,'initialization':{**ground['initialization'],'water_level_m':level}}
        load,_=physical(mesh,frame_ground,active,pressure,True,joint_system=system)
        load+=factor*ref.conditions_load(mesh,case['boundary'],active)+wet['load']
        matrix=np.zeros((nd,nd));seed=np.zeros(nd)
        for e in np.flatnonzero(active):
            op=ref.volume(mesh,e,independent=True);m=mat[mesh['material_ids'][e]];d=ref.elastic_tensor(m['young_kpa'],m['poisson'])
            ns=op['nodes'];ids=(3*ns[:,None]+np.arange(3)).ravel()
            matrix[np.ix_(ids,ids)]+=np.einsum('qai,ab,qbj,q->ij',op['b'],d,op['b'],op['weight'],optimize=True)
            effective=np.asarray(ground['initialization']['geostatic']['seed']['states'][regions[mesh['region_ids'][e]]]['effective_stress_compression_positive_kpa'])
            total=np.tile(effective,(len(op['weight']),1))
            if water.soil[e]:
                if hydraulic_field is None:psi=level-op['position'][:,2]
                else:
                    h=hydraulic_field['relative_head'][ns]
                    psi=hydraulic_field['datum'][ns[0]]-op['position'][:,2]+h[0]+op['shape']@(h-h[0])
                total[:,:3]+=(np.maximum(psi,0)*water.gamma)[:,None]
            seed[ids]-=np.einsum('qai,qa,q->i',op['b'],total,op['weight'])
        live=np.array([active[p['owner_a']] and active[p['owner_b']] for p in system.points])
        nn=np.einsum('qa,qb->qab',system.normal,system.normal)
        stiffness=system.kn[:,None,None]*nn+system.kt[:,None,None]*(np.eye(3)-nn)
        matrix+=np.einsum('qag,qab,qbh,q->gh',system.c,stiffness,system.c,system.weight*live,optimize=True)
        # Input-derived constraints and affine offsets, including frozen orphan
        # displacements after DEATH (from the previous independent solution).
        dummy=dict(active_elements=active.tolist(),boundary_factor=factor,displacement_m=np.zeros_like(xyz).tolist())
        basis,_=ref.staged_constraint_space(mesh,case,dummy);nf=sum(b.shape[1] for b in basis.values());transform=np.zeros((nd,nf));start=0
        for n,b in basis.items():transform[3*n:3*n+3,start:start+b.shape[1]]=b;start+=b.shape[1]
        fixed=np.zeros((len(xyz),3));used=np.unique(ref.cells(mesh)[active]);unused=np.setdiff1d(np.arange(len(xyz)),used);fixed[unused]=old_u.reshape(-1,3)[unused]
        equations={}
        for condition in case['boundary']['constraints']:
            for n in mesh['node_sets'][condition['node_set']]:
                if n not in used:continue
                a,b=equations.setdefault(n,([],[]));a.append(condition['direction']);b.append(factor*(condition['value_m']+np.asarray(condition.get('gradient',[0.,0.,0.]))@(xyz[n]-condition.get('reference_point_m',[0.,0.,0.]))))
        for n,(a,b) in equations.items():fixed[n]=np.linalg.lstsq(a,b,rcond=None)[0]
        fixed=fixed.ravel();rhs=transform.T@(load.ravel()-seed-matrix@fixed);kr=transform.T@matrix@transform
        free=np.linalg.solve(kr,rhs);u=fixed+transform@free
        linear_error=float(np.linalg.norm(kr@free-rhs)/max(np.linalg.norm(rhs),1.));assert linear_error<1e-10
        jump=np.einsum('qag,g->qa',system.c,u);traction=np.einsum('qab,qb->qa',stiffness,jump)
        gn=np.einsum('qa,qa->q',jump,system.normal);tau=traction-np.einsum('qa,qa->q',traction,system.normal)[:,None]*system.normal
        assert np.all(gn[live]<0) and np.all(np.linalg.norm(tau[live],axis=1)<system.mu[live]*system.kn[live]*(-gn[live]))
        old_q[live]=traction[live]
        for e in np.flatnonzero(active):
            op=q[e];m=mat[mesh['material_ids'][e]];d=ref.elastic_tensor(m['young_kpa'],m['poisson']);ids=(3*op['nodes'][:,None]+np.arange(3)).ravel()
            seed_s=np.asarray(ground['initialization']['geostatic']['seed']['states'][regions[mesh['region_ids'][e]]]['effective_stress_compression_positive_kpa'])
            old_s[offset[e]:offset[e+1]]=seed_s-np.einsum('ab,qbi,i->qa',d,op['b'],u[ids])
        state=frame['state'];actual_u=np.asarray(state['displacement_m']).ravel();actual_s=np.asarray([s['effective_stress_compression_positive_kpa'] for s in state['material_states']]);actual_q=np.asarray([s['traction_kpa'] for s in state['joints']['points']])
        ue=float(np.linalg.norm(actual_u-u)/max(np.linalg.norm(u),1e-30));se=float(np.linalg.norm(actual_s-old_s)/max(np.linalg.norm(old_s),1.));qe=float(np.linalg.norm(actual_q-old_q)/max(np.linalg.norm(old_q),1.))
        assert max(ue,se,qe)<1e-7,(frame['name'],ue,se,qe)
        assert all(s['equivalent_plastic_strain']==0 for s in state['material_states'])
        # Direct reduced equilibrium from independently generated fields.
        residual=transform.T@(matrix@actual_u+seed-load.ravel());reassembled=float(np.linalg.norm(residual)/max(np.linalg.norm(load),1.));assert reassembled<1e-10
        assert float(np.max(np.abs(load.ravel()-state['applied_nodal_load_kn'])))<1e-10
        rows.append(dict(frame=frame['name'],displacement_relative_l2=ue,all_gp_stress_relative_l2=se,all_gp_contact_relative_l2=qe,linear_residual=linear_error,reassembled_residual=reassembled))
        old_u=u
    return dict(passed=True,scope='submerged closed stick; elastic MC soil/structure; hydrostatic, steady Darcy and one-side removal',frames=rows,head_max_abs_error_m=max(heads,default=0.))


def main(solver,output,mode='continuous'):
    output.mkdir(parents=True,exist_ok=False)
    _,_,mesh,_=fixture();ground=json.loads((ROOT/'examples/ground-joint-crack-tip.json').read_text());ground['initialization']['water_level_m']=2.
    ground['initialization']['geostatic']['hydrostatic_face_sets']=['exterior']
    case=json.loads((ROOT/'examples/stages-joint-crack-tip.json').read_text());case['joint_system']['joints'][0]['water']={'kind':mode}
    if mode=='soil':
        case['joint_system']['structure_bodies']=['UpperClay']
        case['joint_system']['joints'][0]['water']['soil_side']='a'
    for c in case['boundary']['constraints'][-3:]:c['value_m']=0.
    flow=json.loads((ROOT/'examples/stages-joint-water-flow.json').read_text())['stages'][0]['flow']
    if mode=='soil':flow['heads']=[dict(face_set='bottom',head_m=2.),dict(face_set='Seam:a',head_m=2.1)]
    elif mode=='prescribed':flow['heads']=[dict(face_set='bottom',head_m=2.),dict(face_set='top',head_m=2.)]
    case['stages']=[dict(kind='pressure',name='wet_compression',target_pressure_kpa=20.,initial_increment=.5),dict(kind='steady_seepage',name='darcy',flow=flow,exterior_water_level_m=2.,initial_increment=.5),dict(kind='hydrostatic_water',name='hydrostatic',target_level_m=2.5,initial_increment=.5),dict(kind='excavate',name='wet_remove',bodies=['UpperClay'],remaining_surface_pressure_kpa=0.,initial_increment=.5)]
    case['verification_reject_once_stages']=['wet_compression','wet_remove']
    folder=output/mode;run,inputs=solve(solver,folder,mesh,ground,case)
    result=dict(independent_water_audit=audit(mesh,ground,run),independent_3d=compare(mesh,ground,run),history=restart_and_rejection(solver,folder,run,inputs))
    write(output/'measurements.json',dict(schema='soil-fem-joint-crack-tip-water/1',passed=True,mode=mode,results=result))
    print('Wet crack-tip acceptance passed')


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--mode',choices=['continuous','soil','prescribed'],default='continuous')
    a=p.parse_args();main(a.solver.resolve(),a.output.resolve(),a.mode)
