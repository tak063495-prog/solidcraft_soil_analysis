"""Reaggregate sealed raw values in a fixed order; never manufacture acceptance."""
import argparse
import math
from pathlib import Path
import numpy as np
from b01follow_common import EDGES, dot, metric, mm, norm, read, rotate, write
from b01follow_plan import FIELDS, SHAPES, rotations

ENERGIES = ('membrane_energy_kn_m','bending_energy_kn_m','transverse_shear_energy_kn_m','drilling_penalty_energy_kn_m')
TENSORS = {'membrane_strain':True,'curvature_per_m':True,'membrane_force_kn_m':False,'bending_moment_kn':False}
VECTORS = ('transverse_shear_strain','transverse_shear_force_kn_m')


def stiffness_back(k,q):
    k=np.asarray(k); out=np.empty_like(k)
    for i in range(len(k)//3):
        for j in range(len(k)//3):
            out[3*i:3*i+3,3*j:3*j+3]=mm(mm(q.T,k[3*i:3*i+3,3*j:3*j+3]),q)
    return out


def canonical_recovery(records,q,name):
    values={k:[] for k in (*TENSORS,*VECTORS,'drilling_mismatch_rad')}; energies=[[] for _ in ENERGIES]
    for record in records:
        recovery=next(f['recovery'] for f in record['fields'] if f['name']==name)
        axes=np.asarray(recovery['local_axes']); frame=mm(axes,q)
        for i,key in enumerate(ENERGIES):energies[i].append(recovery[key])
        for sample in recovery['samples']:
            for key,engineering in TENSORS.items():
                a,b,c=sample[key]; c=c/2 if engineering else c
                local=np.array([[a,c,0.],[c,b,0.],[0.,0.,0.]])
                values[key].extend(mm(mm(frame.T,local),frame).ravel(order='C'))
            for key in VECTORS:
                values[key].extend(mm(np.array([[*sample[key],0.]]),frame).ravel())
            values['drilling_mismatch_rad'].append(sample['drilling_mismatch_rad'])
    values={k:np.array(v) for k,v in values.items()}
    values['energy_components_kn_m']=np.array([math.fsum(v) for v in energies])
    return values


def observables(record,recoveries):
    q=np.array(record['rotation_matrix']); result={}
    for row in record['solves']:
        name=row['name']; value=canonical_recovery(recoveries,q,name)
        for field,key in [('displacement_rotation','displacement'),('reaction','reaction'),('internal_force_moment','internal')]:
            v=rotate(np.array(row[field]).reshape(-1,6),q.T)
            value[key+('_translation_m' if key=='displacement' else '_force_kn')]=v[:,:3].ravel()
            value[key+('_rotation_rad' if key=='displacement' else '_moment_kn_m')]=v[:,3:].ravel()
        value['stiffness_energy_kn_m']=np.array([row['strain_energy_kn_m']])
        value['applied_work_kn_m']=np.array([row['applied_work_kn_m']])
        result[name]=value
    return result


def report(run,output,effects):
    records={p.stem:read(p) for p in sorted((run/'assemblies').glob('*.json'))}
    kernel=read(run/'kernel-responses.json'); byid={r['id']:r for r in kernel}
    recovery=read(run/'recovery-responses.json'); recoveries={r['id']:r for r in recovery}
    obs={}
    for name,record in records.items():
        obs[name]=observables(record,[recoveries[name+'--e'+str(i)] for i in range(8)])
    rows=[]; raw_effects=[]
    for shape in SHAPES:
        base=records[shape+'--identity']; baseobs=obs[base['id']]
        for rotation,q in rotations():
            name=shape+'--'+rotation; record=records[name]; value=obs[name]; spectrum=record['spectrum']
            axes=[]; midpoint=[]; areas=[]
            for e,nodes in enumerate(record['elements']):
                element=byid[name+'--e'+str(e)]; ref=byid[base['id']+'--e'+str(e)]
                axes.append(metric(mm(np.array(element['axes']),q),ref['axes'])['maximum_absolute'])
                areas.append(element['area_m2']); x=np.array(record['points'])[nodes]
                scale=max(norm(x[i]-x[j]) for i,j in EDGES)
                midpoint.extend(norm(x[3+e]-(x[a]+x[b])*.5)/(64*np.finfo(float).eps*scale) for e,(a,b) in enumerate(EDGES))
            covariance={}; distortion={}; interaction={}
            for field in FIELDS:
                covariance[field]={}; distortion[field]={}; interaction[field]={}
                for key,v in value[field].items():
                    ref=baseobs[field][key]; reg=obs['regular--'+rotation][field][key]; reg0=obs['regular--identity'][field][key]
                    delta=v-ref; shape_delta=v-reg; cross=delta-(reg-reg0)
                    covariance[field][key]=metric(v,ref); distortion[field][key]=metric(v,reg)
                    interaction[field][key]=metric(cross,np.zeros_like(cross))
                    raw_effects.append(dict(shape=shape,rotation=rotation,field=field,quantity=key,
                        rotation_delta=delta.tolist(),distortion_delta=shape_delta.tolist(),cross_difference=cross.tolist(),
                        correspondence='fixed topology node/element/sample IDs; shapes have different physical positions'))
            rigid_energies={}
            for i in range(6):
                rigid_energies[str(i)]=canonical_recovery([recoveries[name+'--e'+str(e)] for e in range(8)],q,'physical-rigid-'+str(i))['energy_components_kn_m'].tolist()
            rows.append(dict(id=name,shape=shape,rotation=rotation,area_m2=math.fsum(areas),
                q_orthogonality_max_abs=metric(mm(q.T,q),np.eye(3))['maximum_absolute'],
                axes_covariance_max_abs=max(axes),midpoint_maximum_cutoff_ratio=max(midpoint),
                mixed_dof_stiffness_covariance=metric(stiffness_back(record['stiffness'],q),base['stiffness']),
                covariance=covariance,distortion=distortion,cross_difference=interaction,
                algebraic_zero_count=spectrum['modes_below_inherited_threshold'],
                negative_below_minus_threshold=spectrum['negative_below_minus_threshold'],
                physical_rigid_normalized_residual=spectrum['physical_rigid_normalized_residual'],
                physical_rigid_fraction_outside_nullspace=spectrum['physical_rigid_fraction_outside_nullspace'],
                rigid_recovery_energy_components_kn_m=rigid_energies,
                reduced_mixed_dof_condition_number=record['reduced_mixed_dof_condition_number'],
                condition_times_epsilon=record['condition_times_epsilon'],
                global_solution_observations=[{k:v for k,v in s.items() if k not in
                   ('displacement_rotation','applied_force_moment','internal_force_moment','residual','reaction','target')} for s in record['solves']],
                energy_components={f:value[f]['energy_components_kn_m'].tolist() for f in FIELDS}))
    rejected=[dict(id=r['id'],accepted=r['accepted'],error=r.get('error')) for r in kernel if r['id'].startswith('reject--')]
    summary=dict(passed=False,b01_complete=False,diagnostic_only=True,production_scope_expanded=False,
        curved_element_acceptance=None,new_global_solves=192,patches=32,maximum_elements=8,maximum_dofs=150,
        original_zero_mode_threshold=1e-10,algebraic_zero_counts=sorted(set(r['algebraic_zero_count'] for r in rows)),
        same_input_repeat=read(run/'same-input-repeat.json'),rejected_controls=rejected,
        permitted_element_requests_accepted=sum(r['accepted'] for r in kernel if not r['id'].startswith('reject--')),
        axes_covariance_max_abs=max(r['axes_covariance_max_abs'] for r in rows),
        stiffness_covariance_max_relative=max(r['mixed_dof_stiffness_covariance']['relative_l2'] for r in rows),
        physical_rigid_residual_max=max(max(r['physical_rigid_normalized_residual']) for r in rows),
        physical_rigid_outside_nullspace_max=max(max(r['physical_rigid_fraction_outside_nullspace']) for r in rows),
        eigenvalue_and_condition_units='mixed translation/rotation DOF algebraic observations only',
        physical_force_and_moment_are_separate=True,new_acceptance_thresholds=[],extrapolation=False,cases=rows)
    write(effects,dict(diagnostic_only=True,passed=False,b01_complete=False,records=raw_effects))
    write(output,summary)
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--effects',type=Path,required=True)
    a=p.parse_args();report(a.run.resolve(),a.output.resolve(),a.effects.resolve())
