"""Separate native recovery, common-point interpolation and work components."""
import json
import numpy as np
from common import read
from plate6_offset_reference import System
from plate6_reference import SAMPLES
from beam_reference import uniform_load

def fields(folder,coarse):
    m,c=read(folder/'mesh.json'),read(folder/'stages.json');model=System(m,c['beam_system'])
    cm,cc=coarse;coarse_model=System(cm,cc['beam_system'])
    queries=[np.asarray(l)@coarse_model.positions[nodes][:3] for nodes,_,_ in coarse_model.plate_data for l in SAMPLES]
    owners=[]
    for query in queries:
        choices=[]
        for pi,(nodes,_,kernel) in enumerate(model.plate_data):
            p=model.positions[nodes][:3];r=np.linalg.solve((p[1:,:2]-p[0,:2]).T,query[:2]-p[0,:2]);l=np.r_[1-r.sum(),r]
            if l.min()>=-1e-8:choices.append((pi,l))
        assert choices
        owners.append(choices)
    all_rows={}
    for kind,source in [('native',folder/'native/accepted.jsonl'),('independent',folder/'independent.jsonl')]:
        result=[]
        with source.open() as f:
            for line in f:
                row=json.loads(line);si=row['stage_index'];progress=row['progress']
                if progress not in (.25,.5,.75,1.):continue
                s=row['state'] if kind=='native' else row;bs=s['beams'];d=np.c_[bs['node_displacement_m'],bs['node_rotation_rad']]
                u=np.asarray(s['displacement_m']);factor=s['boundary_factor'] if kind=='native' else (progress if si==2 else 1-progress if si==3 else 0.)
                report=row['beam_results'];reports=report.get('plates',[]);reports={v['id']:v for v in reports} if isinstance(reports,list) else reports
                output=dict(stage_index=si,progress=progress,accepted_index=row['accepted_index'])
                if any(bs['active_plates']):
                    source_shear={key:[] for key in ('nodal','other')}
                    recovered=[];direct=[];raw=[];owner_spread=[];raw_fit=[];direct_raw=[];condition=[]
                    for query,choices in zip(queries,owners,strict=True):
                        values=[]
                        for choice_index,(pi,l) in enumerate(choices):
                            nodes,_,kernel=model.plate_data[pi];neutral=kernel.neutral;kin=d[nodes].ravel()-np.asarray(bs['plate_reference_global_kinematics'][pi]);local=neutral.transform@kernel.map@kin
                            x,y=l@neutral.xy[:3];q=neutral.shear*(neutral.operators(l)[2]@local)
                            samples=reports[model.plates[pi]['id']]['samples'];xy=np.array([np.asarray(p['barycentric'])@neutral.xy[:3] for p in samples]);basis=np.vstack([neutral.basis(*p) for p in xy]);samples_q=np.array([p['transverse_shear_force_kn_m'] for p in samples])
                            coeff=np.linalg.lstsq(basis,samples_q.ravel(),rcond=None)[0];fit=neutral.basis(x,y)@coeff
                            values.append(kernel.axes.T@np.r_[q,0.])
                            if choice_index==0:
                                if 'source_probes' in row:
                                    for source_kind in source_shear:
                                        response=np.asarray(row['source_probes'][source_kind+'_structural_response'])
                                        response_local=neutral.transform@kernel.map@response[nodes].ravel()
                                        source_shear[source_kind].extend(kernel.axes.T@np.r_[neutral.shear*(neutral.operators(l)[2]@response_local),0.])
                                recovered.extend(kernel.axes.T@np.r_[fit,0.]);direct.extend(values[-1]);raw.extend(kernel.axes.T@np.r_[neutral.shear*(neutral.raw_shear(x,y)@local),0.])
                                raw_fit.append(float(np.max(np.abs((basis@coeff).reshape(-1,2)-samples_q))))
                                dq=np.array([neutral.shear*(neutral.operators(p['barycentric'])[2]@local) for p in samples])
                                direct_raw.append(float(np.max(np.abs(dq-samples_q))))
                                condition.append(float(np.linalg.cond(basis)))
                        owner_spread.append(float(np.max(np.ptp(values,axis=0))))
                    output.update(shear_recovered_global_kn_m=recovered,shear_direct_global_kn_m=direct,shear_unprojected_global_kn_m=raw,
                                  recovery_fit_maximum_kn_m=max(raw_fit),native_raw_vs_kinematic_maximum_kn_m=max(direct_raw),
                                  fit_condition_maximum=max(condition),owner_spread_maximum_kn_m=max(owner_spread),owner_counts=[len(v) for v in owners],queries_m=[q.tolist() for q in queries])
                    if 'source_probes' in row:
                        output['load_source_probes']=row['source_probes']
                        for source_kind,values in source_shear.items():output[source_kind+'_source_shear_global_kn_m']=values
                kb,rest,fb,effective,*_=model.assemble(bs,factor)
                applied=np.asarray(s['applied_nodal_load_kn']);effective_soil=np.zeros_like(u)
                for n,attached in enumerate(model.attach):
                    if attached is not None:effective_soil[attached]+=effective[n]
                components=dict(soil_and_holding=float((applied.reshape(-1,3)-effective_soil).ravel()@u.ravel()),plate=0.,beam_distributed=0.,nodal=0.)
                for i,(nodes,section,kernel) in enumerate(model.plate_data):
                    if bs['active_plates'][i]:
                        p=model.plates[i];traction=np.array(p.get('constant_traction_kpa',[0.]*3))+factor*np.asarray(p.get('boundary_traction_kpa',[0.]*3));traction[2]-=section['unit_weight_kn_m3']*section['thickness_m']
                        components['plate']+=float(kernel.load(traction)@d[nodes].ravel())
                for i,(nodes,section,length,axes,t,_) in enumerate(model.data):
                    if bs['active_members'][i]:
                        b=model.members[i];line=np.array(b.get('constant_line_load_kn_m',[0.]*3))+factor*np.asarray(b.get('boundary_line_load_kn_m',[0.]*3));line[2]-=section['unit_weight_kn_m3']*section['area_m2']
                        components['beam_distributed']+=float((t.T@uniform_load(length,axes,line))@d[nodes].ravel())
                for p in c['beam_system']['nodal_loads']:
                    node=model.node_ids[p['node_id']];force=np.asarray(p.get('constant_force_kn',[0.]*3))+factor*np.asarray(p.get('boundary_force_kn',[0.]*3));moment=np.asarray(p.get('constant_moment_kn_m',[0.]*3))+factor*np.asarray(p.get('boundary_moment_kn_m',[0.]*3))
                    if model.active_nodes(bs)[node]:components['nodal']+=float(np.r_[force,moment]@d[node])
                direct_work=float((applied.reshape(-1,3)-effective_soil).ravel()@u.ravel()+fb@d.ravel())
                output.update(work_components_kn_m=components,work_sum_kn_m=sum(components.values()),work_recomposition_error_kn_m=sum(components.values())-direct_work,
                              tip_kinematics=d[model.node_ids['tip']].tolist(),beam_attachment_kinematics=d[model.data[0][0][0]].tolist())
                result.append(output)
        all_rows[kind]=result
    return all_rows
