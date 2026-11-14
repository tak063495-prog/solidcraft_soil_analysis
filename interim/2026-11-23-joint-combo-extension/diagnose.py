"""Saved-state/local audits, explicitly separate from fresh 3D solves."""
import argparse
import itertools
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from fixtures import read,write
from reference import Model,contact,contact_work
from fields import load_fields


def diagnose(folder,steps):
    m=Model(*(read(folder/(key+'.json')) for key in ('mesh','ground','case')));s=m.system
    native=read(folder/'native/run.json');assert native['complete']
    audits={a['stage']:a['prepared_state'] for a in native['transition_audits']}
    phases=['geostatic']+[x['name'] for x in m.case['stages']]
    reference=read(folder/'independent-result.json');assert reference['complete']
    old=None;maxima={key:0. for key in ('traction','plastic_slip_increment','accumulated_increment','plastic_dissipation_increment','numerical_dissipation_increment','endpoint_work_increment')}
    mismatch=0;native_counts={};reference_counts={};rows=[];force_max=load_max=balance_max=0.;residual_max=0.;reset_checks=0
    for i,(a,b) in enumerate(itertools.zip_longest((folder/'native/accepted.jsonl').open(),(folder/'native/constitutive.jsonl').open())):
        assert a is not None and b is not None
        state=__import__('json').loads(a);diag=__import__('json').loads(b)
        phase=i//steps;name=phases[phase];progress=(i%steps+1)/steps
        if i%steps==0 and name in audits:
            old=audits[name]
            if name in ('birth','rebirth'):
                for p in old['joints']['points']:
                    assert p['plastic_slip_m']==[0.,0.,0.]
                    assert all(p[key]==0. for key in ('accumulated_plastic_slip_m','plastic_dissipation_kpa_m','numerical_dissipation_kpa_m','endpoint_work_kpa_m'))
                    reset_checks+=1
        points=state['joints']['points'];active=all(p['active'] for p in state['joints']['lifecycle']);u=np.asarray(state['displacement_m']).ravel()
        now=lambda key:np.array([p[key] for p in points])
        counts=np.array(diag['partitions'][len(state['material_states']):]);rc=np.array(reference['rows'][i]['contact_partitions'])
        local_mismatch=int(np.sum(counts!=rc));mismatch+=local_mismatch
        if old is not None and active:
            prior=old['joints']['points'];before=lambda key:np.array([p[key] for p in prior])
            rest=np.array([p.get('birth',{}).get('elastic_offset_m',[0.,0.,0.]) for p in prior])
            for k in counts:native_counts[int(k)]=native_counts.get(int(k),0)+1
            for k in rc:reference_counts[int(k)]=reference_counts.get(int(k),0)+1
            for count in np.unique(counts):
                mask=counts==count;sub=SimpleNamespace(**{key:getattr(s,key)[mask] for key in ('normal','kn','kt','mu')})
                start=before('jump_m')[mask]-rest[mask];target=now('jump_m')[mask]-rest[mask];oq=before('traction_kpa')[mask]
                q,dp,da=contact(sub,start,oq,target,int(count));pd,num,w=contact_work(sub,start,oq,target,int(count))
                observations=[('traction',q[:,0],now('traction_kpa')[mask]),
                    ('plastic_slip_increment',dp[:,0],now('plastic_slip_m')[mask]-before('plastic_slip_m')[mask]),
                    ('accumulated_increment',da[:,0],now('accumulated_plastic_slip_m')[mask]-before('accumulated_plastic_slip_m')[mask]),
                    ('plastic_dissipation_increment',pd,now('plastic_dissipation_kpa_m')[mask]-before('plastic_dissipation_kpa_m')[mask]),
                    ('numerical_dissipation_increment',num,now('numerical_dissipation_kpa_m')[mask]-before('numerical_dissipation_kpa_m')[mask]),
                    ('endpoint_work_increment',w,now('endpoint_work_kpa_m')[mask]-before('endpoint_work_kpa_m')[mask])]
                for key,x,y in observations:maxima[key]=max(maxima[key],float(np.max(np.abs(x-y),initial=0.)))
        # Actual saved stress force audit; this is not input-only new solving.
        force=np.zeros(m.nd);offset=0
        for op,ids,*_ in m.ops:
            n=len(op['weight']);sigma=np.array([p['effective_stress_compression_positive_kpa'] for p in state['material_states'][offset:offset+n]])
            force[ids]-=np.einsum('qai,qa,q->i',op['b'],sigma,op['weight']);offset+=n
        if active:force+=np.einsum('qag,qa,q->g',s.c,now('traction_kpa'),s.weight)
        attached=np.zeros(m.nd)
        if m.frame:
            br=m.frame.solve_interface(state['beams'],u.reshape(-1,3),state['boundary_factor'])
            for node,v in br['interface_residual'].items():force[3*node:3*node+3]+=v
            for node,v in br['attached_load'].items():attached[3*node:3*node+3]+=v;force[3*node:3*node+3]+=v
        if m.tendon:force+=m.tendon.evaluate(state['anchors'],u.reshape(-1,3))[1].ravel()
        actual_internal=np.asarray(diag['internal']);fe=float(np.max(np.abs(force-actual_internal),initial=0.));force_max=max(force_max,fe)
        holding=state.get('geostatic_release_load_kn') if phase==0 else state.get('transition_release_load_kn')
        if not holding:holding=np.zeros(m.nd)
        load=m.base+state['boundary_factor']*m.boundary+attached+holding
        le=float(np.max(np.abs(load-np.asarray(state['applied_nodal_load_kn'])),initial=0.));load_max=max(load_max,le)
        t,_=m.spaces[not active];res=float(np.linalg.norm(t.T@(force-load))/max(np.linalg.norm(load),1.));residual_max=max(residual_max,res)
        q=now('traction_kpa');qn=np.einsum('qa,qa->q',q,s.normal);tau=q-qn[:,None]*s.normal
        energy=qn**2/(2*s.kn)+(tau*tau).sum(axis=1)/(2*s.kt)
        initial=np.array([p.get('birth',{}).get('initial_energy_kpa_m',0.) for p in points])
        energy_error=float(np.max(np.abs(initial+now('endpoint_work_kpa_m')-energy-now('plastic_dissipation_kpa_m')-now('numerical_dissipation_kpa_m')),initial=0.));balance_max=max(balance_max,energy_error)
        rows.append(dict(stage=name,progress=progress,partition_mismatches=local_mismatch,native_partitions=counts.tolist(),
            independent_fresh_partitions=rc.tolist(),all_node_internal_reassembly_maximum_kn=fe,all_node_load_reassembly_maximum_kn=le,
            reassembled_force_residual=res,contact_work_balance_maximum_kpa_m=energy_error))
        old=state
    assert len(rows)==len(phases)*steps
    # Existing consistency gates; no new full-history error threshold.
    assert force_max<1e-10 and load_max<1e-10 and balance_max<1e-10
    assert residual_max<=m.case['settings']['relative_force_tolerance']
    assert max(maxima.values())<1e-8
    return dict(saved_state_audit_only=True,counted_as_fresh_3d_solves=False,accepted_states=len(rows),
        every_gp_local_replay_with_observed_native_partitions=True,local_replay_maximum_errors=maxima,
        contact_partition_mismatches=mismatch,native_partition_histogram=native_counts,independent_partition_histogram=reference_counts,
        all_node_internal_reassembly_maximum_kn=force_max,all_node_load_reassembly_maximum_kn=load_max,
        maximum_reassembled_force_residual=residual_max,maximum_contact_work_balance_error_kpa_m=balance_max,
        exact_new_generation_zero_history_points=reset_checks,rows=rows,
        conclusion='Native-vs-input-only differences include independently selected local partitions. Same-native-jump/same-native-partition replay separately checks law/history/work and nodal assembly. No production defect or full-field acceptance is inferred from small endpoint differences.')


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--steps',type=int,default=8);a=p.parse_args()
    result={}
    for folder in sorted(a.output.iterdir()):
        if (folder/'comparison.json').exists():result[folder.name]=diagnose(folder,a.steps)
    write(a.output/'cause-diagnostics.json',result)


if __name__=='__main__':main()
