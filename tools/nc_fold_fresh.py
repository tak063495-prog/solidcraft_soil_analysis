"""Adapt a fresh nc_history_trace proportional path to the EXISTING arc map.

Only the original midpoint 1/128 and 1/256 initial increments are supported.
No old checkpoint is edited, re-sealed, or substituted for a new one.
"""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'interim/2026-11-09-nc-pseudo-arclength'))
from native_arc import Native
from probe_nc_midpoint import inputs
from nc_fixed_inputs import fixed_inputs_match
from verify_high_order import prepare
from low_preload_nc_reference import ClosedStickSystem


class FreshNative(Native):
    def __init__(self,folder,binary,out):
        self.out=out;out.mkdir(parents=True,exist_ok=False);self.binary=binary;self.call=0;self.processes=[]
        run=json.loads((folder/'run.json').read_text(encoding='utf-8'))
        self.mesh=json.loads((folder/'mesh.json').read_text(encoding='utf-8'))
        self.ground=run['ground'];self.case=run['case']
        self.initial_increment=self.case['stages'][2]['initial_increment']
        if self.initial_increment not in (1/128,1/256):raise ValueError('only original 1/128 or 1/256 proportional series')
        expected=inputs('midpoint');expected[2]['stages'][2]['initial_increment']=self.initial_increment
        if not fixed_inputs_match((self.mesh,self.ground,self.case),expected):raise ValueError('fresh proportional series changes fixed inputs')
        rows=[json.loads(line) for line in (folder/'accepted.jsonl').read_text(encoding='utf-8').splitlines()]
        if len(rows)<2 or rows[-1]['state']!=run['last_committed']:raise ValueError('fresh accepted history does not end at run state')
        self.cp=run['last_committed'];self.active=self.cp['active_elements'];self.lam0=self.cp['stage_progress']
        if self.cp['stage_index']!=3 or not .94<self.lam0<.96:raise ValueError('fresh series has no accepted point near original fold')
        self.stage_origin=next(f['state'] for f in run['frames'] if f['name']=='unload')
        system=ClosedStickSystem(self.mesh,self.ground,self.active,self.case)
        if system.live.any() or np.any(system.kjoint):raise ValueError('fresh arc requires all excavation contacts inactive')
        self.R=np.zeros((system.xyz.size,len(system.free)))
        for j,g in enumerate(system.free):
            v=np.zeros(system.xyz.size);v[g]=1.;self.R[:,j]=system.lift(v)
        self.q,self.offsets=prepare(self.mesh);self.rb=[];self.w=[];self.pi=[];self.par=[]
        parameters={m['id']:{**m['mcc'],'integration':self.case['settings']['mcc_integration']} for m in self.ground['materials']}
        for ei in np.flatnonzero(self.active):
            q=self.q[ei];gd=(q['nodes'][:,None]*3+np.arange(3)).ravel()
            for qi,weight in enumerate(q['weight']):
                self.rb.append(q['b'][qi]@self.R[gd]);self.w.append(weight);self.pi.append(int(self.offsets[ei])+qi);self.par.append(parameters[self.mesh['material_ids'][ei]])
        self.rb=np.asarray(self.rb);self.w=np.asarray(self.w);self.pi=np.asarray(self.pi)
        self.holding=np.asarray(run['excavation_audits'][0]['initial_holding_load_kn']);self.load0=np.asarray(self.cp['applied_nodal_load_kn']);self.L=-self.R.T@self.holding
        previous=rows[-2]['state']
        if previous['active_elements']!=self.active:raise ValueError('preceding accepted point changes topology')
        du=np.linalg.lstsq(self.R,(np.asarray(self.cp['displacement_m'])-previous['displacement_m']).ravel(),rcond=None)[0]
        _,self.k0,replayed=self.evaluate(du,previous['material_states'],'initial-accepted-map-replay')
        # Record every GP discrepancy. Never rewrite the saved accepted state.
        discrepancies=[]
        for i in self.pi:
            a=replayed[i];b=self.cp['material_states'][i]
            discrepancies.append(dict(point=int(i),stress_max_abs_kpa=float(np.max(np.abs(np.asarray(a['effective_stress_compression_positive_kpa'])-b['effective_stress_compression_positive_kpa']))),
                                      exact=a==b))
        (out/'initial-replay.json').write_text(json.dumps(dict(all_gp_discrepancies=discrepancies,initial_increment=self.initial_increment,
            note='geometric reconstruction of the last accepted increment; saved state is retained exactly'),indent=2)+'\n',encoding='utf-8')
        self.source=dict(mesh=self.mesh,ground=self.ground,case=self.case,seed=self.cp,stage_origin=self.stage_origin,
                         accepted_history_sha256=hashlib.sha256((folder/'accepted.jsonl').read_bytes()).hexdigest(),
                         observer_binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest())


class Adapter:
    kind='native'
    def __init__(self,folder,binary,out):
        self.m=FreshNative(folder,binary,out/'material');m=self.m
        self.old=m.cp['material_states'];self.u=np.asarray(m.cp['displacement_m']);self.lam=m.lam0
        self.t=np.r_[np.linalg.solve(m.k0,m.L),1.];self.t/=np.linalg.norm(self.t)
        self.contacts=m.cp['joints'];self.R=m.R;self.L=m.L;self.source=m.source
    def evaluate(self,du,old,label):return self.m.evaluate(du,old,label)
    def audit(self,*args):return self.m.audit(*args)
    def serial(self,state):return state
    def state(self,state):return state
