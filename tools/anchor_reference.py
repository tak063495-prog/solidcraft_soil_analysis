"""Independent host inversion and virtual-work assembly for perfect-bond tendons."""
import copy
import numpy as np
import fem_reference as ref


def locate(mesh,point):
    xyz=np.asarray(mesh['nodes']);cells=ref.cells(mesh);candidates=[]
    bodies={r['id']:r['body'] for r in mesh['regions']}
    for ei,cell in enumerate(cells):
        if point.get('body') is not None and bodies[mesh['region_ids'][ei]]!=point['body']:continue
        x=xyz[cell];target=np.asarray(point['position_m']);a=np.column_stack([x[j]-x[0] for j in (1,2,3)])
        affine=np.linalg.solve(a,target-x[0]);found=None
        for start in (affine,np.full(3,.25)):
            r=start.copy()
            for _ in range(50):
                n,g=ref.basis(np.array([r]),quadratic=len(cell)==10);position=n[0]@x;residual=position-target
                if np.linalg.norm(residual)<1e-12:
                    if r.min()>=-1e-10 and r.sum()<=1+1e-10:found=n[0]
                    break
                j=x.T@g[0];step=np.linalg.solve(j,residual);scale=1.
                for _ in range(20):
                    trial=r-scale*step;nt,_=ref.basis(np.array([trial]),quadratic=len(cell)==10)
                    if np.linalg.norm(nt[0]@x-target)<np.linalg.norm(residual):r=trial;break
                    scale/=2
                else:break
            if found is not None:break
        if found is not None:
            weights=np.zeros(len(xyz));weights[cell]=found;candidates.append((ei,weights))
    assert candidates,('no independent host',point)
    assert all(np.max(np.abs(w-candidates[0][1]))<1e-7 for _,w in candidates)
    return candidates[0][1],[ei for ei,_ in candidates]


class System:
    def __init__(self,mesh,cfg):
        self.mesh=mesh;self.cfg=cfg;self.data=[]
        for d in cfg['anchors']:
            points=[d['head']]+d['bonded_path'];positions=np.asarray([p['position_m'] for p in points])
            bindings=[locate(mesh,p) for p in points];weights=np.array([w for w,_ in bindings]);axis=positions[-1]-positions[0];axis/=np.linalg.norm(axis)
            lengths=np.diff(positions,axis=0)@axis
            c=np.einsum('pn,a->pna',np.diff(weights,axis=0),axis).reshape(len(lengths),-1)
            self.data.append({'definition':d,'positions':positions,'weights':weights,'hosts':[h for _,h in bindings],
                              'axis':axis,'lengths':lengths,'k':d['young_kpa']*d['area_m2']/lengths,'c':c})

    def initial_state(self):
        return {'anchors':[{'active':d['definition']['initially_active'],'reference_extensions_m':[0.]*len(d['k']),
            'force_control_kn':None,'last_lock_off_force_kn':0.,'segment_forces_kn':[0.]*len(d['k']),'segment_taut':[True]*len(d['k'])} for d in self.data]}

    def change(self,state,u,stage):
        result=copy.deepcopy(state);install=stage['kind']=='install_anchors'
        for i,d in enumerate(self.data):
            if d['definition']['id'] not in stage['anchors']:continue
            s=result['anchors'][i];assert s['active']!=install;s['active']=install
            if install:
                s['reference_extensions_m']=(d['c']@np.asarray(u).reshape(-1)).tolist();s['segment_forces_kn']=[0.]*len(d['k']);s['segment_taut']=[True]*len(d['k']);s['last_lock_off_force_kn']=0.
        return result

    def evaluate(self,state,u):
        result=copy.deepcopy(state);internal=np.zeros(np.asarray(u).size);tangent=np.zeros((len(internal),len(internal)));reports=[]
        for i,d in enumerate(self.data):
            s=result['anchors'][i]
            if not s['active']:continue
            extension=d['c']@np.asarray(u).reshape(-1);raw=d['k']*(extension-s['reference_extensions_m']);n=np.maximum(raw,0.);derivative=d['k']*(raw>=0)
            if s['force_control_kn'] is not None:n[0]=s['force_control_kn'];derivative[0]=0.
            internal+=d['c'].T@n;tangent+=d['c'].T@np.diag(derivative)@d['c'];s['segment_forces_kn']=n.tolist();s['segment_taut']=(raw>=0).tolist()
            transfer=np.zeros((len(n)+1,3));transfer[:-1]+=n[:,None]*d['axis'];transfer[1:]-=n[:,None]*d['axis']
            reports.append({'id':d['definition']['id'],'point_positions_m':d['positions'],'point_displacement_m':d['weights']@u,'point_soil_force_kn':transfer,
                            'segment_forces_kn':n,'segment_extension_m':extension,'segment_lengths_m':d['lengths'],
                            'free_length_m':d['lengths'][0],'bonded_length_m':d['lengths'][1:].sum(),'elastic_energy_kn_m':float(np.sum(n*n/(2*d['k'])))})
        return result,internal.reshape(-1,3),tangent,reports


def coupled_forces(system,state,internal,external,report=None):
    result,force,_,reports=system.evaluate(state['anchors'],np.asarray(state['displacement_m']))
    for native,expected in zip(state['anchors']['anchors'],result['anchors']):
        assert np.allclose(native['segment_forces_kn'],expected['segment_forces_kn'],atol=1e-7,rtol=1e-7)
    if report is not None:
        assert len(report['anchors'])==len(reports)
        for native,expected in zip(report['anchors'],reports):
            assert native['id']==expected['id']
            for key,value in expected.items():
                if key!='id':assert np.allclose(native[key],value,atol=1e-7,rtol=1e-7),(key,native[key],value)
        resultant=force.sum(axis=0);moment=np.cross(np.asarray(system.mesh['nodes']),force).sum(axis=0)
        scale=sum(np.linalg.norm(r['point_soil_force_kn'],axis=1).sum() for r in reports)
        assert np.allclose(report['resultant_internal_force_kn'],resultant,atol=1e-8+1e-9*scale,rtol=1e-7)
        tolerance=1e-8+10*system.cfg['embedding_tolerance_m']*scale
        assert np.allclose(report['resultant_internal_moment_kn_m'],moment,atol=tolerance,rtol=1e-7)
        assert np.linalg.norm(resultant)<1e-8+1e-9*scale and np.linalg.norm(moment)<tolerance
    return internal+force,external
