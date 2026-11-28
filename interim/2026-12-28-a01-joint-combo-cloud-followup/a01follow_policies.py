"""Exact existing Rust predictor service and complete local partition policies."""
import copy, json, subprocess
import numpy as np

def serial(v):
    if isinstance(v,np.ndarray):return v.tolist()
    if isinstance(v,np.generic):return v.item()
    raise TypeError(type(v).__name__)

class Predictor:
    def __init__(self, executable):
        self.process=subprocess.Popen([str(executable)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True,encoding='utf-8')
    def close(self):
        self.process.stdin.close();self.process.wait(timeout=10)
        if self.process.returncode:raise RuntimeError('predictor bridge failed')
    def __call__(self,m,matrix,rhs,soil_t,joint,target,q,u,born,counts,iteration,projected_load):
        current=copy.deepcopy(joint);current.update(jump=target,traction=q)
        gap=np.einsum('qa,qa->q',target-joint['rest'],m.system.normal)
        current['status']=np.where(gap>0,0,1)
        rows=[[(j,float(c)) for j,c in enumerate(row) if c!=0.] for row in soil_t]
        stress=np.concatenate([s0-np.einsum('qai,i->qa',op['b'],u[ids])@d.T for op,ids,d,s0,*_ in m.ops])
        stress_scale=max(1.,float(np.sqrt((stress[:,:3]**2).sum(axis=1)+2*(stress[:,3:]**2).sum(axis=1)).max()))
        support=None
        if m.frame:
            k=matrix-soil_t.T@m.k@soil_t;n=m.spaces[not born][0].shape[1]
            cond=k[:n,:n]-k[:n,n:]@np.linalg.pinv(k[n:,n:])@k[n:,:n] if len(k)>n else k
            t=soil_t[:,:n];support=t@cond@t.T
        request=dict(mesh=m.mesh,configuration=m.case['joint_system'],matrix=matrix,rhs=rhs,
            rows=rows,old=joint,current=current,born=born,generations=m._a01follow_generations,displacement=u.reshape(-1,3),counts=counts,
            iteration=iteration,projected_load=projected_load,support_k=support,
            pressure_tolerance=max(m.case['settings']['relative_force_tolerance']*stress_scale,1e-10))
        self.process.stdin.write(json.dumps(request,default=serial,allow_nan=False)+'\n');self.process.stdin.flush()
        line=self.process.stdout.readline()
        if not line:raise RuntimeError('predictor service EOF')
        result=json.loads(line)
        if 'error' in result:raise ArithmeticError(result['error'])
        return np.asarray(result['correction']),result

def required(m,joint,target,u,maximum):
    contact=np.maximum(1,np.ceil(np.linalg.norm(target-joint['jump'],axis=1)/maximum).astype(int))
    bulk=[]
    # All inputs use an elastic bulk; these counts still take part in trust/restart policy.
    for op,ids,*_ in m.ops:
        de=np.einsum('qai,i->qa',op['b'],u[ids]-m._a01follow_old_u[ids])
        bulk.extend(np.maximum(1,np.ceil(np.abs(de).sum(axis=1)/m.case['settings']['material_maximum_substep']).astype(int)))
    return np.r_[np.asarray(bulk,dtype=int),contact]
