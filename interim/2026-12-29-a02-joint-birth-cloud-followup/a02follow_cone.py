"""New ULP inputs; separate native exact gates, summation order and SVD paths."""
import argparse,itertools,json,math,os,platform,subprocess
from pathlib import Path
import numpy as np
from scipy.linalg import lstsq
from a02follow_run import write,sha
from a02follow_binding import study,metric,norm

def ulp(x,k):
    x=np.asarray(x,dtype=float).copy()
    for _ in range(abs(k)):x=np.nextafter(x,np.inf if k>0 else -np.inf)
    return x

def ordered_sum(values):
    result=0.
    for value in values:result=result+float(value)
    return result

def orders(q,n,mu):
    q,n=np.asarray(q),np.asarray(n);n=n/np.linalg.norm(n);records=[]
    for order in itertools.permutations(range(3)):
        qn=ordered_sum(n[i]*q[i] for i in order);shear=q-qn*n
        length=math.sqrt(ordered_sum(shear[i]*shear[i] for i in order));excess=length-mu*(-qn)
        records.append(dict(order=list(order),normal=qn,shear_norm=length,excess=excess,cone=bool(qn<=0 and excess<=0)))
    qn=math.fsum(float(n[i]*q[i]) for i in range(3));shear=q-qn*n
    length=math.sqrt(math.fsum(float(v*v) for v in shear))
    return dict(permutations=records,fixed_fsum=dict(normal=qn,shear_norm=length,excess=length-mu*(-qn),cone=bool(qn<=0 and length<=mu*(-qn))),
                order_changes_gate=len({r['cone'] for r in records})>1,gate_is_observation_only=True)
class Service:
    def __init__(self,exe,raw):
        self.raw=raw;self.p=subprocess.Popen([str(exe)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True,encoding='utf-8',env=dict(os.environ,A02FOLLOW_RAW_PROJECTION=str(raw)))
    def call(self,v):
        if self.raw.exists():self.raw.unlink()
        self.p.stdin.write(json.dumps(v,allow_nan=False)+'\n');self.p.stdin.flush();line=self.p.stdout.readline()
        if not line:raise RuntimeError('cone service EOF')
        r=json.loads(line)
        if self.raw.exists():r['raw_projection']=json.loads(self.raw.read_text())
        return r
    def close(self):
        self.p.stdin.close();self.p.wait(timeout=20);assert self.p.returncode==0

def main(out,build):
    out.mkdir(parents=True,exist_ok=False);mesh,ground,case=study.fixture('transfer_force',4)
    s=study.Model(mesh,ground,case).system
    matrix=np.einsum('qag,q->gqa',s.c,np.sqrt(s.weight)).reshape(3*len(mesh['nodes']),-1)
    services={profile:Service(build/'source/target'/profile/('a02follow_cone_probe'+('.exe' if os.name=='nt' else '')),out/(profile+'-raw.json')) for profile in ('release','debug')}
    material=case['joint_system']['joints'][0]['material'];mu=material['friction_coefficient'];requests=[];rows=[]
    try:
        for name,n,base,axis in [('friction-x',[0.,0.,1.],[30.,0.,-100.],0),('friction-xy',[0.,0.,1.],[18.,24.,-100.],0),('no-tension',[0.,0.,1.],[0.,0.,0.],2),('rotated-friction',[1.,2.,3.],None,0)]:
            if base is None:
                nn=np.array(n)/np.linalg.norm(n);t=np.cross(nn,[0.,0.,1.]);t/=np.linalg.norm(t);base=(-100*nn+30*t).tolist()
            for k in (-8,-4,-2,-1,0,1,2,4,8):
                q=np.asarray(base).copy();q[axis]=ulp(q[axis],k)
                request=dict(kind='local',material=material,normal=n,traction=q.tolist(),target=[1e-5,0.,-1e-5])
                answers={p:v.call(request) for p,v in services.items()}
                rows.append(dict(kind='local',case=name,ulp_steps=k,input_change=metric(q,base),orders=orders(q,n,mu),answers=answers,builds_exact=answers['release']==answers['debug']))
                requests.append(request)
        for name,traction in [('friction-x',[30.,0.,-100.]),('friction-xy',[18.,24.,-100.])]:
            base=np.einsum('qag,a,q->g',s.c,traction,s.weight)
            for perturb in ('force-components','traction-x'):
                for k in (-16,-8,-4,-2,-1,0,1,2,4,8,16):
                    if perturb=='force-components':
                        force=base.copy();mask=force!=0;force[mask]=ulp(force[mask],k)
                    else:
                        q=np.array(traction);q[0]=ulp(q[0],k);force=np.einsum('qag,a,q->g',s.c,q,s.weight)
                    request=dict(kind='projection',mesh=mesh,case=case,residual_kn=force.tolist())
                    answers={p:v.call(request) for p,v in services.items()};ind={}
                    for driver in ('gelsd','gelss'):
                        x,_,rank,sigma=lstsq(matrix,force,cond=1e-12,lapack_driver=driver)
                        q=x.reshape(-1,3)/np.sqrt(s.weight)[:,None];pressure=-np.einsum('qa,qa->q',q,s.normal)
                        excess=np.linalg.norm(q+pressure[:,None]*s.normal,axis=1)-s.mu*pressure
                        ind[driver]=dict(rank=int(rank),traction=q.tolist(),maximum_force_error=float(max(abs(matrix@x-force))),cone=bool(np.all(pressure>=0) and np.all(excess<=0)),maximum_excess=float(max(excess)),orders=[orders(qq,nn,mm) for qq,nn,mm in zip(q,s.normal,s.mu)])
                    rows.append(dict(kind='projection',case=name,perturbation=perturb,ulp_steps=k,input_change=metric(force,base),answers=answers,independent=ind,builds_exact=answers['release']==answers['debug'],linear_path_difference=metric(ind['gelsd']['traction'],ind['gelss']['traction'])))
                    requests.append(request)
    finally:
        for v in services.values():v.close()
    write(out/'requests.json',requests)
    write(out/'cone-results.json',dict(passed=False,a02_complete=False,diagnostic_only=True,rows=rows,
        unchanged_formula_tolerance_clip=True,rejected_work_is_null=True,work_is_local_prescribed_path_not_global_equilibrium=True,
        windows_ci_after_push='未確定',platform=platform.platform(),build_binding=study.read(build/'binding.json'),source_sha256={p.name:sha(p) for p in Path(__file__).parent.glob('a02follow_cone*') if p.is_file()}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--build',type=Path,required=True);a=p.parse_args();main(a.output.resolve(),a.build.resolve())
