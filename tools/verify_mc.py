"""CLI MC paths against iterative complementarity and column closed form."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np
import mc_reference as ref


def verify(solver,output):
    output.mkdir(parents=True,exist_ok=False);rng=np.random.default_rng(20261023)
    initial={'effective_stress_compression_positive_kpa':[100.,100.,100.,0.,0.,0.],
        'plastic_strain_compression_positive':[0.]*6,'void_ratio':.7,'equivalent_plastic_strain':0.,'plastic_work_kpa':0.,'last_return':'elastic'}
    base={'young_kpa':10000.,'poisson':.25,'cohesion_kpa':10.,'friction_deg':30.,'dilation_deg':5.,'tension_cutoff_kpa':0.,'flow':'non_associated'}
    increments=rng.uniform(-.006,.006,(100,6));increments[:,:3]-=increments[:,:3].mean(axis=1)[:,None]
    increments[-10:,:3]-=.012;reports=[]
    configurations=[('nonassociated',base),('zero_dilation',{**base,'dilation_deg':0.}),
        ('associated',{**base,'flow':'associated','dilation_deg':30.}),('davis',{**base,'flow':'davis','dilation_deg':0.}),
        ('tresca',{**base,'friction_deg':0.,'dilation_deg':0.,'tension_cutoff_kpa':None})]
    # A pure-MC hydrostatic tension path with psi>0 also checks its true apex,
    # separately from the Rankine cap present in the other tension paths.
    configurations.append(('pure_mc_apex',{**base,'tension_cutoff_kpa':None}))
    configurations.append(('mixed',base))
    configurations.append(('tension_face',{**base,'cohesion_kpa':50.}))
    modes=set()
    for name,p in configurations:
        path=increments
        if name in ('pure_mc_apex','mixed','tension_face'):
            e=p['young_kpa'];nu=p['poisson'];g=e/(2*(1+nu));lam=e*nu/((1+nu)*(1-2*nu));d=np.full((3,3),lam);np.fill_diagonal(d,lam+2*g)
            trial=np.array([-100.,-110.,-120.])
            if name=='tension_face':trial=np.array([100.,20.,-10.])
            if name=='mixed':
                _,phi,psi=ref.strength(p);m=np.array([1-np.sin(psi),0.,-1-np.sin(psi)])
                trial=np.array([2*p['cohesion_kpa']/np.tan(phi),p['cohesion_kpa']/np.tan(phi),0.])+d@m*.01-d[:,2]*.01
            path=np.array([np.r_[np.linalg.solve(d,trial-100),[0.,0.,0.]]])
        source=output/f'{name}.json';result=output/f'{name}-result.json'
        source.write_text(json.dumps({'schema':'soil-fem-mc-path/1','parameters':p,'initial_state':initial,
            'strain_increments_compression_positive':path.tolist(),'maximum_substep':.001},indent=2))
        subprocess.run([str(solver),'mc-path',str(source),str(result)],check=True)
        run=json.loads(result.read_text());state=initial;errors={k:0. for k in ('stress','plastic_strain','void_ratio','equivalent_plastic_strain','plastic_work')}
        for de,actual in zip(path,run['states'][1:]):
            state=ref.advance(p,state,de,.001);modes.add(actual['last_return'])
            for key,field in [('stress','effective_stress_compression_positive_kpa'),('plastic_strain','plastic_strain_compression_positive'),('void_ratio','void_ratio'),('equivalent_plastic_strain','equivalent_plastic_strain'),('plastic_work','plastic_work_kpa')]:
                error=float(np.linalg.norm(np.asarray(actual[field])-state[field]));errors[key]=max(error,errors[key])
            assert errors['stress']<1e-5 and errors['plastic_strain']<1e-7 and errors['void_ratio']<1e-11 and errors['equivalent_plastic_strain']<1e-7 and errors['plastic_work']<1e-5,(name,errors)
        assert any(s['equivalent_plastic_strain']>1e-6 for s in run['states']),('path never plastic',name)
        reports.append({'model':name,'increments':len(path),'maximum_absolute_errors':errors})
    assert {'elastic','shear_face','compression_edge','extension_edge','shear_apex','tension_face','tension_edge','tension_apex','mixed'}<=modes,modes
    report={'schema':'soil-fem-independent-mc-paths/1','method':'NumPy principal stresses, projected Gauss-Seidel complementarity over all MC/Rankine planes; Rust uses active-set enumeration',
        'paths':reports,'observed_return_modes':sorted(modes),'solver_sha256':hashlib.sha256(solver.read_bytes()).hexdigest(),
        'files_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in output.glob('*.json')}}
    for p in (Path(__file__),Path(ref.__file__)):report['files_sha256'][str(p.resolve())]=hashlib.sha256(p.read_bytes()).hexdigest()
    (output/'verification.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,default=Path('target/release/soil-fem-3d.exe'));p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();verify(a.solver.resolve(),a.output_dir.resolve())
