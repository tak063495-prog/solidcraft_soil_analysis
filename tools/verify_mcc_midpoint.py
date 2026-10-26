"""Compare the production midpoint CLI to independently integrated MCC rates.

These prescribed paths stay on the loading branch. First-yield events,
unloading and global equilibrium have separate Rust closed-form regressions.
The reference integrates all six stresses, log(pc), log(v) and plastic strain
with Radau; it never calls the Rust return or projects the yield surface.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
from scipy.integrate import solve_ivp
from mcc_rate_fem import terms


def coordinates(state):
    return np.r_[state['effective_stress_compression_positive_kpa'],
                 np.log(state['preconsolidation_pressure_kpa']),
                 np.log1p(state['void_ratio']),state['plastic_strain_compression_positive']]


def verify(solver, output):
    output.mkdir(parents=True, exist_ok=False)
    root=Path(__file__).resolve().parents[1]
    watched=[solver,Path(__file__).resolve(),root/'tools/mcc_rate_fem.py']
    digest=lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    binding={str(p):digest(p) for p in watched}
    params=dict(critical_state_m=1.2,lambda_=0.2,kappa=0.04,poisson=0.3,
                reference_pressure_kpa=1.,specific_volume_ncl_at_reference=3.,integration='midpoint')
    params['lambda']=params.pop('lambda_')
    constants=np.array([[1.44,.2,.04,.3]])
    rows=[]
    for label,stress,increment in [
        ('isotropic',[100.,100.,100.,0.,0.,0.],[.004,.004,.004,0.,0.,0.]),
        ('noncoaxial',[130.,85.,85.,8.,-5.,3.],[.002,.001,.003,.004,-.002,.001]),
        ('softening',[80.,20.,20.,0.,0.,0.],[.03,-.015,-.015,.004,-.002,.003]),
    ]:
        stress=np.asarray(stress);de=np.asarray(increment);p=stress[:3].mean()
        dev=stress.copy();dev[:3]-=p
        q2=1.5*(dev[:3]@dev[:3]+2*dev[3:]@dev[3:])
        initial=dict(effective_stress_compression_positive_kpa=stress.tolist(),void_ratio=1.,
                     preconsolidation_pressure_kpa=float(p+q2/(1.44*p)),
                     plastic_strain_compression_positive=[0.]*6)
        old=coordinates(initial)
        path_parameters={**params,'specific_volume_ncl_at_reference':float(2.+.04*np.log(p)+.16*np.log(initial['preconsolidation_pressure_kpa']))}

        def rhs(_time,state):
            yf,d,n,dn,den=terms(state[None],constants)
            multiplier=float(dn[0]@de/den[0])
            if multiplier < -1e-12: raise ValueError('reference path left the loading branch')
            flow=multiplier*n[0]
            return np.r_[d[0]@(de-flow),np.exp(state[7])/.16*sum(flow[:3]),-sum(de[:3]),flow]

        references=[]
        for suffix,rtol,atol in [('coarse',2e-11,2e-12),('fine',2e-12,2e-13)]:
            r=solve_ivp(rhs,(0.,1.),old,method='Radau',rtol=rtol,atol=atol)
            assert r.success,r.message
            references.append(r.y[:,-1])
            (output/f'{label}-{suffix}-reference.json').write_text(
                json.dumps(dict(method='Radau',rtol=rtol,atol=atol,final=r.y[:,-1].tolist(),
                                nfev=r.nfev,no_projection=True),indent=2),encoding='utf-8')
        reference_delta=float(abs(references[1]-references[0]).max())
        assert reference_delta < 1e-8,(label,reference_delta)
        errors=[];partitions=[]
        for count in (128,256,512):
            path=output/f'{label}-{count}.json';result=output/f'{label}-{count}-result.json'
            path.write_text(json.dumps(dict(schema='soil-fem-mcc-path/1',parameters=path_parameters,
                initial_state=initial,strain_increments_compression_positive=[de.tolist()],
                maximum_substep=float(abs(de).sum()/count)),indent=2),encoding='utf-8')
            run=subprocess.run([str(solver),'mcc-path',str(path),str(result)],capture_output=True,text=True)
            (output/f'{label}-{count}.log').write_text(run.stdout+run.stderr,encoding='utf-8')
            assert run.returncode==0,(label,count,run.stderr)
            data=json.loads(result.read_text(encoding='utf-8-sig'))
            assert data['parameters']['integration']=='midpoint'
            actual=coordinates(data['states'][-1]);error=float(abs(actual-references[-1]).max())
            errors.append(error);partitions.append(data['substeps'][0])
            yf=terms(actual[None],constants)[0][0]
            assert abs(yf)<1e-9,(label,count,yf)
        ratios=[errors[i]/errors[i+1] for i in range(2)]
        assert min(ratios)>3.5,(label,errors,partitions)
        assert errors[-1]<2e-5,(label,errors)
        rows.append(dict(path=label,reference_refinement_maximum=reference_delta,
                         maximum_coordinate_errors=errors,partitions=partitions,ratios=ratios))
    assert binding=={str(p):digest(p) for p in watched},'solver/reference changed during comparison'
    report=dict(schema='soil-fem-midpoint-rate-comparison/1',passed=True,paths=rows,
                sha256=binding,no_projection=True,spatial_convergence_accepted=False,p6_complete=False)
    (output/'verification.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--solver',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    verify(args.solver.resolve(),args.output.resolve())
