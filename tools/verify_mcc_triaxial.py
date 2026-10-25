"""Independent axisymmetric MCC rate paths and critical-state benchmarks.

Drained paths hold effective sigma_3, undrained material paths hold volume.
This is a constitutive-point test, not a transient pore-pressure FEM analysis.
Reference formulation: Itasca MCC equations 15,18,21,31; explicit differential
integration here is independent of the Rust implicit return equations.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from verify_checkpoint import invoke


PARAMETERS={'critical_state_m':1.2,'lambda':.2,'kappa':.04,'poisson':.3,'reference_pressure_kpa':1.,'specific_volume_ncl_at_reference':3.5}


def reference(ocr,drained,count,groups=2000,axial_strain=.5):
    assert count%groups==0
    m=PARAMETERS['critical_state_m'];lam=PARAMETERS['lambda'];kap=PARAMETERS['kappa'];nu=PARAMETERS['poisson']
    p=100.;q=0.;pc=ocr*p;v=PARAMETERS['specific_volume_ncl_at_reference']-lam*np.log(pc)+kap*np.log(ocr)
    initial={'effective_stress_compression_positive_kpa':[p,p,p,0.,0.,0.],
        'plastic_strain_compression_positive':[0.]*6,'void_ratio':float(v-1),'preconsolidation_pressure_kpa':pc}
    increments=np.zeros((groups,6));states=[];volumetric=0.;maximum_yield=-1.;maximum_lateral_drift=0.;maximum_plastic_yield=0.
    # Log variables preserve positivity. Monotonic axial loading stays on the
    # plastic branch after its first yield event; no yield clipping is used.
    y=np.array([p,q,np.log(pc),np.log(v)]);plastic=ocr==1.;dz=axial_strain/count
    def yf(y):
        p,q,rho,_=y;pc=np.exp(rho)
        return (q*q+m*m*p*(p-pc))/(m*m*pc*pc)
    def rate(y,plastic):
        p,q,rho,eta=y;pc=np.exp(rho);v=np.exp(eta)
        bulk=v*p/kap;shear=3*bulk*(1-2*nu)/(2*(1+nu));d=np.diag([bulk,3*shear])
        gradient=np.array([m*m*(2*p-pc),2*q]);hardening=m*m*p*pc*v/(lam-kap)*gradient[0]
        denominator=gradient@(d@gradient)+hardening;ep=d
        if plastic:
            assert denominator>0
            a=d@gradient;ep=d-np.outer(a,a)/denominator
        lateral=-nu if drained else -.5
        if drained and plastic:
            row=np.array([1.,-1/3])@ep
            lateral=-(row@np.array([1.,2/3]))/(row@np.array([2.,-2/3]))
        strain=np.array([1+2*lateral,2*(1-lateral)/3])
        multiplier=float(gradient@(d@strain)/denominator) if plastic else 0.
        assert multiplier>=-1e-12
        stress=d@(strain-multiplier*gradient)
        return np.array([*stress,v/(lam-kap)*multiplier*gradient[0],-strain[0]]),np.array([lateral,lateral,1.,0.,0.,0.])
    def step(y,h,plastic):
        first,_=rate(y,plastic);derivative,strain=rate(y+first*h/2,plastic)
        return y+h*derivative,h*strain
    for i in range(count):
        elastic,de=step(y,dz,False) if not plastic else (None,None)
        if not plastic and yf(elastic)>0:
            low=0.;high=1.
            for _ in range(45):
                fraction=(low+high)/2;candidate,_=step(y,dz*fraction,False)
                if yf(candidate)>0:high=fraction
                else:low=fraction
            y,de=step(y,dz*low,False);plastic=True
            y,rest=step(y,dz*(1-low),True);de+=rest
        elif not plastic:y=elastic
        else:y,de=step(y,dz,True)
        p,q,rho,eta=y;pc=np.exp(rho);v=np.exp(eta);volumetric+=de[:3].sum()
        assert p>0 and pc>0 and v>1 and np.isfinite(y).all()
        increments[i//(count//groups)]+=de
        maximum_yield=max(maximum_yield,yf(y));maximum_lateral_drift=max(maximum_lateral_drift,abs(p-q/3-100.))
        if plastic:maximum_plastic_yield=max(maximum_plastic_yield,abs(yf(y)))
        if (i+1)%(count//groups)==0:states.append([p,q,pc,v-1,volumetric])
    return {'initial_state':initial,'strain_increments':increments.tolist(),'states':states,
        'maximum_normalized_yield':float(maximum_yield),'maximum_absolute_plastic_yield':float(maximum_plastic_yield),
        'maximum_lateral_stress_drift_kpa':maximum_lateral_drift}


def verify(solver,output,steps):
    output.mkdir(parents=True,exist_ok=False);reports=[]
    with (output/'reference-code.py').open('xb') as f:f.write(Path(__file__).read_bytes())
    for ocr in (1.,4.):
        for drained in (False,True):
            name=('nc' if ocr==1 else 'oc')+('-drained' if drained else '-undrained')
            coarse=reference(ocr,drained,steps);fine=reference(ocr,drained,2*steps)
            for suffix,data in (('coarse',coarse),('fine',fine)):
                with (output/(name+'-'+suffix+'-reference.json')).open('x',encoding='utf-8') as f:json.dump(data,f)
            expected=np.asarray(fine['states']);refinement=np.linalg.norm(np.asarray(coarse['states'])[:,:3]-expected[:,:3])/np.linalg.norm(expected[:,:3])
            assert refinement<1e-6 and fine['maximum_absolute_plastic_yield']<1e-6
            if drained:assert fine['maximum_lateral_stress_drift_kpa']<1e-7
            errors=[]
            for suffix,h in (('coarse',1e-4),('fine',1e-5)):
                path=output/(name+'-'+suffix+'-path.json');result=output/(name+'-'+suffix+'-result.json')
                data={'schema':'soil-fem-mcc-path/1','parameters':PARAMETERS,'initial_state':fine['initial_state'],
                    'strain_increments_compression_positive':fine['strain_increments'],'maximum_substep':h}
                with path.open('x',encoding='utf-8') as f:json.dump(data,f)
                invoke([str(solver),'mcc-path',str(path),str(result)],output/(name+'-'+suffix+'.log'))
                native=json.loads(result.read_text());states=native['states'][1:];actual=[]
                for s in states:
                    stress=np.asarray(s['effective_stress_compression_positive_kpa']);p=stress[:3].mean();dev=stress.copy();dev[:3]-=p
                    q=np.sqrt(1.5*(dev[:3]@dev[:3]+2*(dev[3:]@dev[3:])))
                    actual.append([p,q,s['preconsolidation_pressure_kpa'],s['void_ratio']])
                actual=np.asarray(actual);error=np.linalg.norm(actual[:,:3]-expected[:,:3])/np.linalg.norm(expected[:,:3]);errors.append(float(error))
                assert error<2e-3 and np.max(np.abs(actual[:,3]-expected[:,3]))<1e-6,(name,suffix,error)
                ncl=1+actual[:,3]+PARAMETERS['lambda']*np.log(actual[:,2])-PARAMETERS['kappa']*np.log(actual[:,2]/actual[:,0])
                assert np.max(np.abs(ncl-PARAMETERS['specific_volume_ncl_at_reference']))<2e-4
                if drained:assert np.max(np.abs(actual[:,0]-actual[:,1]/3-100))<.5
                else:assert np.max(np.abs(actual[:,3]-fine['initial_state']['void_ratio']))<1e-10
            assert errors[1]<errors[0] and errors[1]<2e-5,(name,errors)
            m=PARAMETERS['critical_state_m'];exponent=(PARAMETERS['lambda']-PARAMETERS['kappa'])/PARAMETERS['lambda']
            critical_p=100/(1-m/3) if drained else 100*(ocr/2)**exponent
            critical_q=m*critical_p;last=actual[-1]
            assert abs(last[0]/critical_p-1)<.01 and abs(last[1]/critical_q-1)<.01,(name,last,critical_p,critical_q)
            assert abs(last[2]/(2*last[0])-1)<.01
            reports.append({'case':name,'ocr':ocr,'reference_steps':[steps,2*steps],'reference_refinement_error':float(refinement),
                'native_relative_p_q_pc_error_coarse_fine':errors,'final_p_q_pc_e':last.tolist(),
                'theoretical_critical_p_kpa':critical_p,'theoretical_critical_q_kpa':critical_q,
                'maximum_reference_absolute_plastic_yield':fine['maximum_absolute_plastic_yield'],
                'maximum_native_ncl_specific_volume_error':float(np.max(np.abs(ncl-PARAMETERS['specific_volume_ncl_at_reference'])))})
    report={'schema':'soil-fem-mcc-triaxial-verification/1','passed':True,'solver_sha256':hashlib.sha256(solver.read_bytes()).hexdigest(),
        'reference':'independent axisymmetric elastoplastic differential equations; explicit midpoint with elastic-to-plastic event localization; constant effective lateral stress or constant volume; critical-state closed forms',
        'formulation_source':'https://docs.itascacg.com/itasca920/common/models/camclay/doc/modelcamclay.html',
        'limitations':'Material point paths up to axial strain 0.5 test asymptotic constitutive behavior, not admissible strain for small-strain 3D FEM. Undrained material paths do not solve pore pressure, drainage or consolidation time.',
        'cases':reports,'sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in output.glob('*') if p.is_file()}}
    report['reference_code_sha256']=hashlib.sha256((output/'reference-code.py').read_bytes()).hexdigest()
    with (output/'verification.json').open('x',encoding='utf-8') as f:json.dump(report,f,indent=2)
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--steps',type=int,default=20000)
    a=p.parse_args();verify(a.solver.resolve(),a.output.resolve(),a.steps)
