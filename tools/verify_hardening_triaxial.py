"""Independent drained IDH triaxial curves and synthetic stiffness identification.

The scalar lateral-stress controller and quadratic/cap-bracket material update
are independent of Rust. Its fine strain path is prescribed to Rust, as in the
MCC triaxial check. This is synthetic parameter recovery, not field calibration.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import numpy as np
from hardening_reference import integrate,invariants
from make_hardening_examples import initial
from verify_checkpoint import invoke


def write(path,value):
    with path.open('x',encoding='utf-8') as stream:json.dump(value,stream,indent=2)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def reference(parameters,cell,count):
    state=initial([cell,cell,cell,0.,0.,0.],.7,1.2*cell);seed=copy.deepcopy(state)
    increments=[];states=[];axial=[];total=0.;maximum_drift=0.
    def step(axial_increment):
        nonlocal state,total,maximum_drift
        old=state;lateral=-parameters['poisson']*axial_increment
        def evaluate(value):return integrate(parameters,old,[value,value,axial_increment,0.,0.,0.],32)
        for iteration in range(30):
            trial=evaluate(lateral);residual=trial['effective_stress_compression_positive_kpa'][0]-cell
            if abs(residual)<1e-10*max(cell,1.):break
            h=2e-8
            slope=(evaluate(lateral+h)['effective_stress_compression_positive_kpa'][0]-evaluate(lateral-h)['effective_stress_compression_positive_kpa'][0])/(2*h)
            correction=-residual/slope;accepted=False
            for power in range(18):
                candidate=lateral+correction/2**power
                try:test=evaluate(candidate)
                except (AssertionError,ValueError,OverflowError):continue
                if abs(test['effective_stress_compression_positive_kpa'][0]-cell)<abs(residual):lateral=candidate;accepted=True;break
            assert accepted,('independent drained controller',cell,iteration,residual)
        else:raise AssertionError('independent drained controller iteration limit')
        state=trial;total+=axial_increment
        maximum_drift=max(maximum_drift,abs(residual))
        increments.append([lateral,lateral,axial_increment,0.,0.,0.]);states.append(state);axial.append(total)
    for _ in range(count):step(.08/count)
    p,q=invariants(state['effective_stress_compression_positive_kpa'])
    euro=parameters['unloading_young_ref_kpa']*(state['cap_pressure_kpa']/parameters['reference_pressure_kpa'])**parameters['stiffness_exponent']
    excursion=.2*q/euro;portion=count//4
    for _ in range(portion):step(-excursion/portion)
    for _ in range(portion):step(excursion/portion)
    return {'initial_state':seed,'states':states,'increments':increments,'axial_strain':axial,
        'loading_increments':count,'unloading_increments':portion,'maximum_lateral_stress_drift_kpa':maximum_drift}


def verify(solver,output,steps):
    assert steps>=100 and steps%4==0
    output.mkdir(parents=True,exist_ok=False);root=Path(__file__).resolve().parents[1]
    par=json.loads((root/'examples/hardening-parameters.json').read_text());records=[]
    for cell in (25.,100.,400.):
        folder=output/f'cell-{cell:g}';folder.mkdir()
        coarse=reference(par,cell,steps);fine=reference(par,cell,2*steps)
        for label,data in [('coarse',coarse),('fine',fine)]:write(folder/(label+'-reference.json'),data)
        target=np.asarray([s['effective_stress_compression_positive_kpa'] for s in fine['states']]);earlier=np.asarray([s['effective_stress_compression_positive_kpa'] for s in coarse['states']])
        ref_error=float(np.linalg.norm(earlier-target[1::2])/np.linalg.norm(target[1::2]));assert ref_error<2e-4,(cell,ref_error)
        errors=[];actual_fine=None
        for label,h in [('coarse',1e-4),('fine',1e-5)]:
            spec={'schema':'soil-fem-hardening-path/1','parameters':par,'initial_state':fine['initial_state'],'maximum_substep':h,'increments':fine['increments']}
            ip=folder/(label+'-input.json');op=folder/(label+'-native.json');write(ip,spec)
            invoke([str(solver),'hardening-path',str(ip),str(op)],folder/(label+'.log'))
            actual=json.loads(op.read_text())['states'][1:];stress=np.asarray([s['effective_stress_compression_positive_kpa'] for s in actual])
            error=float(np.linalg.norm(stress-target)/np.linalg.norm(target));errors.append(error)
            assert error<.002,(cell,label,error)
            if label=='fine':actual_fine=actual
        assert errors[1]<errors[0]*.6 and errors[1]<5e-4,(cell,errors)
        states=actual_fine;stress=np.asarray([s['effective_stress_compression_positive_kpa'] for s in states]);q=stress[:,2]-stress[:,0]
        load=fine['loading_increments'];unload=fine['unloading_increments'];eps=np.asarray(fine['axial_strain'])
        ultimate=par['shear_ratio_ultimate']*(cell+par['shear_pressure_shift_kpa'])/(1-par['shear_ratio_ultimate']/3)
        assert np.all(np.diff(q[:load])>=-1e-8) and q[load-1]>ultimate/2
        eps50=float(np.interp(ultimate/2,np.r_[0.,q[:load]],np.r_[0.,eps[:load]]));e50=ultimate/(2*eps50)
        eur=float((q[load-1]-q[load+unload-1])/(eps[load-1]-eps[load+unload-1]));pc=states[load-1]['cap_pressure_kpa']
        expected=par['unloading_young_ref_kpa']*(pc/par['reference_pressure_kpa'])**par['stiffness_exponent']
        assert abs(eur/expected-1)<1e-5 and e50<eur
        for s in states[load:]:
            assert abs(s['cap_pressure_kpa']-pc)<1e-6
            assert abs(s['shear_plastic_strain']-states[load-1]['shear_plastic_strain'])<1e-8
        record={'confining_stress_kpa':cell,'reference_refinement_relative_stress_error':ref_error,
            'native_coarse_fine_relative_stress_errors':errors,'maximum_fine_lateral_stress_drift_kpa':float(np.max(np.abs(stress[:,0]-cell))),
            'q_half_asymptotic_strength_kpa':ultimate/2,'axial_strain_at_half_strength':eps50,'measured_e50_kpa':e50,
            'measured_unloading_young_kpa':eur,'preload_cap_pressure_kpa':pc}
        records.append(record);print('Drained IDH curve and unloading identification verified:',cell,flush=True)
    x=np.log([r['preload_cap_pressure_kpa']/par['reference_pressure_kpa'] for r in records]);y=np.log([r['measured_unloading_young_kpa'] for r in records])
    exponent,intercept=np.polyfit(x,y,1);ref=float(np.exp(intercept))
    assert abs(exponent-par['stiffness_exponent'])<1e-5 and abs(ref/par['unloading_young_ref_kpa']-1)<1e-5
    report={'schema':'soil-fem-hardening-drained-triaxial/1','passed':True,'solver_sha256':sha(solver),'cases':records,
        'reference_loading_increments':[steps,2*steps],'reference_subincrements_per_trial':32,'maximum_material_point_axial_strain':.08,
        'synthetic_identification':{'stiffness_exponent':float(exponent),'unloading_young_ref_kpa':ref},
        'scope':'Three confining pressures, independent drained stress control; primary E50 measured at half the DP asymptotic strength, elastic unload/reload, synthetic recovery of unloading stiffness and cap-pressure exponent. 8% material-point strain is not a validated 8% small-strain FEM deformation. No measured site data, variable-dilatancy or liquefaction calibration.',
        'sha256':{str(p.relative_to(output)):sha(p) for p in sorted(output.rglob('*')) if p.is_file()}}
    write(output/'verification.json',report)
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--steps',type=int,default=200)
    a=p.parse_args();verify(a.solver.resolve(),a.output_dir.resolve(),a.steps)
