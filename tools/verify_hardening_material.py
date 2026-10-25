"""Input-only independent IDH paths, discrete return checks and refinement."""
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from hardening_reference import integrate,audit,invariants
from make_hardening_examples import initial
from verify_checkpoint import invoke


def write(path,value):
    with path.open('x',encoding='utf-8') as stream:json.dump(value,stream,indent=2)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(solver,output):
    root=Path(__file__).resolve().parents[1];output.mkdir(parents=True,exist_ok=False)
    parameters=json.loads((root/'examples/hardening-parameters.json').read_text())
    definitions=[]
    for exponent in (0.,.5,1.):
        par={**parameters,'stiffness_exponent':exponent}
        definitions.append((f'isotropic-m{exponent}',par,100.,[[.0005,.0005,.0005,0.,0.,0.]]*8+[[-.0003,-.0003,-.0003,0.,0.,0.]]*4))
    definitions.append(('shear-reversal',{**parameters,'dilation_ratio':0.},10000.,
        [[.0008,-.0004,-.0004,0.,0.,0.]]*8+[[-.0008,.0004,.0004,0.,0.,0.]]*10))
    mixed=[[.001,.0003,.0002,.0006,-.0003,.0002]]*6+[[0.,0.,0.,-.0008,.0004,-.0005]]*5+[[-.0003,-.0003,-.0003,0.,0.,0.]]*3
    definitions.append(('mixed-corner',parameters,120.,mixed))
    # Keep the valid high-nu unloading inside the declared compression domain.
    # The original excessive expansion is retained below as a negative case,
    # with an independent upper bound proving that positive p is impossible.
    bounded=mixed[:-3]+[[-.00006,-.00006,-.00006,0.,0.,0.]]*3
    definitions.append(('high-poisson-narrow-cap',{**parameters,'poisson':.4,'cap_axis_ratio':.8},120.,bounded))
    results={};returns=set()
    for name,par,pc,increments in definitions:
        folder=output/name;folder.mkdir();levels=[];summaries=[]
        for maximum in (1e-4,5e-5,2.5e-5):
            label=f'{maximum:.1e}';inp=folder/(label+'-input.json');native=folder/(label+'-native.json');reference=folder/(label+'-reference.json')
            spec={'schema':'soil-fem-hardening-path/1','parameters':par,'initial_state':initial([100.,100.,100.,0.,0.,0.],.7,pc),'maximum_substep':maximum,'increments':increments}
            write(inp,spec);invoke([str(solver),'hardening-path',str(inp),str(native)],folder/(label+'.log'))
            actual=json.loads(native.read_text());assert actual['schema']=='soil-fem-hardening-path-result/1'
            state=copy.deepcopy(spec['initial_state']);states=[state];substeps=0
            for de in increments:
                count=max(1,math.ceil(sum(abs(x) for x in de)/maximum));substeps+=count;state=integrate(par,state,de,count);states.append(state)
            assert actual['material_substeps']==substeps
            write(reference,{'states':states,'method':'Independent analytic elimination of p/q and quadratic shear consistency, scalar cap bracket; same prescribed strain increments and subincrement count. No native state generates a reference.'})
            assert len(actual['states'])==len(states)
            errors={}
            for field in ('effective_stress_compression_positive_kpa','plastic_strain_compression_positive',
                          'shear_plastic_strain','cap_plastic_volumetric_strain','cap_pressure_kpa','void_ratio','plastic_work_kpa'):
                a=np.asarray([s[field] for s in actual['states']]);b=np.asarray([s[field] for s in states])
                err=float(np.max(np.abs(a-b)/np.maximum(1.,np.abs(b))));assert err<2e-8,(name,label,field,err);errors[field]=err
            assert np.max(audit(par,actual['states']))<2e-8
            returns.update(s['last_return'] for s in actual['states']);levels.append(states)
            summaries.append({'maximum_substep':maximum,'maximum_normalized_component_errors':errors})
        refinement=[]
        for field in ('effective_stress_compression_positive_kpa','cap_pressure_kpa','plastic_work_kpa'):
            values=[np.asarray([s[field] for s in states]) for states in levels]
            a=float(np.linalg.norm(values[0]-values[1]));b=float(np.linalg.norm(values[1]-values[2]))
            if a>1e-9:assert b<a*.8,(name,field,a,b)
            refinement.append({'field':field,'coarse_medium_l2_difference':a,'medium_fine_l2_difference':b})
        results[name]={'levels':summaries,'refinement':refinement}
        print('Independent IDH material path verified:',name,flush=True)
    assert returns>={'elastic','shear','cap','shear_cap'}
    bad=copy.deepcopy(spec);bad['initial_state']['cap_pressure_kpa']+=1.
    excessive=dict(spec,increments=mixed);prefix=copy.deepcopy(spec['initial_state'])
    for de in mixed[:-3]:prefix=integrate(spec['parameters'],prefix,de,max(1,math.ceil(sum(abs(x) for x in de)/spec['maximum_substep'])))
    par=spec['parameters'];mean,q=invariants(prefix['effective_stress_compression_positive_kpa'])
    young=par['unloading_young_ref_kpa']*(prefix['cap_pressure_kpa']/par['reference_pressure_kpa'])**par['stiffness_exponent']
    bulk=young/(3*(1-2*par['poisson']));shear=young/(2*(1+par['poisson']))
    volume=sum(sum(v[:3]) for v in mixed[-3:]);upper=float(mean+bulk*volume+bulk*par['dilation_ratio']*q/(3*shear))
    assert upper<par['minimum_mean_effective_pressure_kpa']
    certificate={'maximum_possible_final_mean_pressure_kpa':upper,'minimum_declared_pressure_kpa':par['minimum_mean_effective_pressure_kpa'],
        'prefix_state':prefix,'prescribed_volumetric_expansion':volume,
        'derivation':'For purely isotropic expansion, q is nonincreasing. Integral 3G dlambda_s <= q0; K/G is constant and K never decreases because pc hardens. Hence p_end <= p0+K0*delta_eps_v+(K/G)*dilation*q0/3, even allowing cap plasticity. The negative bound proves this prescribed path leaves the compression domain; it is not a collapse/Fs classification.'}
    write(output/'excessive-expansion-certificate.json',certificate)
    for name,invalid in [('corrupt-cap-history',bad),('invalid-pressure',dict(spec,increments=[[-.1,-.1,-.1,0.,0.,0.]])),('invalid-control',dict(spec,maximum_substep=0.,increments=[])),('excessive-expansion',excessive)]:
        inp=output/(name+'.json');out=output/(name+'-result.json');write(inp,invalid);before=sha(inp)
        invoke([str(solver),'hardening-path',str(inp),str(out)],output/(name+'.log'),1)
        assert not out.exists() and sha(inp)==before
    report={'schema':'soil-fem-hardening-material-verification/1','passed':True,'solver_sha256':sha(solver),
        'cases':results,'native_return_types':sorted(returns),
        'excessive_expansion_domain_certificate':certificate,
        'scope':'Invariant DP shear + elliptic cap with separate hardening; distinct from HS/PH/CYSoil and MC. No empirical sand calibration or small-strain cyclic modulus claim.',
        'sha256':{str(p.relative_to(output)):sha(p) for p in sorted(output.rglob('*')) if p.is_file()}}
    write(output/'verification.json',report)
    print('Independent IDH paths and refinement passed:',len(results))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();verify(a.solver.resolve(),a.output_dir.resolve())
