"""Read the new full-equilibrium evidence in fixed physical/C order with fsum."""
import argparse, itertools, json, math
from pathlib import Path
import numpy as np
from a03shear_binding import study, fields, load, values, order, compare, raw_observations, norm, metric
from a03shear_fixture import HEIGHTS,STRESSES,plan
from a03shear_run import write
SELECTED=('u','stress','traction','pore_pressure','reaction','endpoint_work','numerical_dissipation','elastic_energy','load_work_cumulative')

def assess(run,out):
    out.mkdir(parents=True,exist_ok=False)
    frozen=study.read(run/'plan.json');assert frozen==plan()
    driver=study.read(run/'driver-result.json')
    result=dict(passed=False,a03_complete=False,diagnostic_only=True,new_limited_acceptance=False,
        formal_paths=driver['formal_paths'],formal_states=driver['formal_states'],failures=driver['failures'],
        geometry={},raw={},work_differences={},directional={},pre_solve_seed_forces={},prepared_imbalance={},scope={},implementation_gates={},
        water_accounting=dict(flux=None,storage=None,birth_water_transfer=None,balance_acceptance=None),
        fixed_order='phase / exact i/4 / physical connectivity / component; C-order; math.fsum; left-right',
        finest_is_exact=False,acceptance_tolerance_changed=False,s1_s5_changed=False)
    if driver['failures']:raise ValueError('incomplete matrix; do not fabricate comparisons')
    vectors={}
    specs=frozen['specs']
    def find(h,s,n,gl=5):return next(t for t in specs if (t['height'],t['stress'],t['n'],t['volume'])==(h,s,n,gl))
    def diff(label,left,right):
        vectors[label]={k:(values(left,k)-values(right,k)).tolist() for k in SELECTED}
        return compare(left,right)
    # Save every accepted residual and iteration count, not only common points.
    for spec in specs:
        folder=run/spec['name']
        for engine in spec['engines']:
            rows=raw_observations(run,spec,engine)
            if engine=='native':
                native=study.read(folder/'native/run.json')
                attempts=native['geostatic']['attempts']+native['attempts']
                assert len(attempts)==len(rows) and all(a['accepted'] for a in attempts)
                for r,a in zip(rows,attempts):
                    assert r['row']['progress']==a['progress_to']
                    r.update(newton_iterations=a['newton_iterations'],linear_solve_count=len(a['linear_solves']),
                             solver_relative_force_residual=a['relative_force_residual'])
            else:
                it=study.read(folder/'independent-result.json')['rows'];assert len(it)==len(rows)
                for r,a in zip(rows,it):r.update(newton_iterations=a['newton_iterations'],direct_linear_solves=a['direct_linear_solves'],method=a['method'])
            result['raw'][spec['name']+'/'+engine]=rows
        if len(spec['engines'])==2:
            left,right=(result['raw'][spec['name']+'/'+e] for e in spec['engines'])
            result['work_differences'][spec['name']]=[dict(row=a['row'],**{k:a[k]-b[k] for k in ('endpoint_work_kn_m','numerical_dissipation_kn_m','elastic_energy_kn_m','load_work_cumulative_kn_m')}) for a,b in zip(left,right)]
            comparison=study.read(folder/'implementation-comparison.json')['fields']
            force_keys=('internal','applied_load','effective_internal','pore_internal','contact_internal','gravity_load','surcharge_load','exterior_water_load','joint_water_load','holding_load')
            result['implementation_gates'][spec['name']]=dict(
                u_stress_traction_below_original_1e7=all(comparison[k]['relative_l2']<1e-7 for k in ('u','stress','traction')),
                force_absolute_below_original_1e10=all(comparison[k]['maximum_absolute']<1e-10 for k in force_keys),
                categorical_mismatches=sum(v['mismatches'] for v in comparison.values() if 'mismatches' in v))
    for h,s in itertools.product(HEIGHTS,STRESSES):
        label=f'h{int(h*10)}-{s}';entry={}
        for engine in ('native','independent'):
            rows=[load(run,find(h,s,n),engine,4) for n in (4,8,16)]
            entry[engine]=dict(n4_to_8=diff(label+'/'+engine+'/4-8',rows[0],rows[1]),
                n8_to_16=diff(label+'/'+engine+'/8-16',rows[1],rows[2]),
                orders={k:order(*(values(r,k) for r in rows)) for k in SELECTED})
        entry['gl5_gl7']=diff(label+'/gl5-gl7',load(run,find(h,s,8),'independent'),load(run,find(h,s,8,7),'independent'))
        entry['implementation']={str(n):diff(label+'/implementation/'+str(n),load(run,find(h,s,n),'native'),load(run,find(h,s,n),'independent')) for n in (4,8,16)}
        result['geometry'][label]=entry
    # Baseline-centered sign reversal, for every N and each implementation.
    for engine,n,axis,percent in itertools.product(('native','independent'),(4,8,16),('xz','yz'),(1,5)):
        baseline=load(run,find(.5,'zero',n),engine,4)
        plus=load(run,find(.5,f'{axis}{percent}p',n),engine,4)
        minus=load(run,find(.5,f'{axis}{percent}m',n),engine,4)
        label=f'{engine}/n{n}/{axis}{percent}'
        diff('shear-plus-minus/'+label,plus,minus)
        diff('shear-plus-zero/'+label,plus,baseline)
        diff('shear-minus-zero/'+label,minus,baseline)
        components={}; raw={}
        for key in SELECTED:
            positive=values(plus,key)-values(baseline,key)
            negative=values(minus,key)-values(baseline,key)
            defect=positive-(-negative)
            raw[key]=dict(positive_minus_zero=positive.tolist(),negative_minus_zero=negative.tolist(),
                          positive_minus_negative_effect=defect.tolist())
            if key in ('u','traction','reaction','stress'):
                components[key]=[dict(positive_effect_l2=norm(positive[...,i]),negative_effect_l2=norm(negative[...,i]),
                    positive_effect_maximum=float(np.max(np.abs(positive[...,i]))),
                    negative_effect_maximum=float(np.max(np.abs(negative[...,i]))),
                    antisymmetry=metric(positive[...,i],-negative[...,i])) for i in range(positive.shape[-1])]
        vectors['baseline-centered-antisymmetry/'+label]=raw
        result['directional'][label]=dict(components=components,
             definition='Delta+=Y(+tau)-Y(0), Delta-=Y(-tau)-Y(0); defect Delta+ - (-Delta-)',
             even_work_effect_is_not_a_failure_of_horizontal_antisymmetry=True)
    for spec in specs:
        folder=run/spec['name']
        result['pre_solve_seed_forces'][spec['name']]=study.read(folder/'pre-solve-seed-forces.json')
        for engine in spec['engines']:
            prepared=list(fields.load_fields(folder/(engine+'-prepared.jsonl')))
            result['prepared_imbalance'][spec['name']+'/'+engine]=[
                dict(row=row,holding_load_kn=f['holding_load'].tolist(),holding_l2_kn=norm(f['holding_load']),
                     resultant_kn=[math.fsum(map(float,f['holding_load'][:,i])) for i in range(3)],
                     internal_minus_physical_kn=(f['internal']-(f['applied_load']-f['holding_load'])).tolist())
                for row,f in prepared]
            from mc_reference import normalized_yield
            ground=study.read(folder/'ground.json'); mesh=study.read(folder/'mesh.json')
            mapping=study.read(folder/(engine+'-layout.json'))
            mats={v['id']:v['mohr_coulomb'] for v in ground['materials']}
            yields=[];pressures=[];frictions=[]
            for row,f in load(run,spec,engine):
                for e,original in enumerate(mapping['element_order']):
                    yields.extend(normalized_yield(mats[mesh['material_ids'][original]],f['stress'][e*27:(e+1)*27]).tolist())
                if f['joint_active'][0]:
                    pressures.extend(f['normal_pressure'].tolist())
                    tau=np.sqrt(np.sum(f['shear_traction']**2,axis=1))
                    frictions.extend((tau/(.3*f['normal_pressure'])).tolist())
            result['scope'][spec['name']+'/'+engine]=dict(maximum_normalized_mc_yield=max(yields),
                minimum_active_normal_pressure_kpa=min(pressures),maximum_active_friction_ratio=max(frictions),
                all_mc_elastic=max(yields)<=0.,all_contact_closed_stick=min(pressures)>0 and max(frictions)<1.)
    write(out/'difference-vectors.json',vectors);write(out/'measurements.json',result)
    print('formal',result['formal_paths'],result['formal_states'],flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();assess(a.run,a.output)
