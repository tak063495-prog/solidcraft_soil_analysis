"""Read the new full-equilibrium evidence in fixed physical/C order with fsum."""
import argparse, itertools, json, math
from pathlib import Path
import numpy as np
from a03aniso_binding import study, fields, load, values, order, compare, raw_observations, norm, metric
from a03aniso_fixture import HEIGHTS,STRESSES,plan
from a03aniso_run import write
SELECTED=('u','stress','traction','pore_pressure','reaction','endpoint_work','numerical_dissipation','elastic_energy','load_work_cumulative')

def assess(run,out):
    out.mkdir(parents=True,exist_ok=False)
    frozen=study.read(run/'plan.json');assert frozen==plan()
    driver=study.read(run/'driver-result.json')
    result=dict(passed=False,a03_complete=False,diagnostic_only=True,new_limited_acceptance=False,
        formal_paths=driver['formal_paths'],formal_states=driver['formal_states'],failures=driver['failures'],
        geometry={},raw={},work_differences={},directional={},height_effects={},implementation_gates={},
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
    # Both horizontal components and directional stress/reaction changes remain visible,
    # including roundoff-size u/traction changes under the unchanged side rollers.
    for h,engine in itertools.product(HEIGHTS,('native','independent')):
        label=f'h{int(h*10)}/{engine}'
        left,right=(load(run,find(h,s,8),engine,4) for s in ('k58','k85'))
        comparison=diff('directional/'+label,left,right)
        component={}
        for k in ('u','traction','reaction','stress'):
            a,b=values(left,k),values(right,k)
            component[k]=[metric(a[...,i],b[...,i]) for i in range(a.shape[-1])]
        result['directional'][label]=dict(k58_minus_k85=comparison,components=component,
            each_input_components={s:{k:[dict(l2=norm(values(load(run,find(h,s,8),engine,4),k)[...,i]),maximum_absolute=float(np.max(np.abs(values(load(run,find(h,s,8),engine,4),k)[...,i])))) for i in range(values(load(run,find(h,s,8),engine,4),k).shape[-1])] for k in ('u','traction')} for s in STRESSES})
    for h,s,engine in itertools.product((.4,.6),STRESSES,('native','independent')):
        label=f'h{int(h*10)}-h5/{s}/{engine}'
        result['height_effects'][label]=diff('height/'+label,load(run,find(h,s,8),engine,4),load(run,find(.5,s,8),engine,4))
    write(out/'difference-vectors.json',vectors);write(out/'measurements.json',result)
    print('formal',result['formal_paths'],result['formal_states'],flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();assess(a.run,a.output)
