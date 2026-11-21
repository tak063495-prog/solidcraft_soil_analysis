"""Measure the h=0.4 paths with the same field order, coordinates and work rule."""
import argparse
import itertools
import math
from pathlib import Path

from a03h04_run import CASES, engine, plan, study
from a03shear_binding import common, compare, fields, load, metric, norm, order, raw_observations, values
from mc_reference import normalized_yield
import numpy as np

SELECTED = ('u', 'stress', 'traction', 'pore_pressure', 'reaction', 'endpoint_work',
            'numerical_dissipation', 'elastic_energy', 'load_work_cumulative')
FORCES = ('internal', 'applied_load', 'effective_internal', 'pore_internal', 'contact_internal',
          'gravity_load', 'surcharge_load', 'exterior_water_load', 'joint_water_load', 'holding_load')


def assess(run, out):
    out.mkdir(parents=True, exist_ok=False)
    specs = plan()
    assert study.read(run/'plan.json') == specs
    driver = study.read(run/'driver-result.json')
    assert driver['formal_paths'] == 21 and driver['formal_states'] == 768 and not driver['failures']
    report = dict(passed=False, formal_paths=21, formal_states=768, paths={}, comparisons={},
                  orders={}, directional={}, implementation_gates={})
    common_rows = {}
    with (out/'signed-differences.jsonl').open('x', encoding='utf-8', newline='\n') as vectors:
        import json

        def difference(label, a, b):
            compared = compare(a, b)
            report['comparisons'][label] = compared
            for key in SELECTED:
                record = dict(pair=label, field=key, left_minus_right=(values(a,key)-values(b,key)).tolist())
                print(json.dumps(record, allow_nan=False, separators=(',', ':')), file=vectors)

        for spec in specs:
            s, n, gl = spec['stress'], spec['n'], spec['volume']
            folder = run/spec['name']
            mesh, ground, _ = [study.read(path) for path in study.paths(folder)]
            materials = {m['id']:m['mohr_coulomb'] for m in ground['materials']}
            for e in spec['engines']:
                rows = load(run, spec, e)
                common_rows[s,n,gl,e] = common(rows, 4)
                raw = raw_observations(run, spec, e)
                if e == 'native':
                    solved = study.read(folder/'native/run.json')
                    attempts = solved['geostatic']['attempts'] + solved['attempts']
                    for r, a in zip(raw, attempts, strict=True):
                        assert a['accepted'] and r['row']['progress'] == a['progress_to']
                        r.update(newton_iterations=a['newton_iterations'], linear_solves=len(a['linear_solves']))
                else:
                    for r, a in zip(raw, study.read(folder/'independent-result.json')['rows'], strict=True):
                        r.update(newton_iterations=a['newton_iterations'], direct_linear_solves=a['direct_linear_solves'])
                layout = study.read(folder/(e+'-layout.json'))
                maximum_yield, minimum_pressure, maximum_ratio = -math.inf, math.inf, 0.
                for _, f in rows:
                    for i, original in enumerate(layout['element_order']):
                        mat = materials[mesh['material_ids'][original]]
                        maximum_yield = max(maximum_yield, max(normalized_yield(mat,f['stress'][i*27:(i+1)*27])))
                    if f['joint_active'][0]:
                        assert np.all(f['contact_status'] == 1)
                        pressure = f['normal_pressure']
                        minimum_pressure = min(minimum_pressure, min(pressure))
                        maximum_ratio = max(maximum_ratio, max(np.sqrt(np.sum(f['shear_traction']**2, axis=1))/(.3*pressure)))
                    assert not np.any(f['material_return']) and not np.any(f['plastic_work']) and not np.any(f['plastic_dissipation'])
                prepared = [dict(row=row, holding_load_kn=f['holding_load'].tolist(),
                                 holding_l2_kn=norm(f['holding_load']),
                                 internal_minus_physical_kn=(f['internal']-(f['applied_load']-f['holding_load'])).tolist())
                            for row,f in fields.load_fields(folder/(e+'-prepared.jsonl'))]
                report['paths'][spec['name']+'/'+e] = dict(
                    rows=raw, prepared=prepared, maximum_normalized_mc_yield=float(maximum_yield),
                    minimum_active_pressure_kpa=float(minimum_pressure), maximum_friction_ratio=float(maximum_ratio))
                assert maximum_yield <= 0 and minimum_pressure > 0 and maximum_ratio < 1
                if gl == 7:
                    std = next(p for p in specs if (p['stress'],p['n'],p['volume'])==(s,n,5))
                    difference(s+'/gl5-gl7/all32',load(run,std,e),rows)
            if len(spec['engines']) == 2:
                matched = study.read(folder/'implementation-comparison.json')['fields']
                report['implementation_gates'][spec['name']] = dict(
                    fields=all(matched[k]['relative_l2'] < 1e-7 for k in ('u','stress','traction')),
                    forces=all(matched[k]['maximum_absolute'] < 1e-10 for k in FORCES),
                    categorical_mismatches=sum(v.get('mismatches',0) for v in matched.values()))
        for s in CASES:
            difference(s+'/gl5-gl7/common16',common_rows[s,8,5,'independent'],common_rows[s,8,7,'independent'])
            for n in (4,8,16):
                difference(f'{s}/n{n}/native-independent',common_rows[s,n,5,'native'],common_rows[s,n,5,'independent'])
            for e in ('native','independent'):
                a,b,c = [common_rows[s,n,5,e] for n in (4,8,16)]
                label = s+'/'+e
                difference(label+'/4-8',a,b)
                difference(label+'/8-16',b,c)
                report['orders'][label] = {k:order(*(values(r,k) for r in (a,b,c))) for k in SELECTED}
        for n, e in itertools.product((4,8,16),('native','independent')):
            base, minus, plus = (common_rows[s,n,5,e] for s in CASES)
            label = f'n{n}/{e}'
            for name,a,b in [('plus-minus',plus,minus),('plus-zero',plus,base),('minus-zero',minus,base)]:
                difference(label+'/'+name,a,b)
            components = {}
            for key in SELECTED:
                dp, dm = (values(r,key)-values(base,key) for r in (plus,minus))
                print(json.dumps(dict(pair=label+'/antisymmetry', field=key, defect=(dp+dm).tolist()),
                                 allow_nan=False,separators=(',', ':')),file=vectors)
                if key in ('u','traction','stress','reaction'):
                    components[key] = [dict(plus_l2=norm(dp[...,i]), minus_l2=norm(dm[...,i]),
                                            defect=metric(dp[...,i],-dm[...,i])) for i in range(dp.shape[-1])]
            report['directional'][label] = components
    engine.write(out/'measurements.json',report)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    assess(args.run,args.output)
