"""Compare only actual shared coordinates; quadrature comparisons use area integrals."""
import argparse
from decimal import Decimal, localcontext
from math import fsum
from pathlib import Path
import c01_cause_diagnostics as diag


def load(folder, old):
    import numpy as np
    rows, fields = [], {k:[] for k in diag.FIELDS + ('u',)}
    for row, f in old.states(folder):
        rows.append((row['stage'], row['progress']))
        for k in fields:
            fields[k].append(f[k])
    if len(set(rows)) != len(rows):
        raise ValueError('duplicate actual coordinate')
    return dict(zip(rows, range(len(rows)))), {k:np.stack(v) for k,v in fields.items()}


def select(data, coordinates):
    index, fields = data
    if not set(coordinates) <= set(index):
        raise ValueError('missing actual coordinate; no interpolation permitted')
    return {k:v[[index[c] for c in coordinates]] for k,v in fields.items()}


def comparison(a, b, n):
    coords = [(s,i/n) for s in diag.STAGES for i in range(1,n+1)]
    fa, fb = select(a, coords), select(b, coords)
    if set(fa) != set(fb):
        raise ValueError('field set differs')
    return dict(common_states=len(coords), interpolation_used=False,
                fields={k:diag.metric(fa[k],fb[k]) for k in fa},
                stages={s:{k:diag.metric(fa[k][j*n:(j+1)*n],fb[k][j*n:(j+1)*n]) for k in fa} for j,s in enumerate(diag.STAGES)})


def directions(data, n):
    coords = [(s,i/n) for s in diag.STAGES for i in range(1,n+1)]
    fields = [select(x,coords) for x in data]
    result = {}
    for key in diag.FIELDS:
        a, b = fields[0][key]-fields[1][key], fields[1][key]-fields[2][key]
        with localcontext() as context:
            context.prec = 80
            na = Decimal.from_float(fsum(float(x)*float(x) for x in a.flat)).sqrt()
            nb = Decimal.from_float(fsum(float(x)*float(x) for x in b.flat)).sqrt()
            cosine = float(Decimal.from_float(fsum(float(x)*float(y) for x,y in zip(a.flat,b.flat)))/(na*nb)) if na and nb else None
            order = float((na/nb).ln()/Decimal(2).ln()) if na and nb else None
        result[key] = dict(coarse_absolute_norm=float(na),fine_absolute_norm=float(nb),difference_direction_cosine=cosine,observed_order=order)
    return result


def saved_integrals(folder, system, old):
    import numpy as np
    rows, fields = [], {k:[] for k in ('plastic_dissipation','elastic_energy','plastic_slip_integral','traction_integral','u')}
    for row, f in old.states(folder):
        rows.append((row['stage'],row['progress']))
        for k in ('plastic_dissipation','elastic_energy'):
            fields[k].append(fsum(float(x)*float(w) for x,w in zip(f[k],system.weight)))
        for k in ('plastic_slip','traction'):
            fields[k+'_integral'].append([fsum(float(x)*float(w) for x,w in zip(f[k][:,i],system.weight)) for i in range(3)])
        fields['u'].append(f['u'])
    return dict(zip(rows,range(len(rows)))), {k:np.asarray(v) for k,v in fields.items()}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--saved-switch',type=Path,required=True)
    p.add_argument('--saved-path',type=Path,required=True)
    p.add_argument('--experiments',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--reference-levels',type=int,nargs=3,default=(128,256,512))
    args = p.parse_args()
    old,_ = diag.dependencies()
    exps = {}
    summaries = {}
    for folder in sorted(args.experiments.iterdir()):
        if not (folder/'summary.json').exists():
            continue
        s = diag.read(folder/'summary.json')
        exps[(s['control'],s['global_steps'],s['contact_order'])] = folder
        summaries[folder.name] = s
    reference, orders, quadrature, reproduction, unavailable = {}, {}, {}, {}, {}
    for n in (32,64):
        candidate = exps.get(('switch',n,4))
        if candidate is not None:
            reproduction[str(n)] = comparison(load(args.saved_switch/f'g{n}/independent',old),load(candidate/'independent',old),n)
    levels = tuple(args.reference_levels)
    if levels[1] != 2*levels[0] or levels[2] != 2*levels[1]:
        p.error('reference levels must double')
    for control, saved in (('switch',args.saved_switch),('path',args.saved_path)):
        folders = [exps.get((control,n,4)) for n in levels]
        paths = [(folder/'independent') if folder is not None else saved/f'g{n}/independent' for n,folder in zip(levels,folders)]
        if not all((path/'manifest.json').exists() for path in paths):
            unavailable[control] = [n for n,path in zip(levels,paths) if not (path/'manifest.json').exists()]
            continue
        a,b,c = [load(path,old) for path in paths]
        mid,fine = levels[1:]
        reference[f'{control}/independent{mid}__{fine}'] = comparison(b,c,mid)
        native = saved/f'g{mid}/native'
        if (native/'manifest.json').exists():
            reference[f'{control}/native{mid}__independent{mid}'] = comparison(load(native,old),b,mid)
            reference[f'{control}/native{mid}__independent{fine}'] = comparison(load(native,old),c,mid)
        orders[control+'/independent'+'__'.join(map(str,levels))] = directions([a,b,c],levels[0])
    if not reference:
        raise ValueError('no complete reference triple for requested levels')
    for n in (32,64):
        for a,b in ((4,5),(5,6),(4,6)):
            pa,pb = exps.get(('switch',n,a)),exps.get(('switch',n,b))
            if pa is None or pb is None:
                continue
            mesh,ground,case = [diag.read(pa/name) for name in ('mesh.json','ground.json','case.json')]
            from types import SimpleNamespace
            import numpy as np
            sa,sb = [SimpleNamespace(weight=np.asarray([p['weight_m2'] for p in diag.read(path/'geometry.json')['points']])) for path in (pa,pb)]
            quadrature[f'switch/N{n}/q{a}__{b}'] = comparison(saved_integrals(pa/'independent',sa,old),saved_integrals(pb/'independent',sb,old),n)
    report = dict(schema='c01-cause-reference-discretization/1',starting_commit=diag.BASE,passed=False,c01_decision='未達',original_acceptance_retained=True,extrapolation_used=False,interpolation_used=False,experiments=summaries,reference_levels=list(levels),unavailable_reference_sequences=unavailable,reference_comparisons=reference,reference_directions=orders,quadrature_integral_comparisons=quadrature,same_rule_saved_reference_reproduction=reproduction)
    args.output.mkdir(parents=True,exist_ok=False)
    diag.write(args.output/'comparison.json',report)


if __name__ == '__main__':
    main()
