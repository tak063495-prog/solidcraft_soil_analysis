"""Bounded local diagnostics. Run directly; writes only to --out."""
import argparse
import json
import platform
import time
from pathlib import Path
import numpy as np
import scipy
from mc import Material, implicit, explicit, yield_value, wet_update
from column import solve


def paths():
    zero = np.zeros((3, 3))
    triax = [np.diag([.002*a, .002*a, -.01*a]) for a in np.linspace(0,1,9)]
    shear = []
    reload = []
    for a in np.linspace(0,1,13):
        e = np.zeros((3,3))
        e[0,2] = e[2,0] = .005*a
        shear.append(e)
    for a in np.concatenate((np.linspace(0,1,7),np.linspace(1,-.6,9)[1:],np.linspace(-.6,1.2,10)[1:])):
        e = np.diag([.0007*a, -.0002*a, -.001*a])
        e[0,2] = e[2,0] = .004*a
        reload.append(e)
    return dict(triaxial_compression=triax, simple_shear=shear, unload_reload=reload)


def integrate(targets, n, method, initial=None):
    s = np.diag([-60., -100., -140.]) if initial is None else initial.copy()
    ep = np.zeros((3, 3))
    stresses, plastics, works = [], [], []
    branches = {}
    previous = targets[0]
    max_drift = 0.
    for e in targets[1:]:
        de = e-previous
        if method == 'explicit':
            r = explicit(s, de, n)
            s, ep = r['stress'], ep+r['dep']
            works.append(r['work'])
            max_drift = max(max_drift, r['max_yield_drift'])
        else:
            for _ in range(n):
                r = implicit(s, de/n)
                s, ep = r['stress'], ep+r['dep']
                works.append(r['work'])
                branches[r['branch']] = branches.get(r['branch'],0)+1
        stresses.append(s.copy())
        plastics.append(ep.copy())
        previous = e
    return dict(stress=np.array(stresses), plastic=np.array(plastics),
                minimum_plastic_work=min(works), work=sum(works),
                max_yield_drift=max_drift, branches=branches)


def tangent_checks():
    old = -100*np.eye(3)
    # Distinct trial eigenvalues, chosen to exercise classical branch interiors.
    cases = dict(elastic=np.diag([.0001,-.0001,0.]),
                 face=np.array([[.0003,.003,0],[.003,-.0002,0],[0,0,-.0008]]),
                 edge=np.diag([.002,.0021,-.01]),
                 apex=np.diag([.01,.008,.009]))
    d = np.array([[.3,.2,-.1],[.2,-.4,.17],[-.1,.17,.2]])
    rows = []
    for expected, de in cases.items():
        r = implicit(old,de,tangent=True)
        for h in (1e-7,1e-8,1e-9):
            plus, minus = implicit(old,de+h*d), implicit(old,de-h*d)
            fd = (plus['stress']-minus['stress'])/(2*h)
            exact = r['tangent'](d)
            rows.append(dict(expected=expected, branch=r['branch'], h=h,
                             same_branch=plus['branch']==minus['branch']==r['branch'],
                             absolute_error_kpa=float(np.linalg.norm(fd-exact)),
                             error_divided_by_E=float(np.linalg.norm(fd-exact)/Material().E),
                             scaled_error=float(np.linalg.norm(fd-exact)/max(np.linalg.norm(exact),1.))))
    return rows


def main(out, levels):
    start = time.perf_counter()
    result = dict(schema='py-c03-diagnostics/1', accepted=None,
                  tolerance='未定義', C03='未達', native_comparison_executed=False,
                  versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
                  paths={}, tangent=tangent_checks())
    point_histories={}
    for name, targets in paths().items():
        rows=[]
        previous = None
        previous_solutions = None
        previous_adjacent = None
        for n in levels:
            initial = -100*np.eye(3) if name=='triaxial_compression' else np.diag([-60.,-100.,-140.])
            a, b = integrate(targets,n,'implicit',initial), integrate(targets,n,'explicit',initial)
            err = float(np.linalg.norm(a['stress']-b['stress']))
            adjacent = (None if previous_solutions is None else
                        [float(np.linalg.norm(a['stress']-previous_solutions[0])),
                         float(np.linalg.norm(b['stress']-previous_solutions[1]))])
            own_orders = ([None,None] if previous_adjacent is None else
                          [None if min(previous_adjacent[i],adjacent[i])<1e-9 else
                           float(np.log2(previous_adjacent[i]/adjacent[i])) for i in range(2)])
            rows.append(dict(substeps_per_increment=n, stress_difference_kpa=err,
                             initial_effective_stress_kpa=initial.tolist(),
                             adjacent_implicit_explicit_stress_differences_kpa=adjacent,
                             adjacent_observed_orders_implicit_explicit=own_orders,
                             plastic_strain_difference=float(np.linalg.norm(a['plastic']-b['plastic'])),
                             observed_order=(None if previous is None or min(err,previous[1])<1e-9 else
                                             float(np.log(previous[1]/err)/np.log(n/previous[0]))),
                             implicit_min_plastic_work_kpa=a['minimum_plastic_work'],
                             explicit_min_plastic_work_kpa=b['minimum_plastic_work'],
                             implicit_work_kpa=a['work'],explicit_work_kpa=b['work'],
                             explicit_max_yield_drift_kpa=b['max_yield_drift'], branches=a['branches']))
            previous=(n,err)
            previous_solutions=(a['stress'],b['stress'])
            previous_adjacent=adjacent
        result['paths'][name]=rows
        point_histories[name]=dict(substeps_per_increment=levels[-1],
                                   initial_effective_stress_kpa=initial.tolist(),
                                   target_strains=[e.tolist() for e in targets],
                                   implicit_stress=a['stress'].tolist(),explicit_stress=b['stress'].tolist(),
                                   implicit_plastic_strain=a['plastic'].tolist(),explicit_plastic_strain=b['plastic'].tolist())
    wet=[]
    for z in (-.25,0.,.5,.75,1.):
        r=wet_update(-100*np.eye(3),np.diag([.002,.0021,-.01]),z)
        wet.append(dict(z=z,pore=r['pore'],
                        separation_error_kpa=float(np.linalg.norm(r['total']-r['effective']+r['pore']*np.eye(3))),
                        skeleton_vs_dry_error_kpa=float(np.linalg.norm(r['effective']-implicit(-100*np.eye(3),np.diag([.002,.0021,-.01]))['stress']))))
    result['wet']=wet
    all_hist={}
    for m,c in ((1,1),(2,1),(4,1),(1,2),(1,4),(8,8)):
        all_hist[f'm{m}-c{c}']=solve(m,c)
    result['column']={}
    for name,hist in all_hist.items():
        result['column'][name]=dict(final=hist[-1],
                                   max_free_force_residual=max(h['free_force_residual'] for h in hist),
                                   max_all_node_force_sum=max(h['all_node_force_sum'] for h in hist),
                                   max_energy_defect=max(abs(h['energy_defect']) for h in hist),
                                   max_local_energy_identity_error=max(abs(h['local_energy_identity_error']) for h in hist),
                                   minimum_bulk_dissipation_increment=min(h['minimum_bulk_dissipation_increment'] for h in hist),
                                   minimum_friction_dissipation_increment=min(h['minimum_friction_dissipation_increment'] for h in hist))
    result['partition_differences']={}
    for axis,pairs in dict(material=[('m1-c1','m2-c1'),('m2-c1','m4-c1')],contact=[('m1-c1','m1-c2'),('m1-c2','m1-c4')]).items():
        rows=[]
        for left,right in pairs:
            a,b=all_hist[left],all_hist[right]
            rows.append(dict(pair=[left,right],
                             plastic_strain_l2=float(np.linalg.norm(np.array([h['plastic_strain'] for h in a])-np.array([h['plastic_strain'] for h in b]))),
                             plastic_slip_l2_m=float(np.linalg.norm(np.array([h['plastic_slip'] for h in a])-np.array([h['plastic_slip'] for h in b]))),
                             energy_defect_difference=float(abs(a[-1]['energy_defect']-b[-1]['energy_defect']))))
        result['partition_differences'][axis]=rows
    result['elapsed_seconds']=time.perf_counter()-start
    out.mkdir(parents=True,exist_ok=True)
    (out/'diagnostics.json').write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
    (out/'column-history.json').write_text(json.dumps(all_hist,separators=(',',':'),allow_nan=False)+'\n',encoding='utf-8')
    (out/'stress-point-history.json').write_text(json.dumps(point_histories,separators=(',',':'),allow_nan=False)+'\n',encoding='utf-8')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--levels',type=int,nargs='+',default=[8,16,32,64])
    args=p.parse_args()
    if len(args.levels)<3 or args.levels[0]<1 or any(b!=2*a for a,b in zip(args.levels,args.levels[1:])):
        p.error('three or more positive, exactly doubling substep counts required')
    result=main(args.out,args.levels)
    print(json.dumps(dict(elapsed_seconds=result['elapsed_seconds'],accepted=None,C03='未達'),ensure_ascii=False))
