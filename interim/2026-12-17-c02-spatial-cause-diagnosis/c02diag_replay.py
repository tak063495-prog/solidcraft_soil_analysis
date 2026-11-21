"""Finite saved-polyline constitutive reference. Never a global equilibrium solve."""
import argparse
import importlib.util
from pathlib import Path
import sys
import time
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import c02diag_driver as d
from c02diag_mapping import ContactPlan
from c02diag_labels import catalog,masks,budget
from c02diag_transfer import norm,dot

FIELDS=('plastic_slip','traction','plastic_dissipation','elastic_energy')


def run(source,out):
    out.mkdir(parents=True,exist_ok=False);c=d.legacy()
    meshes,ground,case,systems,positions,roots,weights=d.common(source,c)
    path=d.ROOT/'interim/2026-12-15-c01-local-path-decomposition/c01_path_replay.py'
    spec=importlib.util.spec_from_file_location('c02diag_c01_constitutive_replay',path)
    local=importlib.util.module_from_spec(spec);spec.loader.exec_module(local)
    routes=[(h,g) for h in (0,1) for g in (32,64,128,256)]+[(2,g) for g in (32,64)]
    plans={(h,m):ContactPlan(meshes[h],systems[h],h,positions,roots,m) for h in range(3) for m in ('cubic','nearest')}
    raw={};rows={};reference={};catalogs={};timings=[]
    for h,g in routes:
        r,fields=d.load_contact(c,source/f'h{h}-g{g}'/'native');rows[h,g]=r
        raw[h,g]={k:plans[h,'cubic'].history(fields[k]) for k in FIELDS}
        catalogs[h,g]=catalog(r,plans[h,'nearest'].history(fields['gap']))
        system=systems[h]
        for multiple in (16,32):
            started=time.monotonic();jump=np.zeros_like(system.normal);q=jump.copy();plastic=jump.copy();diss=np.zeros(len(jump));saved={k:[] for k in FIELDS};counts=[];minimum=0.;balance=0.;work=diss.copy();numerical=diss.copy()
            for target in fields['jump']:
                lower=np.maximum(1,np.ceil(np.sqrt(local.dot3(target-jump,target-jump))/5e-5)).astype(int)
                count=lower*multiple
                q,dp,dd,dw,dn=local.integrate(jump,q,target,system.normal,system.kn,system.kt,system.mu,count)
                plastic+=dp;diss+=dd;work+=dw;numerical+=dn
                energy=local.energy(q,system.normal,system.kn,system.kt)
                for k,value in zip(FIELDS,(plastic,q,diss,energy),strict=True):saved[k].append(value.copy())
                minimum=min(minimum,float(np.min(dd)));balance=max(balance,float(np.max(abs(work-energy-diss-numerical))));counts.append(count)
                jump=target
            saved={k:np.stack(v) for k,v in saved.items()}
            np.savez_compressed(out/f'h{h}-g{g}-m{multiple}.npz',**saved,counts=np.stack(counts))
            reference[h,g,multiple]={k:plans[h,'cubic'].history(v) for k,v in saved.items()}
            timing=dict(h=h,steps=g,multiple=multiple,elapsed_seconds=time.monotonic()-started,minimum_increment_dissipation=minimum,maximum_discrete_work_balance_residual=balance)
            timings.append(timing);print(timing,flush=True)
            d.write(out/f'h{h}-g{g}-m{multiple}-audit.json',timing)
    result=dict(passed=False,local_replay_is_equilibrium=False,finite_reference_not_exact_or_error_bound=True,reference_count_definition='ceil(norm(accepted Delta jump)/5e-5) times 16 or 32; not observed frozen Newton counts',c01_integrator_sha256=d.sha(path),timings=timings,time_spatial_change={})
    for (ha,hb),ga,gb in [((0,1),32,64),((1,2),32,64),((0,1),64,256),((0,1),128,256)]:
        route=[(ha,ga),(hb,ga),(ha,gb),(hb,gb)]
        mask,duration=masks([catalogs[x] for x in route],c.STAGES,ga)
        f=[d.sampled(rows[h,g],raw[h,g],ga)[1] for h,g in route]
        refs={m:[d.sampled(rows[h,g],reference[h,g,m],ga)[1] for h,g in route] for m in (16,32)}
        key=f'h{ha}{hb}-g{ga}-{gb}';result['time_spatial_change'][key]={}
        for field in FIELDS:
            def change(values):return (values[2][field]-values[3][field])-(values[0][field]-values[1][field])
            total=change(f);r16=change(refs[16]);r32=change(refs[32]);ss=dot(total.ravel(),total.ravel())
            _,parts=budget(r32,mask,c.STAGES,ga)
            components=dict(local_integration=total-r32,opening_label=parts['opening_label'],other_saved_path=parts['other_path'])
            result['time_spatial_change'][key][field]=dict(total_absolute_l2=norm(total),components={name:dict(absolute_l2=norm(value),signed_projection=dot(value.ravel(),total.ravel())/ss if ss else None) for name,value in components.items()},cross_inner_products={a+'/'+b:dot(x.ravel(),y.ravel()) for a,x in components.items() for b,y in components.items()},finite_16_to_32_change_over_total=norm(r32-r16)/max(norm(total),1e-12),identity_residual_l2=norm(total-sum(components.values())))
    d.write(out/'local-budgets.json',result)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.source.resolve(),a.output.resolve())
