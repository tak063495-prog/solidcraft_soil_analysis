"""Explicit, observation-only C-02 cause diagnosis. Production defaults untouched."""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
import time
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE))
from c02diag_mapping import ContactPlan
from c02diag_metrics import difference,triple,time_identity
from c02diag_transfer import norm
from c02diag_discretization import water_factory,alternate_system
from c02diag_manufactured import manufactured
from c02diag_path import primitives,path_summary,compare_primitives

FIELDS=('traction','plastic_slip','accumulated','plastic_dissipation','numerical_dissipation','endpoint_work','elastic_energy','jump','gap')
METHODS=('nearest','cubic','p2_l2')
BASE='026a4302fc5ed9b65eead8a493b1e11850b8b9ee'


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def sha(path):
    digest=hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):digest.update(chunk)
    return digest.hexdigest()


def write(path,value):
    with path.open('x',encoding='utf-8',newline='\n') as stream:
        json.dump(value,stream,ensure_ascii=False,indent=2,allow_nan=False,default=lambda x:x.tolist() if isinstance(x,np.ndarray) else x.item() if isinstance(x,np.generic) else str(x) if isinstance(x,Path) else (_ for _ in ()).throw(TypeError(type(x).__name__)));stream.write('\n')


def legacy():
    path=ROOT/'interim/2026-11-19-wet-joint-spatial-convergence/study_c02.py'
    spec=importlib.util.spec_from_file_location('c02diag_legacy19',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def load_contact(c,folder):
    manifest=read(folder/'manifest.json');rows=[];fields={k:[] for k in FIELDS}
    for chunk in manifest['chunks']:
        path=folder/chunk['file']
        if sha(path)!=chunk['sha256']:raise ValueError('changed accepted field')
        if chunk['start']!=len(rows) or chunk['count']!=len(chunk['rows']):raise ValueError('history chunk order')
        with np.load(path,allow_pickle=False) as data:
            for key in FIELDS:
                value=c.old.codec.decode(data[key],manifest['encoding']) if manifest.get('encoding') else data[key]
                if len(value)!=len(chunk['rows']) or not np.isfinite(value).all():raise ValueError('invalid contact field')
                fields[key].append(value)
        rows.extend(chunk['rows'])
    if len(rows)!=manifest['count']:raise ValueError('incomplete contact history')
    return rows,{k:np.concatenate(v) for k,v in fields.items()}


def common(source,c):
    meshes=[read(source/f'mesh-h{h}.json') for h in range(3)]
    ground=read(source/'ground.json');case=read(source/'stages.json')
    systems=[c.old.Water(m,ground,case).contact for m in meshes]
    positions=np.asarray([p['position_m'] for s in systems for p in s.points])
    roots=np.asarray([p['face_index']//4**h for h,s in enumerate(systems) for p in s.points])
    weights=np.asarray([p['weight_m2']/3 for s in systems for p in s.points])
    return meshes,ground,case,systems,positions,roots,weights


def sampled(rows,fields,steps):
    ids=[]
    for i,r in enumerate(rows):
        v=float(r['progress'])*steps
        if abs(v-round(v))<1e-10:ids.append(i)
    selected=[rows[i] for i in ids]
    keys=[(r['stage'],round(float(r['progress'])*steps)) for r in selected]
    if len(keys)!=7*steps or len(set(keys))!=len(keys):raise ValueError('missing common accepted times')
    return keys,{k:v[ids] for k,v in fields.items()}


def area_shape(weights,values):
    return weights.reshape((1,len(weights))+(1,)*(values.ndim-2))


def label_brackets(rows,gaps):
    result={}
    for stage in dict.fromkeys(r['stage'] for r in rows):
        ids=[i for i,r in enumerate(rows) if r['stage']==stage]
        flags=gaps[ids]>0;events=[]
        for i in range(1,len(ids)):
            for point in np.flatnonzero(flags[i]!=flags[i-1]):
                events.append(dict(point=int(point),left=float(rows[ids[i-1]]['progress']),right=float(rows[ids[i]]['progress']),to_open=bool(flags[i,point])))
        result[stage]=events
    return result


def measure(source,out,variants=None):
    out.mkdir(parents=True,exist_ok=False);c=legacy()
    meshes,ground,case,systems,positions,roots,weights=common(source,c)
    plans={(h,m):ContactPlan(meshes[h],systems[h],h,positions,roots,m) for h in range(3) for m in METHODS}
    cache={};source_hashes={}
    routes=[(h,g) for h in (0,1) for g in (32,64,128,256)]+[(2,g) for g in (32,64)]
    for h,g in routes:
        folder=source/f'h{h}-g{g}'/'native'
        cache[h,g]=load_contact(c,folder)
        source_hashes[f'h{h}-g{g}/native/manifest.json']=sha(folder/'manifest.json')
    rows={};mapped={};labels={};primitive_cache={}
    for (h,g),(r,f) in cache.items():
        rows[h,g]=r
        primitive_cache[h,g]=primitives(f['jump'],f['traction'],f['plastic_slip'],systems[h].normal,systems[h].kt)
        labels[h,g]=plans[h,'nearest'].history(f['gap'])>0
        for m in METHODS:
            mapped[h,g,m]={k:plans[h,m].history(v) for k,v in f.items()}
    output=dict(passed=False,status='unmet',base_commit=BASE,definition_changed_for_acceptance=False,source_manifests_sha256=source_hashes,common_points=len(positions),methods=list(METHODS),spatial={},time_path={},reference={},interpolation_sensitivity={},manufactured=manufactured())
    output['path_primitive']={f'h{h}-g{g}':path_summary(primitive_cache[h,g],f['plastic_slip']) for (h,g),(r,f) in cache.items()}
    output['path_primitive_time']={}
    for h in range(3):
        for ga,gb in ((32,64),(64,128),(128,256)):
            if (h,gb) not in cache:continue
            _,fine=sampled(rows[h,gb],primitive_cache[h,gb],ga)
            _,pfine=sampled(rows[h,gb],{'plastic_slip':cache[h,gb][1]['plastic_slip']},ga)
            output['path_primitive_time'][f'h{h}-g{ga}-{gb}']=compare_primitives(primitive_cache[h,ga],fine,cache[h,ga][1]['plastic_slip'],pfine['plastic_slip'])
    for g in (32,64):
        for m in METHODS:
            values=[mapped[h,g,m] for h in range(3)]
            output['spatial'][f'fixed{g}-{m}']={k:dict(point=triple([v[k] for v in values]),area=triple([v[k] for v in values],area_shape(weights,values[0][k]))) for k in FIELDS}
    for (ha,hb),ga,gb in [((0,1),32,64),((1,2),32,64),((0,1),64,256),((0,1),128,256)]:
        for m in METHODS:
            selections=[];label_rows=[];keys=[]
            for h,g in [(ha,ga),(hb,ga),(ha,gb),(hb,gb)]:
                ky,v=sampled(rows[h,g],mapped[h,g,m],ga);keys.append(ky);selections.append(v)
                _,ll=sampled(rows[h,g],{'label':labels[h,g]},ga);label_rows.append(ll['label'])
            if any(k!=keys[0] for k in keys):raise ValueError('time history alignment')
            output['time_path'][f'h{ha}{hb}-g{ga}-{gb}-{m}']={k:time_identity(*(s[k] for s in selections),label_rows) for k in FIELDS}
    for (h,g),(r,f) in cache.items():
        key=f'h{h}-g{g}';output['interpolation_sensitivity'][key]={}
        for m in ('nearest','p2_l2'):
            output['interpolation_sensitivity'][key][m]={k:difference(mapped[h,g,m][k],mapped[h,g,'cubic'][k]) for k in FIELDS}
        for m in METHODS:
            v=mapped[h,g,m];gap=v['gap'];traction=v['traction'];normal=np.asarray([systems[0].normal[0]]*len(positions))
            pressure=-np.sum(traction*normal,axis=-1);tau=traction+pressure[:,:,None]*normal
            limit=systems[0].mu[0]*np.maximum(pressure,0)
            output['interpolation_sensitivity'][key][m+'_law_observation']=dict(open_point_times=int((gap>0).sum()),traction_in_open_max_kpa=float(np.max(np.linalg.norm(traction[gap>0],axis=-1),initial=0)),friction_excess_max_kpa=float(np.max(np.linalg.norm(tau,axis=-1)-limit,initial=0)),negative_plastic_dissipation_max=float(np.max(-v['plastic_dissipation'],initial=0)))
        if h<2:
            ir,independent=load_contact(c,source/f'h{h}-g{g}'/'independent')
            if [(x['stage'],x['progress']) for x in ir]!=[(x['stage'],x['progress']) for x in r]:raise ValueError('independent schedule mismatch')
            output['reference'][key]={k:difference(f[k],independent[k]) for k in FIELDS}
    switches={f'h{h}-g{g}':label_brackets(r,plans[h,'nearest'].history(cache[h,g][1]['gap'])) for (h,g),r in rows.items()}
    write(out/'switch-brackets.json',switches)
    np.savez_compressed(out/'common-points.npz',positions=positions,roots=roots,area_weights=weights)
    # Preserve aligned arrays for inexpensive inspection without any new solve.
    for g in (32,64):
        np.savez_compressed(out/f'point-differences-g{g}.npz',**{f'{m}_{k}_h{h}{h+1}':mapped[h,g,m][k]-mapped[h+1,g,m][k] for m in METHODS for k in FIELDS for h in (0,1)})
    write(out/'measurements.json',output)
    return output


def solve_variant(source,out,h,g,rule,interpolation):
    if h not in (0,1) or g not in (32,64):raise ValueError('cloud diagnostic limited to h0/h1, g32/g64')
    out.mkdir(parents=True,exist_ok=False);c=legacy()
    mesh=read(source/f'mesh-h{h}.json');ground=read(source/'ground.json');case=read(source/'stages.json')
    reference=c.old.probe
    default=reference.Water
    reference.Water=water_factory(default,rule,interpolation)
    writer=c.Writer(out/'independent');started=time.monotonic()
    try:
        result=reference._reference(mesh,ground,case,g,1,7,False,True,writer,maximum_m=5e-5,maximum_dofs=15000)
        writer.close()
        system=reference.Water(mesh,ground,case).contact
        write(out/'result.json',result)
        write(out/'layout.json',dict(points=system.points,rule=rule,interpolation=interpolation,h=h,steps=g,input_only=True,production_changed=False))
        write(out/'timing.json',dict(elapsed_seconds=time.monotonic()-started,complete=True))
    except BaseException as error:
        write(out/'failure.json',dict(passed=False,complete=False,error=repr(error),elapsed_seconds=time.monotonic()-started))
        raise
    finally:reference.Water=default


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',required=True,choices=('measure','solve-variant','manufactured'))
    parser.add_argument('--source',type=Path)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--h',type=int,default=0);parser.add_argument('--steps',type=int,default=32)
    parser.add_argument('--rule',choices=('jacobi4','legendre4'),default='jacobi4')
    parser.add_argument('--interpolation',choices=('tri6','corner3'),default='tri6')
    args=parser.parse_args();source=args.source.resolve() if args.source else None;out=args.output.resolve()
    if args.mode=='manufactured':
        out.mkdir(parents=True,exist_ok=False);write(out/'manufactured.json',manufactured())
    elif args.mode=='measure':measure(source,out)
    else:solve_variant(source,out,args.h,args.steps,args.rule,args.interpolation)

if __name__=='__main__':main()
