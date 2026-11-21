"""A-04 fresh soil construction solves, separate process reference and checkpoints.

Production restrictions stay intact. All logs/failures are retained. Never use a
saved solution as an input to the independent solver or count reaudits as solves.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import argparse
import copy
import hashlib
import itertools
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
import zipfile
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE));sys.path.insert(0,str(ROOT/'tools'))
from soil_reference import Model,SAMPLE_KEYS

BASE='04dd721ab82dd2c044569edf9e0b775b605cf2ef'
CATEGORY={'cell_ids','face_keys','active_elements','material_return'}


def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def convert(v):
    if isinstance(v,np.ndarray):return v.astype(float).tolist() if v.dtype.kind=='f' else v.tolist()
    if isinstance(v,np.generic):return v.item()
    raise TypeError(type(v).__name__)
def write(p,v):
    with Path(p).open('x',encoding='utf-8') as f:json.dump(v,f,default=convert,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
def exact(v):
    def strip(x):
        if isinstance(x,dict):return {k:strip(a) for k,a in x.items() if k!='elapsed_seconds'}
        if isinstance(x,list):return [strip(a) for a in x]
        return x
    return json.dumps(strip(v),sort_keys=True,separators=(',',':'),allow_nan=False)
def invoke(cmd,log,timeout=900):
    cmd=list(map(str,cmd));start=time.monotonic()
    with Path(log).open('x',encoding='utf-8') as f:
        f.write(repr(cmd)+'\n');f.flush();p=subprocess.Popen(cmd,stdout=f,stderr=subprocess.STDOUT);timed=False
        try:code=p.wait(timeout=int(os.environ.get('A04_PROCESS_BUDGET_SECONDS',timeout)))
        except subprocess.TimeoutExpired:timed=True;p.kill();code=p.wait()
    result=dict(command=cmd,pid=p.pid,exit_code=code,timed_out=timed,elapsed_seconds=time.monotonic()-start)
    write(Path(log).with_suffix('.process.json'),result);return result


def fixture(n,saturated):
    mesh=read(ROOT/'examples/seepage-column-tet4-mesh.json');ground=read(ROOT/'examples/ground-one-way-mc.json');base=read(ROOT/'examples/one-way-mc.json')
    for mat in base['flow']['materials']:
        mat['saturated_conductivity_m_s']=(np.eye(3)*.02).tolist();mat['retention']['saturated_specific_storage_per_m']=.02
    levels=[2.,2.2,2.2,2.1] if saturated else [.3,.5,.1,.4]
    base['flow']['initial_head']={'head_m':levels[0]}
    base['flow']['analysis']={'kind':'transient','periods':[dict(name=f'p{i}',end_time_s=float(i+1),boundaries=[dict(face_set=face,condition=dict(kind='head',value=dict(head_m=level))) for face in ('bottom','top' if i<2 else 'floor')]) for i,level in enumerate(levels)]}
    base['flow']['settings'].update(initial_time_step_s=1/n,maximum_time_step_s=1/n,minimum_time_step_s=1/(1024*n),maximum_steps=10000,maximum_free_faces=1000,store_integration_points=True,linear_solver=dict(method='dense_lu'))
    base['mechanics']['settings'].update(kinematics='side_rollers',maximum_free_dofs=1000,relative_force_tolerance=1e-11,relative_linear_tolerance=1e-12,tangent_perturbation=1e-6,initial_increment=.25,minimum_increment=1/1024,linear_solver=dict(method='dense_lu'))
    base['mechanics']['pressure_on_original_surface_only']=True
    base['exterior_periods']=[dict(name=f'p{i}',reservoirs=[dict(face_set='floor',level_m=levels[i])] if saturated and i>=2 else []) for i in range(4)]
    case=dict(schema='soil-fem-one-way-construction/1',one_way=base,excavations=[dict(name='cut',before_period='p2',bodies=['UpperClay'],water_boundary_policy='explicit_period_boundaries',cut_traction_policy='explicit_period_reservoirs',initial_increment=.25)])
    return mesh,ground,case


def sink_file(path):
    stream=path.open('x',encoding='utf-8')
    def sink(row,fields):
        stream.write(json.dumps(dict(row=row,fields=fields),default=convert,allow_nan=False,separators=(',',':'))+'\n');stream.flush()
    return stream,sink


def independent(folder):
    m=Model(*[read(folder/(k+'.json')) for k in ('mesh','ground','case')]);stream,sink=sink_file(folder/'independent-fields.jsonl')
    try:result=m.solve(sink)
    finally:stream.close()
    write(folder/'independent-result.json',result)


def native_fields(folder):
    mesh,ground,case=[read(folder/(k+'.json')) for k in ('mesh','ground','case')];m=Model(mesh,ground,case)
    run=read(folder/'native.json');assert run['complete'] and all(a['accepted'] for a in run['attempts']+run['initialization']['attempts'])
    assert len(run['frames'])==4*round(1/case['one_way']['flow']['settings']['maximum_time_step_s'])+4
    initial=run['initial_water_m3'];stream,sink=sink_file(folder/'native-fields.jsonl');rows=[];previous_u=np.zeros_like(m.xyz)
    try:
        allframes=[(f,True) for f in run['initialization']['frames'][1:]]+[(f,False) for f in run['frames']]
        strain=np.zeros_like(m.seed)
        for frame,init in allframes:
            s=frame['state'];h=s['hydraulic'];c=s['mechanical'];native=s['fluid'];active=np.array(c['active_elements']);hm=m.full if active.all() else m.cut
            elems=hm.elements;raw=np.array([[p[k] for k in SAMPLE_KEYS] for p in native['hydraulic_points']]);mp=np.array([[p[k] for k in SAMPLE_KEYS] for p in native['mechanical_points']])
            assert [(p['element_index'],p['local_index']) for p in native['hydraulic_points']]==[(int(e),j) for e in elems for j in range(len(hm.weight[0]))]
            assert [p['element_index'] for p in native['mechanical_points']]==elems.tolist()
            faceinfo=run['initialization']['hydraulic_faces'] if init else frame['hydraulic_faces']
            lookup={tuple(sorted(f['corners'])):j for j,f in enumerate(faceinfo)}
            v=dict(h=np.array([h['relative_cell_head_m'][e] for e in elems]),trace=np.array([h['relative_face_head_m'][lookup[k]] for k in hm.keys]),storage=np.array([h['cell_storage_offset_m3'][e] for e in elems]),flux=np.array([h['cell_outflow_m3_s'][e] for e in elems]))
            assert all(h['relative_cell_head_m'][e] is None for e in np.flatnonzero(~active))
            for e in np.flatnonzero(~active):assert h['cell_storage_offset_m3'][e] is None and h['cell_outflow_m3_s'][e] is None
            fluid=dict(hydraulic_samples=raw,mechanical_samples=mp,hydraulic_positions=np.array([p['position_m'] for p in native['hydraulic_points']]),hydraulic_weights=np.array([p['weight_m3'] for p in native['hydraulic_points']]),
                fluid_internal=np.array(native['fluid_internal_nodal_force_kn']).reshape(-1,3),water_body=np.array(native['water_body_nodal_load_kn']).reshape(-1,3),exterior=np.array(native['exterior_nodal_load_kn']).reshape(-1,3),
                physical=np.array(native['physical_nodal_load_kn']).reshape(-1,3),rhs=np.array(native['equivalent_mechanical_rhs_kn']).reshape(-1,3),water=native['water_equivalent_mass_volume_m3'])
            u=np.array(c['displacement_m']);states=c['material_states'];stress=np.array([a['effective_stress_compression_positive_kpa'] for a in states]);void=np.array([a['void_ratio'] for a in states])
            for e in np.flatnonzero(active):strain[e]=m.q[e]['b'][0]@u.ravel()[m.ids[e]]
            load=np.array(c['applied_nodal_load_kn']).reshape(-1,3)
            if init:
                progress=frame['initialization_progress'];holding=(m.force(m.seed,active).reshape(-1,3)-fluid['rhs'])*(1-progress);removed=inflow=0.
                row=dict(phase='initialization',name=frame['name'],time_s=h['time_s'],progress=progress)
            else:
                holding=np.array(frame['holding_load_kn']).reshape(-1,3);removed=frame['cumulative_removed_water_m3'];inflow=frame['cumulative_net_water_inflow_m3']
                row=dict(phase=frame['phase'],name=frame['name'],time_s=h['time_s'],progress=frame['progress'])
                assert frame['mechanical_audit']['relative_effective_rhs_residual']<=1e-8
                if frame['hydraulic_balance'] is not None:assert frame['hydraulic_balance']['maximum_scaled_residual']<=1.
            f=m.fields(hm,v,fluid,u,stress,strain,void,active,load,holding,removed,inflow,initial)
            f['plastic_strain']=np.array([a['plastic_strain_compression_positive'] for a in states]);f['equivalent_plastic_strain']=np.array([a['equivalent_plastic_strain'] for a in states]);f['plastic_work']=np.array([a['plastic_work_kpa'] for a in states]);f['material_return']=np.array([0 if a['last_return']=='elastic' else 1 for a in states])
            assert not np.any(f['plastic_strain']) and not np.any(f['material_return'])
            if not init:assert abs(float(f['mass_defect_m3'])-frame['signed_global_water_budget_error_m3'])<1e-10
            sink(row,f);rows.append(row);previous_u=u
    finally:stream.close()
    return dict(accepted_states=len(rows),all_coordinates=rows,native_water_fields='production Fluid',native_effective_internal='separate reassembly of saved production stress; diagnostic, not a fresh solve',contact=None,hydraulic_joint_mapping=None)


def fields(path):
    with Path(path).open() as f:
        for line in f:
            v=json.loads(line);yield v['row'],{k:np.asarray(x) for k,x in v['fields'].items()}

def compare(left,right,common=False,phase=None):
    a=list(fields(left));b=list(fields(right));target={exact(r):f for r,f in b};stats={};points=[]
    for row,fa in a:
        if phase is not None and row['phase']!=phase:continue
        key=exact(row)
        if common and key not in target:continue
        assert key in target,row
        fb=target[key];assert set(fa)==set(fb);point={}
        for k,x in fa.items():
            y=fb[k];assert x.shape==y.shape and np.isfinite(x).all() and np.isfinite(y).all(),(row,k,x.shape,y.shape)
            if k in CATEGORY:
                v=dict(mismatches=int(np.count_nonzero(x!=y)),components=x.size)
                old=stats.setdefault(k,dict(mismatches=0,components=0));old['mismatches']+=v['mismatches'];old['components']+=x.size
            else:
                d=np.asarray(x,dtype=float)-y;ds=float(np.square(d).sum());rs=float(np.square(y).sum());maximum=float(np.max(np.abs(d),initial=0.))
                v=dict(maximum_absolute=maximum,relative_l2=ds**.5/max(rs**.5,1e-12),difference_squared=ds,reference_squared=rs,components=x.size)
                old=stats.setdefault(k,dict(maximum_absolute=0.,difference_squared=0.,reference_squared=0.,components=0));old['maximum_absolute']=max(old['maximum_absolute'],maximum);old['difference_squared']+=ds;old['reference_squared']+=rs;old['components']+=x.size
            point[k]=v
        points.append(dict(row=row,fields=point))
    assert points
    for v in stats.values():
        if 'difference_squared' in v:v['relative_l2']=v['difference_squared']**.5/max(v['reference_squared']**.5,1e-12)
    return dict(accepted_pairs=len(points),field_count=len(stats),all_fields=stats,all_points=points,interpolation=False,numerical_acceptance=False)


def restart(folder,solver,n,rehash):
    out=folder/'restarts';out.mkdir();baseline=read(folder/'native.json');rows=[]
    # Count excludes initialization in construction; all transient/release states are frames.
    for count in (n+n//2,2*n+2,3*n+4+n//2):
        cp=out/f'checkpoint-{count}.json';dest=out/f'resumed-{count}.json'
        save=invoke([solver,'checkpoint-one-way-construction',folder/'mesh.json',folder/'ground.json',folder/'case.json',count,cp],out/f'save-{count}.log')
        resume=invoke([solver,'resume-one-way-construction',folder/'mesh.json',folder/'ground.json',folder/'case.json',cp,dest],out/f'resume-{count}.log')
        assert save['exit_code']==resume['exit_code']==0 and not save['timed_out'] and not resume['timed_out']
        assert save['pid']!=resume['pid'] and exact(read(dest))==exact(baseline)
        saved=read(cp);assert len(saved['run']['frames'])==count
        assert exact(read(dest)['frames'][count:])==exact(baseline['frames'][count:])
        rows.append(dict(accepted_count=count,saved_time_s=saved['run']['last_committed']['hydraulic']['time_s'],phase=saved['run']['frames'][-1]['phase'],checkpoint_sha256=sha(cp),save=save,resume=resume,suffix_states_exact=len(baseline['frames'])-count,whole_run_exact_except_elapsed=True))
    tamper=folder/'tamper';tamper.mkdir();cp=out/f'checkpoint-{2*n+2}.json';original=read(cp);mutations=[]
    for mode in ('head','trace','flux','storage','time','active','fluid-load','ledger'):
        bad=copy.deepcopy(original);s=bad['run']['last_committed']
        if mode=='head':s['hydraulic']['relative_cell_head_m'][0]+=.01
        elif mode=='trace':s['hydraulic']['relative_face_head_m'][0]+=.01
        elif mode=='flux':s['hydraulic']['cell_outflow_m3_s'][0][0]+=.01
        elif mode=='storage':s['hydraulic']['cell_storage_offset_m3'][0]+=.01
        elif mode=='time':s['hydraulic']['time_s']+=.01
        elif mode=='active':s['mechanical']['active_elements'][-1]=True
        elif mode=='fluid-load':s['fluid']['equivalent_mechanical_rhs_kn'][0]+=.01
        else:bad['run']['runtime']['removed_water_m3']+=.01
        raw=tamper/(mode+'-raw.json');sealed=tamper/(mode+'-rehashed.json');write(raw,bad)
        seal=invoke([rehash,raw,sealed],tamper/(mode+'-rehash.log'));assert seal['exit_code']==0
        assert read(sealed)['run_sha256']!=original['run_sha256']
        assert read(sealed)['run']['solver_binary_sha256']==sha(solver)
        target=tamper/(mode+'-result.json');p=invoke([solver,'resume-one-way-construction',folder/'mesh.json',folder/'ground.json',folder/'case.json',sealed,target],tamper/(mode+'.log'))
        diagnostic=(tamper/(mode+'.log')).read_text();assert p['exit_code']==1 and not p['timed_out'] and not target.exists()
        assert 'construction checkpoint input, executable, domain, archive or history changed' in diagnostic,diagnostic
        mutations.append(dict(mode=mode,outer_hash_recomputed=True,producer_binary_sha256=sha(solver),process=p,reason=diagnostic))
    return dict(restarts=rows,tampering=mutations)


def snapshot(out):
    names=subprocess.check_output(['git','ls-files'],cwd=ROOT,text=True).splitlines()
    selected=[f for f in names if Path(f).suffix in ('.py','.rs','.toml') or f in ('Cargo.lock','requirements-reference.lock.txt')]
    selected+=[str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix=='.py' or p.name=='SPEC.ja.md'];selected+=['examples/a04_checkpoint_rehash.rs']
    fp={f:sha(ROOT/f) for f in sorted(set(selected))}
    write(out/'source-fingerprints.json',fp)
    with zipfile.ZipFile(out/'source-snapshot.zip','x',zipfile.ZIP_DEFLATED) as z:
        for f in fp:z.write(ROOT/f,f)
    return fp


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--steps',type=int,nargs='+',default=[4,8,16]);p.add_argument('--independent',type=Path);p.add_argument('--development',action='store_true');a=p.parse_args()
    if a.independent:independent(a.independent);return
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);assert all(n>=2 and n&(n-1)==0 for n in a.steps)
    fp=snapshot(out);binpath=out/'bin';binpath.mkdir();exe='.exe' if os.name=='nt' else ''
    solver=binpath/('solver'+exe);rehash=binpath/('rehash'+exe)
    shutil.copy2(ROOT/('target/release/soil-fem-3d'+exe),solver);shutil.copy2(ROOT/('target/release/examples/a04_checkpoint_rehash'+exe),rehash)
    write(out/'binaries.json',dict(solver_sha256=sha(solver),rehash_sha256=sha(rehash)))
    report=dict(passed=False,a04_complete=False,base=BASE,cases={},development=a.development)
    try:
        for saturated in (False,True):
            label='saturated' if saturated else 'VG'
            for n in a.steps:
                folder=out/f'{label}-n{n}';folder.mkdir()
                for name,value in zip(('mesh','ground','case'),fixture(n,saturated)):write(folder/(name+'.json'),value)
                proc=invoke([solver,'solve-one-way-construction',folder/'mesh.json',folder/'ground.json',folder/'case.json',folder/'native.json'],folder/'native.log');assert proc['exit_code']==0 and not proc['timed_out'],proc
                independent_proc=invoke([sys.executable,Path(__file__),'--output',out,'--independent',folder],folder/'independent.log');assert independent_proc['exit_code']==0 and not independent_proc['timed_out'],independent_proc
                observation=native_fields(folder);comparison=compare(folder/'native-fields.jsonl',folder/'independent-fields.jsonl');write(folder/'comparison.json',comparison)
                allf=comparison['all_fields'];assert all(v.get('mismatches',0)==0 for v in allf.values())
                assert all(allf[k]['relative_l2']<1e-7 for k in ('u','stress'))
                assert all(allf[k]['maximum_absolute']<1e-10 for k in ('effective_internal','fluid_internal','water_body','exterior','physical_load','equivalent_rhs','applied_load'))
                result=dict(native_process=proc,independent_process=independent_proc,observation=observation,comparison=comparison,existing_bounded_gates_passed=True)
                if n==a.steps[-1] and not a.development:result.update(restart(folder,solver,n,rehash))
                report['cases'][folder.name]=result;print('completed',folder.name,comparison['accepted_pairs'],flush=True)
        sensitivity={}
        for label in ('VG','saturated'):
            for n1,n2 in zip(a.steps,a.steps[1:]):
                for kind in ('native','independent'):
                    sensitivity[f'{label}/{kind}/{n1}__{n2}']=compare(out/f'{label}-n{n1}'/(kind+'-fields.jsonl'),out/f'{label}-n{n2}'/(kind+'-fields.jsonl'),True)
        transient_sensitivity={}
        for label in ('VG','saturated'):
            for n1,n2 in zip(a.steps,a.steps[1:]):
                for kind in ('native','independent'):
                    transient_sensitivity[f'{label}/{kind}/{n1}__{n2}']=compare(out/f'{label}-n{n1}'/(kind+'-fields.jsonl'),out/f'{label}-n{n2}'/(kind+'-fields.jsonl'),True,'transient')
        write(out/'transient-time-sensitivity.json',transient_sensitivity);report['transient_time_sensitivity']=transient_sensitivity
        write(out/'time-sensitivity.json',sensitivity);report['time_sensitivity']=sensitivity
        assert all(sha(ROOT/f)==v for f,v in fp.items());report['sources_unchanged']=True
        report['fresh_native_paths']=report['fresh_independent_paths']=len(report['cases']);report['contact']=None
        write(out/'measurements.json',report)
    except BaseException as e:
        write(out/'failure.json',dict(passed=False,complete=False,error=repr(e),completed_cases=list(report['cases'])));raise

if __name__=='__main__':main()
