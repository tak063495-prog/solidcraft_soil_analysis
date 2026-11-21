"""Explicit opt-in bounded diagnostics. No production input or default changes."""
import argparse
import json
import os
from pathlib import Path
import shutil
import sys
import traceback
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from a02cause_binding import ROOT, study, reference, fields, geometry
from a02cause_metrics import metric, compare, common
from a02cause_projection import analyze as projection
CASES = [('transfer_force','transfer_force',False), ('stress_free','stress_free',False), ('stress_free-cycle','stress_free',True)]


def write(path, value):
    with Path(path).open('x', encoding='utf-8', newline='\n') as f:
        json.dump(value, f, ensure_ascii=False, indent=2, allow_nan=False);f.write('\n')


def fresh_reference(folder, steps, order, maximum):
    original = reference.fem.volume
    def volume(*args, **kw):
        if kw.get('independent') and kw.get('independent_order') == 5:kw = dict(kw, independent_order=order)
        return original(*args, **kw)
    reference.fem.volume = volume
    try:model = reference.Model(*(study.read(folder/(k+'.json')) for k in ('mesh','ground','case')))
    finally:reference.fem.volume = original
    with (folder/'independent-fields.jsonl').open('x',encoding='utf-8',newline='\n') as f:
        result = reference.run(model,steps,fields.Sink(f),contact_maximum=maximum)
    result.update(volume_gauss_legendre_order=order,diagnostic_contact_maximum_m=maximum,
                  uses_native_fields=False,new_global_equilibrium_solve=True)
    write(folder/'independent-result.json',result)


def local_replay(rows, system, parts):
    nq=len(system.points);jump=np.zeros((nq,3));traction=jump.copy();rest=jump.copy()
    total=np.zeros(nq);result=[];generation=None
    for row,f in rows:
        if generation is None or not np.array_equal(generation,f['generation']):
            jump=f['reference_jump'].copy();traction=f['initial_traction'].copy()
            rest=f['elastic_offset'].copy();total=np.zeros(nq);generation=f['generation'].copy()
        if np.all(f['joint_active']):
            after=f['jump'];total+=reference.contact_work(system,jump-rest,traction,after-rest,parts)[2]
            traction=reference.contact(system,jump-rest,traction,after-rest,parts)[0][:,0];jump=after.copy()
        result.append(total.copy())
    return np.array(result)


def collect(out, steps):
    summary=dict(passed=False,a02_complete=False,original_acceptance_all_met=False,new_limited_acceptance=False,
        production_defaults_unchanged=True,refusal_scope_unchanged=True,no_initial_force_clipping=True,
        local_replay_is_global_equilibrium=False,pc_executed=False,cloud_diagnostics_remaining=True,cases={},
        base='7b254b1cdbe9c288989631944cd698aa861f83ac')
    for name,policy,cycle in CASES:
        r=dict(same_grid={},native_grid_sensitivity={},reference_grid_sensitivity={},independent_controls={},local_replay={});by={}
        for n in steps:
            folder=out/f'n{n}'/name
            by[n]={side:list(fields.load_fields(folder/(side+'-fields.jsonl'))) for side in ('native','independent')}
            r['same_grid'][str(n)]=compare(by[n]['native'],by[n]['independent'],fields.CATEGORIES)
            r['same_grid'][str(n)]['legacy_comparison']=study.read(folder/'comparison.json')
        coarse=min(steps)
        for a,b in zip(steps,steps[1:]):
            for side,dest in [('native','native_grid_sensitivity'),('independent','reference_grid_sensitivity')]:
                r[dest][f'{a}-to-{b}']=compare(common(by[a][side],coarse),common(by[b][side],coarse),fields.CATEGORIES)
        if 4 in steps:
            for control in ('contact-quarter','volume-gl7'):
                r['independent_controls'][control]=compare(by[4]['independent'],list(fields.load_fields(out/'controls'/name/control/'independent-fields.jsonl')),fields.CATEGORIES)
            model=reference.Model(*study.fixture(policy,4,cycle))
            for parts in (1,2,4,8,16):
                a=local_replay(by[4]['native'],model.system,parts);b=local_replay(by[4]['independent'],model.system,parts)
                r['local_replay'][str(parts)]=dict(native_vs_independent=metric(a,b),
                    native_vs_saved_work=metric(a,np.stack([f['endpoint_work'] for _,f in by[4]['native']])),
                    independent_vs_saved_work=metric(b,np.stack([f['endpoint_work'] for _,f in by[4]['independent']])))
        summary['cases'][name]=r
    write(out/'measurements.json',summary);return summary


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--steps',type=int,nargs='+',default=[2,4,8]);p.add_argument('--reference',type=Path)
    p.add_argument('--order',type=int,choices=[5,7],default=5);p.add_argument('--maximum',type=float,default=.05)
    p.add_argument('--collect',action='store_true');a=p.parse_args()
    if any(n<2 or n>64 or n%2 for n in a.steps) or sorted(set(a.steps))!=a.steps:p.error('increasing even grids 2..64 only')
    if a.reference:
        if len(a.steps)!=1:p.error('one reference grid required')
        fresh_reference(a.reference,a.steps[0],a.order,a.maximum);return
    out=a.output.resolve()
    if a.collect:collect(out,a.steps);return
    out.mkdir(parents=True,exist_ok=False);binaries=out/'bin';binaries.mkdir()
    sources=sorted(p for folder in ('src','tools','examples','tests','interim') for p in (ROOT/folder).rglob('*')
                   if p.is_file() and p.suffix in ('.py','.rs') and '__pycache__' not in p.parts)
    fingerprints={p.relative_to(ROOT).as_posix():study.sha(p) for p in sources};write(out/'source-fingerprints.json',fingerprints)
    exe='.exe' if os.name=='nt' else '';fixed={}
    for name,original in [('solver',ROOT/'target/release'/('soil-fem-3d'+exe)),
        ('probe',ROOT/'target/release/examples'/('joint_a02_history_probe'+exe)),
        ('projection',ROOT/'target/release/examples'/('joint_a02_projection_probe'+exe)),
        ('rehash',ROOT/'target/release/examples'/('joint_ext_checkpoint_rehash'+exe))]:
        fixed[name]=binaries/(name+exe);shutil.copy2(original,fixed[name])
    write(out/'binaries.json',{k:study.sha(v) for k,v in fixed.items()});processes=[]
    def independent(folder,n,order=5,maximum=.05):
        process=study.invoke([sys.executable,Path(__file__).resolve(),'--output',out,'--reference',folder,
                              '--steps',n,'--order',order,'--maximum',maximum],folder/'independent.log',timeout=900)
        write(folder/'independent-process.json',process);processes.append(process)
        if process['exit_code']!=0 or process['timed_out']:raise RuntimeError(process)
    try:
        for n in a.steps:
            level=out/f'n{n}';level.mkdir()
            for name,policy,cycle in CASES:
                folder=level/name;study.make_inputs(folder,study.fixture(policy,n,cycle))
                print('native',n,name,flush=True);study.native(folder,fixed['probe'],n)
                print('reference',n,name,flush=True);independent(folder,n)
                write(folder/'comparison.json',fields.compare(folder/'native-fields.jsonl',folder/'independent-fields.jsonl'))
                write(folder/'equilibrium-audit.json',study.audits(folder));write(folder/'geometry-audit.json',geometry.audit(folder))
        if 4 in a.steps:
            controls=out/'controls';controls.mkdir()
            for name,policy,cycle in CASES:
                parent=controls/name;parent.mkdir()
                for kind,order,maximum in [('contact-quarter',5,.0125),('volume-gl7',7,.05)]:
                    folder=parent/kind;study.make_inputs(folder,study.fixture(policy,4,cycle))
                    print('reference-control',name,kind,flush=True);independent(folder,4,order,maximum)
            base=out/'n4';write(out/'refusals.json',study.refuse_inputs(base,fixed['solver'],4))
            for name,_,_ in CASES:
                write(base/name/'restart-observations.json',study.restart(base/name,fixed['probe'],4))
                write(base/name/'tamper-observations.json',study.tamper(base/name,fixed['probe'],fixed['rehash'],4))
            projection(out/'projection-observations',fixed['projection'])
        for path,digest in fingerprints.items():assert study.sha(ROOT/path)==digest,path
        write(out/'source-unchanged.json',dict(source_count=len(fingerprints),all_sha256_identical=True))
        collect(out,a.steps);write(out/'completion.json',dict(diagnostic_completed=True,a02_passed=False,processes=processes))
    except BaseException as e:
        write(out/'failure.json',dict(diagnostic_completed=False,a02_passed=False,error=repr(e),traceback=traceback.format_exc()));raise


if __name__=='__main__':main()
