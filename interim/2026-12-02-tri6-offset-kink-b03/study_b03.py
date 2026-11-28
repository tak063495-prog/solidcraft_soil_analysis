"""Fresh production/public-API and input-only independent runs, all accepted points."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import argparse,shutil,sys,json
from pathlib import Path
from b03_cases import HERE,B02,ROOT,make,read,write
from study import command,sha,exact
from compare import compare


def run(out,binary,kinds,level_max=1):
    out.mkdir(parents=True,exist_ok=False)
    frozen=out/'solver'/binary.name;frozen.parent.mkdir();shutil.copy2(binary,frozen)
    sources=[p for base in ('src','tools','tests','examples',str(HERE.relative_to(ROOT)),str(B02.relative_to(ROOT)))
             for p in (ROOT/base).rglob('*') if p.suffix in ('.rs','.py','.md') or (base=='examples' and p.suffix=='.json')]
    sources += [ROOT/p for p in ('Cargo.toml','Cargo.lock','rust-toolchain.toml')]
    before={str(p.relative_to(ROOT)):sha(p) for p in sources};write(out/'source-binding-before.json',before)
    plan=[('planar',0,4)]
    for kind in kinds:
        plan += [(kind,0,4)]+[(kind,l,4) for l in range(1,level_max+1)]+[(kind,0,s) for s in (8,16)]
    records=[];restarts=[]
    for kind,level,steps in plan:
        name=f'{kind}-h{level}-s{steps}';folder=out/name;folder.mkdir()
        m,g,c,a=make(kind,level,steps)
        for k,v in [('mesh',m),('ground',g),('stages',c),('refinement',a)]:write(folder/(k+'.json'),v)
        base=[folder/(k+'.json') for k in ('mesh','ground','stages')]
        stops=[2*steps+steps//2]
        if level==0 and steps==4 and kind=='kink-zero':stops=[steps//2,steps+steps//2,2*steps+steps//2,3*steps+steps//2,4*steps+steps//2]
        native=command([frozen,'solve',*base,folder/'native',','.join(map(str,stops))],folder/'native.log')
        independent=command([sys.executable,B02/'reference.py',folder,folder/'independent.jsonl'],folder/'independent.log')
        cmp=compare(folder)
        # B-02 comparator reuses its implementation; this wrapper labels scope.
        full=read(folder/'comparison.json');full.pop('b02_complete',None);full['b03_complete']=False;write(folder/'comparison.json',full)
        cmp.pop('b02_complete',None)
        plain=command([frozen,'plain',*base,folder/'plain'],folder/'plain.log')
        original=read(folder/'native/run.json');assert original['complete'] and exact(original,read(folder/'plain/run.json'))
        obs=[json.loads(s) for s in (folder/'native/accepted.jsonl').read_text().splitlines()]
        for stop in stops:
            checkpoint=folder/f'native/checkpoint-{stop}.json';target=folder/f'resumed-{stop}'
            resumed=command([frozen,'resume',*base,target,checkpoint],folder/f'resumed-{stop}.log')
            later=[json.loads(s) for s in (target/'accepted.jsonl').read_text().splitlines()]
            row=dict(case=name,checkpoint_accepted_index=stop,producer_pid=native['pid'],resume_pid=resumed['pid'],
                     checkpoint_sha256=sha(checkpoint),later_accepted_states=len(later),
                     run_exact_except_elapsed=exact(original,read(target/'run.json')),
                     all_later_states_exact=exact(later,[s for s in obs if s['accepted_index']>stop]))
            row['passed']=row['run_exact_except_elapsed'] and row['all_later_states_exact'] and native['pid']!=resumed['pid']
            restarts.append(row);write(out/'restarts.json',restarts);assert row['passed']
        rec=dict(case=name,kind=kind,level=level,load_steps=steps,nodes=len(m['nodes']),tet10=len(m['tetrahedra']),
                 plates=len(c['beam_system']['plates']),native_process=native,independent_process=independent,
                 plain_process=plain,solver_sha256=sha(frozen),input_sha256={p.name:sha(p) for p in base},
                 all_run_equal_to_unobserved=True,comparison=cmp)
        records.append(rec);write(out/'measurements.json',dict(passed=False,b03_complete=False,cases=records))
        print(name,cmp['accepted_points'],'existing fields',cmp['existing_field_gates_passed'],
              'reassembly',cmp['maximum_native_reassembled_residual'],cmp['maximum_independent_reassembled_residual'],flush=True)
    after={str(p.relative_to(ROOT)):sha(p) for p in sources};write(out/'source-binding-after.json',after);assert before==after
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);p.add_argument('--binary',type=Path,required=True)
    p.add_argument('--kinds',nargs='+',default=['kink-zero','normal-sign-pair','step-split','kink-split','height-step']);p.add_argument('--level-max',type=int,default=1)
    a=p.parse_args();run(a.output.resolve(),a.binary.resolve(),a.kinds,a.level_max)
