"""Fresh 768/8,16 full solves; observe the five saved physical frame endpoints.

This avoids quadratic checkpoint replay observation cost. The full solve still
accepts every original increment; intermediate field observations are not invented.
"""
import argparse
import json
import shutil
import sys
from pathlib import Path
from common import ROOT,HERE,inputs,write,read,command,sha
from compare import compare
from study import exact

def run(out,cli,observer):
    out.mkdir(parents=True,exist_ok=False);solver=out/'solver';solver.mkdir()
    frozen=solver/cli.name;helper=solver/observer.name
    shutil.copy2(cli,frozen);shutil.copy2(observer,helper)
    sources=[p for tree in ('src','tools','tests','examples') for p in (ROOT/tree).rglob('*') if p.suffix in ('.py','.rs','.json','.toml')]
    sources += [p for folder in (HERE,ROOT/'interim/2026-12-01-tri6-offset-composite-b02') for p in folder.iterdir() if p.suffix in ('.py','.rs')]
    before={str(p.relative_to(ROOT)):sha(p) for p in sources};write(out/'source-binding-before.json',before)
    records=[]
    for steps in (8,16):
        name=f'mesh2-load{steps}';folder=out/name;folder.mkdir()
        m,g,c,refinements=inputs(2,steps,.05,True)
        for key,value in [('mesh',m),('ground',g),('stages',c),('refinement',refinements)]:write(folder/(key+'.json'),value)
        paths=[folder/(s+'.json') for s in ('mesh','ground','stages')];(folder/'plain').mkdir()
        native=command([frozen,'solve-stages',*paths,folder/'plain/run.json'],folder/'native.log')
        assert read(folder/'plain/run.json')['complete']
        observation=command([helper,'inspect',*paths,folder/'native',folder/'plain/run.json'],folder/'observation.log')
        assert exact(read(folder/'native/run.json'),read(folder/'plain/run.json'))
        independent=command([sys.executable,HERE/'reference_probe.py',folder,folder/'independent-all.jsonl'],folder/'independent.log')
        with (folder/'independent-all.jsonl').open() as source,(folder/'independent.jsonl').open('x') as selected:
            all_count=0;selected_count=0
            for line in source:
                all_count+=1
                if json.loads(line)['progress']==1.:selected.write(line);selected_count+=1
        result=compare(folder)
        records.append(dict(name=name,level=2,steps=steps,tet10=len(m['tetrahedra']),native=native,independent=independent,observation=observation,
                            solver_sha256=sha(frozen),observer_sha256=sha(helper),input_sha256={p.name:sha(p) for p in paths},unobserved_run_equal=True,comparison=result,
                            native_field_observation='five saved frame endpoints; intermediate fields not observed',native_accepted_increments=sum(a['accepted'] for a in read(folder/'native/run.json')['attempts'])+sum(a['accepted'] for a in read(folder/'native/run.json')['geostatic']['attempts']),
                            independent_accepted_increments=all_count,independent_comparison_observations=selected_count))
        write(out/'fresh-solves.json',dict(passed=False,b02_complete=False,cases=records))
        print(name,'full accepted',all_count,'observed',selected_count,'native',result['maximum_native_reassembled_residual'],'independent',result['maximum_independent_reassembled_residual'],flush=True)
    after={str(p.relative_to(ROOT)):sha(p) for p in sources};write(out/'source-binding-after.json',after);assert before==after

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);p.add_argument('--cli',type=Path,required=True);p.add_argument('--observer',type=Path,required=True);a=p.parse_args();run(a.output.resolve(),a.cli.resolve(),a.observer.resolve())
