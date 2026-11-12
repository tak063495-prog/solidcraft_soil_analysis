"""Fresh solve orchestration; evidence is never silently replaced."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import copy
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
from cases import ROOT, inputs, read, write
from compare import compare
HERE = Path(__file__).resolve().parent


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def semantic(v):
    if isinstance(v, dict): return {k:semantic(x) for k,x in v.items() if k != 'elapsed_seconds'}
    if isinstance(v, list): return [semantic(x) for x in v]
    return v


def exact(a, b):
    return json.dumps(semantic(a), sort_keys=True, separators=(',', ':'), allow_nan=False) == json.dumps(semantic(b), sort_keys=True, separators=(',', ':'), allow_nan=False)


def command(cmd, log, expect=0):
    log.parent.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    with log.open('x', encoding='utf-8') as f:
        p = subprocess.Popen(list(map(str, cmd)), cwd=ROOT, stdout=f, stderr=subprocess.STDOUT, text=True)
        rc = p.wait()
    meta = dict(command=list(map(str,cmd)), pid=p.pid, exit_code=rc,
                elapsed_seconds=time.monotonic()-start, log_sha256=sha(log))
    write(log.with_suffix('.process.json'), meta)
    if expect is not None and rc != expect: raise RuntimeError(f'command rc={rc}: {log}')
    return meta


def run(output, binary, names=None):
    output.mkdir(parents=True, exist_ok=False)
    binary = binary.resolve()
    frozen = output/'solver'/binary.name; frozen.parent.mkdir(); shutil.copy2(binary, frozen)
    sources = [p for base in ('src','tools','tests','examples') for p in (ROOT/base).rglob('*') if p.suffix in ('.py','.rs','.json','.toml')]
    sources += [p for p in HERE.iterdir() if p.suffix in ('.py','.rs','.md')]
    binding = {str(p.relative_to(ROOT)):sha(p) for p in sources}
    write(output/'source-binding-before.json', binding)
    plan = [(f'{"beam" if b else "plate"}-{sign}',0,4,e,b) for b in (False,True) for sign,e in (('negative',-.05),('zero',0.),('positive',.05))]
    plan += [('mesh-1',1,4,.05,True), ('mesh-2',2,4,.05,True), ('load-8',0,8,.05,True), ('load-16',0,16,.05,True)]
    if names: plan = [p for p in plan if p[0] in names]
    records = []; restarts = []
    for name, level, steps, offset, beam in plan:
        folder = output/name; folder.mkdir()
        m,g,c,refinements = inputs(level,steps,offset,beam)
        for key,v in [('mesh',m),('ground',g),('stages',c),('refinement',refinements)]: write(folder/(key+'.json'),v)
        base = [folder/(s+'.json') for s in ('mesh','ground','stages')]
        stops = [2*steps+steps//2]
        if name in ('beam-negative','beam-positive'): stops += [steps+steps//2,3*steps+steps//2,4*steps+steps//2]
        native = command([frozen,'solve',*base,folder/'native',','.join(map(str,stops))], folder/'native.log')
        independent = command([sys.executable,HERE/'reference.py',folder,folder/'independent.jsonl'],folder/'independent.log')
        comparison = compare(folder)
        direct = command([frozen,'plain',*base,folder/'plain'],folder/'plain.log')
        original = read(folder/'native/run.json'); plain = read(folder/'plain/run.json')
        assert original['complete'] and exact(original,plain), name+' diagnostic stop/resume vs unobserved Run'
        observations = [json.loads(s) for s in (folder/'native/accepted.jsonl').read_text().splitlines()]
        for stop in stops:
            saved = folder/f'native/checkpoint-{stop}.json'
            target = folder/f'resumed-{stop}'
            resumed = command([frozen,'resume',*base,target,saved],folder/f'resumed-{stop}.log')
            rr = read(target/'run.json')
            later = [json.loads(s) for s in (target/'accepted.jsonl').read_text().splitlines()]
            expected = [s for s in observations if s['accepted_index']>stop]
            assert native['pid'] != resumed['pid']
            match = exact(rr,original) and exact(later,expected)
            item = dict(case=name, checkpoint_accepted_index=stop, producer_pid=native['pid'], resume_pid=resumed['pid'],
                        checkpoint_sha256=sha(saved), later_accepted_states=len(later), run_exact_except_elapsed=exact(rr,original),
                        all_later_states_exact=exact(later,expected), passed=match)
            restarts.append(item); write(output/'restarts.json', restarts)
            assert match, item
        records.append(dict(case=name, level=level, load_steps=steps, offset_m=offset, beam=beam,
                            nodes=len(m['nodes']), tet10=len(m['tetrahedra']), plates=len(c['beam_system']['plates']),
                            native_process=native, independent_process=independent, plain_process=direct,
                            solver_sha256=sha(frozen), input_sha256={s.name:sha(s) for s in base},
                            all_run_equal_to_unobserved=True, comparison=comparison))
        write(output/'measurements.json',dict(passed=False,b02_complete=False,cases=records,
              note='B-02 unsupported combinations remain unsupported; per-field existing gates separately reported'))
        print(name, comparison['accepted_points'], 'existing gates', comparison['existing_field_gates_passed'], flush=True)
    after = {str(p.relative_to(ROOT)):sha(p) for p in sources}
    write(output/'source-binding-after.json', after)
    assert binding == after
    return records


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);p.add_argument('--binary',type=Path,required=True);p.add_argument('--cases',nargs='+')
    a=p.parse_args();run(a.output.resolve(),a.binary,a.cases)
