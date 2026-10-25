"""Actual-CAD separate-process restarts and rejected converged trials."""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import shutil


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact(a, b):
    def stable(v):
        if isinstance(v, dict):return {k: stable(x) for k, x in v.items() if k != 'elapsed_seconds'}
        if isinstance(v, list):return [stable(x) for x in v]
        return v
    return stable(a) == stable(b)


def verify(base, output, reuse_restarts=None):
    root = Path(__file__).resolve().parents[1]
    solver = base/'solver/soil-fem-3d.exe'
    mesh, ground, case = [base/'inputs'/name for name in ('mesh.json', 'ground.json', 'case.json')]
    native_file = base/'terrain-job/result.json';native = read(native_file)
    assert native['complete'] and native['schema'] == 'soil-fem-projected-mcc-stage-result/1'
    inputs = [solver, mesh, ground, case, native_file, Path(__file__), root/'tools/run_logged.py']
    assert sha(solver) == native['provenance']['solver_binary_sha256']
    output.mkdir(exist_ok=False)

    def execute(tag, args, success=True):
        command = [sys.executable, '-X', 'utf8', str(root/'tools/run_logged.py'), str(output/(tag+'.log')), str(solver), *map(str, args)]
        process = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding='utf-8')
        if (process.returncode == 0) != success:raise AssertionError(tag+': '+process.stdout[-3000:])
        print(tag+(' completed' if success else ' refused as expected'), flush=True)

    geo = sum(a['accepted'] for a in native['geostatic']['attempts'])
    counts = [sum(a['accepted'] and a['stage'] == stage['name'] for a in native['attempts']) for stage in native['case']['stages']]
    assert len(counts) == 3 and all(n >= 2 for n in [geo]+counts)
    stops = [geo//2, geo+counts[0]//2, geo+counts[0]+counts[1]+counts[2]//2]
    if reuse_restarts:
        inputs.extend(reuse_restarts/name for stop in stops for name in
            (f'checkpoint-{stop}.json', f'resumed-{stop}.json', f'checkpoint-{stop}.log', f'resume-{stop}.log'))
    fingerprints = {str(p): sha(p) for p in inputs}
    restarts = []
    for stop in stops:
        checkpoint = output/f'checkpoint-{stop}.json';result = output/f'resumed-{stop}.json'
        if reuse_restarts:
            for name in (checkpoint.name, result.name, f'checkpoint-{stop}.log', f'resume-{stop}.log'):
                source=reuse_restarts/name
                if source.suffix=='.log':assert source.read_text(encoding='utf-8').rstrip().endswith('Exit code: 0')
                shutil.copy2(source, output/name)
        else:
            execute(f'checkpoint-{stop}', ['checkpoint-stages', mesh, ground, case, stop, checkpoint])
            execute(f'resume-{stop}', ['resume-stages', mesh, ground, case, checkpoint, result])
        restored = read(result)
        assert restored['complete'] and exact(native, restored), stop
        sealed = read(checkpoint)
        assert sealed['identity']==native['provenance']
        state = sealed['run']['last_committed']
        offsets = state['mixed_mcc']['projected_volume']['compatibility_offsets']
        restarts.append(dict(accepted_increments=stop, full_result_except_elapsed_seconds_exact=True,
            execution_evidence_origin=str(reuse_restarts.resolve()) if reuse_restarts else 'new separate-process executions in this verification',
            maximum_frozen_compatibility_offset_m3=max(abs(d['integral_m3']) for d in offsets)))
    control = copy.deepcopy(read(case));injected = copy.deepcopy(control)
    injected['verification_reject_once_stages'] = [s['name'] for s in injected['stages']]
    for stage in control['stages']:stage['initial_increment'] = control['settings']['initial_increment']/2
    # The first 1/16 reload on this overconsolidated slope is elastic. A larger
    # converged trial is required to exercise discarded plastic history.
    injected['stages'][0]['initial_increment']=1.0
    control['stages'][0]['initial_increment']=0.5
    for path, value in ((output/'control-case.json', control), (output/'injected-case.json', injected)):
        with path.open('x', encoding='utf-8') as stream:json.dump(value, stream, indent=2)
    execute('control', ['solve-stages', mesh, ground, output/'control-case.json', output/'control.json'])
    execute('injected', ['solve-stages', mesh, ground, output/'injected-case.json', output/'injected.json'])
    a, b = read(output/'control.json'), read(output/'injected.json')
    assert a['complete'] and b['complete'] and exact(a['frames'], b['frames']) and exact(a['last_committed'], b['last_committed'])
    rejected = [row for row in b['attempts'] if not row['accepted']]
    assert len(rejected) == 3 and all(row['committed_state_unchanged_on_rejection'] for row in rejected)
    assert rejected[0]['rejected_trial_plastic_changed'] and rejected[-1]['rejected_trial_active_set_changed']
    changed = copy.deepcopy(read(case));changed['stages'][0]['target_factor'] += .1
    with (output/'changed-load-case.json').open('x', encoding='utf-8') as stream:json.dump(changed, stream, indent=2)
    execute('changed-load-refused', ['resume-stages', mesh, ground, output/'changed-load-case.json', output/f'checkpoint-{stops[-1]}.json', output/'must-not-exist.json'], False)
    assert not (output/'must-not-exist.json').exists()
    assert fingerprints == {str(p): sha(p) for p in inputs}
    report = dict(schema='soil-fem-projected-mcc-lifecycle-verification/1', passed=True,
        goal_complete=False, p6_complete=False, production_backend_accepted=False,
        scope='Actual CAD slope projected-volume MCC: geostatic/load/local-DEATH midpoint restart, pressure/volume/offset/material histories and rejected converged trials',
        restarts=restarts, rejected_converged_trials=3, rejected_plastic_and_domain_changes=True,
        full_control_frames_and_states_exact=True, changed_load_checkpoint_refused=True,
        input_and_source_sha256=fingerprints,
        evidence_sha256={p.name: sha(p) for p in output.iterdir() if p.is_file()})
    with (output/'verification.json').open('x', encoding='utf-8') as stream:json.dump(report, stream, indent=2)
    print(json.dumps(dict(passed=True, restarts=restarts, rejected_converged_trials=3)), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser();p.add_argument('base', type=Path);p.add_argument('output', type=Path)
    p.add_argument('--reuse-restarts',type=Path)
    a = p.parse_args();verify(a.base.resolve(), a.output.resolve(), a.reuse_restarts.resolve() if a.reuse_restarts else None)
