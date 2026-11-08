"""Bounded controls and one minimum trial; never produces a resumable Run."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import time
import zipfile
from prepare_experiment import prepare, BASE, ROOT, PRIOR


def read(p): return json.loads(p.read_text(encoding='utf-8'))
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p, v):
    with p.open('x', encoding='utf-8') as f:
        json.dump(v, f, indent=2, ensure_ascii=False, allow_nan=False); f.write('\n')
def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module

def matrices(path, label):
    return [dict(label=label+':'+str(e['iteration']), matrix_columns=e['matrix_columns'],
                 reduced_residual=e['reduced_residual'])
            for e in map(json.loads, path.read_text().splitlines()) if e['kind'] == 'tangent']

def stripped(v):
    if isinstance(v, dict): return {k:stripped(x) for k,x in v.items() if k != 'elapsed_seconds'}
    if isinstance(v, list): return [stripped(x) for x in v]
    return v

def capture(out):
    out.mkdir(parents=True, exist_ok=False)
    subprocess.run(['git','merge-base','--is-ancestor',BASE,'HEAD'], cwd=ROOT, check=True)
    manifest = read(PRIOR/'archive-sha256.json')
    assert sha(PRIOR/'evidence.zip') == manifest['archive_sha256']
    with zipfile.ZipFile(PRIOR/'evidence.zip') as z:
        assert set(z.namelist()) == set(manifest['entries'])
        for name, digest in manifest['entries'].items():
            p = Path(name); assert not p.is_absolute() and '..' not in p.parts
            data = z.read(name); assert hashlib.sha256(data).hexdigest() == digest
            dest = out/'prior'/p; dest.parent.mkdir(parents=True, exist_ok=True); dest.write_bytes(data)
    write(out/'prior-binding.json', manifest)
    records = []
    def run(args, name, timeout=240):
        args = [str(a) for a in args]; start = time.monotonic()
        with (out/(name+'.log')).open('x') as f:
            result = subprocess.run(args, cwd=ROOT, stdout=f, stderr=subprocess.STDOUT, timeout=timeout)
        records.append(dict(command=args, exit_code=result.returncode, elapsed_seconds=time.monotonic()-start))
        (out/'processes.json').write_text(json.dumps(records, indent=2)+'\n')
        result.check_returncode()
    target = out/'target'
    ext = '.exe' if os.name == 'nt' else ''
    binaries = {}
    for label, retain in [('control', False), ('retain', True)]:
        source = out/label; prepare(source, retain)
        run(['cargo','fmt','--manifest-path',source/'Cargo.toml'], label+'-fmt')
        run(['cargo','build','--locked','--release','--manifest-path',source/'Cargo.toml',
             '--target-dir',target,'--example','nc_stop_probe','--example','linear_probe'], label+'-build')
        binary = out/('nc_stop_probe-'+label+ext)
        shutil.copyfile(target/('release/examples/nc_stop_probe'+ext), binary); binary.chmod(0o755)
        binaries[label] = binary
    linear = target/('release/examples/linear_probe'+ext)
    write(out/'binary-sha256.json', {**{k:sha(v) for k,v in binaries.items()}, 'linear_probe':sha(linear)})
    req = []
    for name in ('last-accepted','terminal-rejected'):
        req += matrices(out/'prior/evidence/diagnostic'/name/'trace.jsonl', name)
    write(out/'gate-matrices.json', req)
    run([linear,out/'gate-matrices.json',out/'gate-linear-results.json'], 'gate-linear-probe')
    audit = load('prior_nc_audit', PRIOR/'audit_linear.py').audit
    gate = audit(req, read(out/'gate-linear-results.json')); write(out/'gate-linear-audit.json', gate)
    if not gate['all_saved_matrices_pass_independent_linear_gate']:
        write(out/'bounded-result.json', dict(passed=False, history_converged=False, nonlinear_trial_performed=False))
        return
    old = out/'prior/evidence/prior/nc-resumed'
    mesh = old/'input/mesh.json'
    for label, variant, cp, progress in [('control-trial','control',0,'0.9514002799987793'),
                                         ('retain-accepted-control','retain',1,'0.9514000415802002'),
                                         ('retain-trial','retain',0,'0.9514002799987793')]:
        run([binaries[variant],mesh,old/'job'/f'recovery-checkpoint-{cp}.json',progress,out/label],label)
    controls = dict(
        compensated_control_trace_exact_except_elapsed=stripped([json.loads(l) for l in (out/'control-trial/trace.jsonl').read_text().splitlines()]) ==
            stripped([json.loads(l) for l in (out/'prior/candidate-trial/trace.jsonl').read_text().splitlines()]),
        preceding_accepted_state_exact=read(out/'retain-accepted-control/diagnostic.json')['outcome'].get('state') ==
            read(old/'job/recovery-checkpoint-0.json')['run']['last_committed'])
    write(out/'controls.json',controls)
    req = matrices(out/'retain-trial/trace.jsonl','retain-trial'); write(out/'trial-matrices.json',req)
    run([linear,out/'trial-matrices.json',out/'trial-linear-results.json'],'trial-linear-probe')
    write(out/'trial-linear-audit.json',audit(req,read(out/'trial-linear-results.json')))
    outcome = read(out/'retain-trial/diagnostic.json')
    write(out/'bounded-result.json', dict(passed=False, history_converged=False, production_adopted=False,
        full_removal_reached=False, target_trial_converged=outcome['outcome']['converged'],
        checkpoint_resume=False, new_continuation_created=False, accepted_new_increments=0,
        controls=controls, note='single trial diagnostic; a converged trial is not a committed increment'))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();capture(args.output.resolve())
