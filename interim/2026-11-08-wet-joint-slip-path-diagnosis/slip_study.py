"""Bounded input-only diagnostics of the EXISTING 54-node wet JOINT cycle.

All-history acceptance is deliberately false. No thresholds are fitted here.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import time
import traceback
import zipfile
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from diagnose_cut_volume_history import inputs, validate_path, summary, FIELDS
from wet_joint_cycle_reference import serializable, differences
from reference_observed import cut_volume_reference

BASE = 'd2b3dd4da9a1e92bbb6bc27b8d03acfc81f0903d'
OBSERVED = ('u', 'jump', 'traction', 'plastic_slip', 'accumulated')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    Path(path).write_text(json.dumps(serializable(value), ensure_ascii=False,
                                   allow_nan=False, indent=2)+'\n', encoding='utf-8')


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def encode(array):
    """Same reversible uint64 predictor as the existing accepted-field archive."""
    raw = np.ascontiguousarray(array, dtype='<f8').view('<u8')
    delta = raw.copy()
    delta[1:2] = raw[1:2] ^ raw[:1]
    delta[2:] = raw[2:] ^ (raw[1:-1]+(raw[1:-1]-raw[:-2]))
    return np.ascontiguousarray(np.moveaxis(delta.view('u1').reshape(*raw.shape, 8), -1, 0))


def decode(array):
    assert array.dtype == np.uint8 and array.shape[0] == 8
    shape = array.shape[1:]
    raw = np.ascontiguousarray(np.moveaxis(array, 0, -1)).view('<u8').reshape(shape)
    if len(raw) > 1:
        raw[1] ^= raw[0]
    for i in range(2, len(raw)):
        raw[i] ^= raw[i-1]+(raw[i-1]-raw[i-2])
    return raw.view('<f8')


class Sink:
    def __init__(self, folder):
        self.folder = folder
        folder.mkdir()
        self.rows = []
        self.values = []
        self.chunks = []
        self.count = 0

    def __call__(self, row, fields):
        if self.rows and row['stage'] != self.rows[-1]['stage']:
            self.flush()
        self.rows.append(row)
        self.values.append({k: fields[k].copy() for k in OBSERVED})
        if len(self.rows) >= 256:
            self.flush()

    def flush(self):
        if not self.rows:
            return
        path = self.folder/f'{len(self.chunks):04d}.npz'
        arrays = {k: np.stack([v[k] for v in self.values]) for k in OBSERVED}
        coded = {k: encode(v) for k, v in arrays.items()}
        for k, v in coded.items():
            assert np.array_equal(decode(v.copy()).view('<u8'), arrays[k].view('<u8'))
        np.savez_compressed(path, **coded)
        self.chunks.append(dict(file=path.name, sha256=sha(path), start=self.count,
                                count=len(self.rows), rows=self.rows))
        self.count += len(self.rows)
        write(self.folder/'manifest.json', dict(encoding='uint64-linear-xor-byte-shuffle',
              fields=list(OBSERVED), count=self.count, chunks=self.chunks))
        print('accepted', self.count, self.rows[-1]['stage'], self.rows[-1]['progress'], flush=True)
        self.rows = []
        self.values = []


def solve_reference(folder, steps, local):
    folder.mkdir(parents=True, exist_ok=False)
    values = inputs()
    for name, value in zip(('mesh', 'ground', 'stages'), values):
        write(folder/(name+'.json'), value)
    sink = Sink(folder/'accepted')
    linear = dict(count=0, maximum_true_relative_residual=0.)

    def check_linear(row):
        r = row['true_relative_residual']
        assert np.isfinite(r) and r <= 1e-10, ('true linear residual', row)
        linear['count'] += 1
        if r > linear['maximum_true_relative_residual']:
            linear.update(maximum_true_relative_residual=r, maximum_at=row)

    start = time.monotonic()
    try:
        result = cut_volume_reference(*values, steps, local, 7, sink, check_linear)
        sink.flush()
        validate_path(result, *values, steps, local)
        assert sink.count == steps*7
        write(folder/'result.json', result)
        status = dict(complete=True, steps=steps, contact_steps=local,
                      accepted=sink.count, linear=linear, elapsed_seconds=time.monotonic()-start)
    except Exception as exc:
        sink.flush()
        status = dict(complete=False, steps=steps, contact_steps=local, accepted=sink.count,
                      failure=repr(exc), linear=linear, elapsed_seconds=time.monotonic()-start)
        write(folder/'status.json', status)
        raise
    write(folder/'status.json', status)
    return status


def bounded_reference(folder, steps, local, seconds):
    """The outer wall clock bound does not change Newton/tolerances/step sizes."""
    assert not folder.exists()
    log = folder.parent/(folder.name+'.log')
    command = [sys.executable, str(Path(__file__).resolve()), 'reference',
               '--output', str(folder), '--steps', str(steps), '--local', str(local)]
    start = time.monotonic()
    with log.open('xb') as stream:
        child = subprocess.Popen(command, stdout=stream, stderr=subprocess.STDOUT)
        try:
            code = child.wait(timeout=seconds)
            reason = 'process_exited'
        except subprocess.TimeoutExpired:
            child.terminate()
            try:
                code = child.wait(timeout=10)
            except subprocess.TimeoutExpired:
                child.kill()
                code = child.wait()
            reason = 'wall_clock_limit'
    result = dict(command=command, returncode=code, reason=reason, budget_seconds=seconds,
                  elapsed_seconds=time.monotonic()-start,
                  complete=code == 0 and (folder/'result.json').exists(),
                  passed=False, history_converged=False)
    write(folder/'terminal.json', result)
    return result


def pair(a, b):
    """Seven endpoint comparisons with the actual denominator and point causes."""
    rows = differences(a, b)
    detail = []
    previous_delta = np.zeros((32, 3))
    for x, y, row in zip(a['frames'], b['frames'], rows):
        p = np.asarray(x['plastic_slip']); q = np.asarray(y['plastic_slip'])
        acc = np.asarray(y['accumulated']); delta = p-q
        p2 = np.sum(delta**2, axis=1)
        total = float(p2.sum())
        order = np.argsort(-p2)
        detail.append(dict(frame=row['frame'], numerator_l2_m=float(np.linalg.norm(delta)),
            denominator_l2_m=float(np.linalg.norm(q)), denominator_floor_m=1e-12,
            accumulated_l2_m=float(np.linalg.norm(acc)),
            net_over_accumulated=float(np.linalg.norm(q)/max(np.linalg.norm(acc), 1e-30)),
            difference_over_accumulated=float(np.linalg.norm(delta)/max(np.linalg.norm(acc), 1e-30)),
            error_change_since_previous_frame_l2_m=float(np.linalg.norm(delta-previous_delta)),
            summed_vector_a_m=p.sum(axis=0), summed_vector_b_m=q.sum(axis=0),
            maximum_point_vector_difference_m=float(np.sqrt(p2.max())),
            all_points=[dict(point=int(i), a_m=p[i], b_m=q[i], delta_m=delta[i],
                squared_difference_fraction=float(p2[i]/total) if total else 0.) for i in order]))
        previous_delta = delta
    return dict(fields=rows, maxima=summary(rows), plastic_slip_detail=detail,
                passed=False, history_converged=False)


def native(solver, folder):
    from verify_checkpoint import invoke
    from verify_wet_joint_cycle import restarts, path_evidence
    prod = load_module('cut_production_study', ROOT/'interim/2026-11-08-wet-joint-production-cut-volume/study.py')
    folder.mkdir(parents=True, exist_ok=False)
    # Exact original physical AND numerical input; no new geometry or NC.
    for name, file in zip(('mesh', 'ground', 'stages'),
            ('joint-water-column-tet10-mesh', 'ground-joint-water-mc', 'stages-joint-water-cycle')):
        (folder/(name+'.json')).write_bytes((ROOT/'examples'/(file+'.json')).read_bytes())
    start = time.monotonic()
    invoke([str(solver), 'solve-stages', *[str(folder/(n+'.json')) for n in ('mesh','ground','stages')],
            str(folder/'stage-results.json')], folder/'solve.log')
    run = json.loads((folder/'stage-results.json').read_text())
    gates = prod.native_gates(run)
    assert run['case']['settings']['relative_force_tolerance'] == 1e-8
    assert run['case']['settings']['relative_linear_tolerance'] == 1e-10
    rows = restarts(solver, folder, run)
    fields, audits = prod.fields_and_audit(*inputs()[:2], run, None)
    report = dict(complete=run['complete'], frames=path_evidence(run), accepted_gates=gates,
        restarts=rows, fields=fields, independent_force_audit=audits,
        holding_audit=prod.holding_audit(*inputs()[:2], run), solver_sha256=sha(solver),
        elapsed_seconds=time.monotonic()-start, numerical_route='unchanged original input',
        full_history_acceptance=False)
    write(folder/'checks.json', report)
    return report


def saved(folder):
    """Authenticate the old raw results; do not silently relabel comparisons."""
    prod = load_module('cut_production_saved', ROOT/'interim/2026-11-08-wet-joint-production-cut-volume/study.py')
    refs, provenance = prod.saved_references()
    old = ROOT/'interim/2026-11-08-wet-joint-production-cut-volume'
    manifest = json.loads((old/'archive-sha256.json').read_text())
    assert sha(old/'evidence.zip') == manifest['archive_sha256']
    with zipfile.ZipFile(old/'evidence.zip') as z:
        name = next(n for n in z.namelist() if n.endswith('g512-m2.5e-6/stage-results.json'))
        raw = z.read(name)
        assert hashlib.sha256(raw).hexdigest() == manifest['entries'][name]
        run = json.loads(raw)
    native_fields, _ = prod.fields_and_audit(*inputs()[:2], run, None)
    refs['native-g512-m2.5e-6'] = native_fields
    provenance['native512'] = dict(archive_sha256=manifest['archive_sha256'], entry=name,
        sha256=manifest['entries'][name], rerun=False)
    folder.mkdir(exist_ok=False)
    for key, value in refs.items():
        write(folder/(key+'.json'), value)
    write(folder/'provenance.json', provenance)
    return refs, provenance


def measure(folder):
    refs = {p.stem: json.loads(p.read_text()) for p in (folder/'saved').glob('*.json') if p.stem != 'provenance'}
    routes = {}
    for p in sorted((folder/'paths').glob('g*')):
        if not p.is_dir():
            continue
        terminal = json.loads((p/'terminal.json').read_text())
        manifest = json.loads((p/'accepted/manifest.json').read_text())
        terminal['flushed_accepted_states'] = manifest['count']
        terminal['last_flushed_accepted'] = manifest['chunks'][-1]['rows'][-1]
        terminal['unflushed_tail_certified'] = False
        routes[p.name] = terminal
        if terminal['complete']:
            result = json.loads((p/'result.json').read_text())
            validate_path(result, *inputs(), result['steps'], result['contact_steps'])
            terminal['solve_checks'] = json.loads((p/'status.json').read_text())
            terminal['maximum_relative_force_residual'] = result['maximum_relative_force_residual']
            refs[p.name] = result
        else:
            terminal['whole_path_linear_summary_available'] = False
    pairs = [('g2048-l64','g4096-l64'), ('native-g512-m2.5e-6','g4096-l64'),
             ('g512-l64','g4096-l64'), ('native-g512-m2.5e-6','g512-l64'),
             ('g512-l64','g1024-l64'), ('g1024-l64','g2048-l64'),
             ('g1024-l64','g4096-l64'),
             ('g512-l64','g512-l128'), ('g512-l128','g4096-l64'),
             ('g4096-l64','g4096-l128'), ('g4096-l64','g8192-l64'),
             ('g2048-l64','g8192-l64')]
    results = {a+'__'+b: pair(refs[a],refs[b]) for a,b in pairs if a in refs and b in refs}
    return dict(schema='wet-joint-slip-path-diagnosis/1', starting_commit=BASE,
        passed=False, history_converged=False, production_adopted=False,
        comparison_identity='5.083166% is native512(m2.5e-6) vs independent4096/64; 0.953507% is independent2048/64 vs4096/64. They are different pairs.',
        routes=routes, comparisons=results,
        provenance=json.loads((folder/'saved/provenance.json').read_text()),
        native=json.loads((folder/'native-original/checks.json').read_text()),
        all_history_acceptance_criterion_available=False,
        nc=dict(executed=False, retained_progress=0.9514000415802002, accepted_increments=266,
                next_minimum_increment_newton_limit=60, passed=False),
        large_runs_executed=False, windows_ci='unverified',
        p3b_p4_rows_complete=False, structure_water_joint_complete=False, p6_complete=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('reference','bounded','native','saved','measure'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--steps', type=int, default=512)
    parser.add_argument('--local', type=int, default=64)
    parser.add_argument('--seconds', type=float, default=600)
    parser.add_argument('--solver', type=Path)
    args = parser.parse_args()
    out = args.output.resolve()
    if args.action == 'reference':
        result = solve_reference(out, args.steps, args.local)
    elif args.action == 'bounded':
        result = bounded_reference(out, args.steps, args.local, args.seconds)
    elif args.action == 'native':
        result = native(args.solver.resolve(), out)
    elif args.action == 'saved':
        result = saved(out)[1]
    else:
        result = measure(out)
        write(out/'measurements.json', result)
    print(json.dumps({k:v for k,v in result.items() if k in ('complete','elapsed_seconds','passed','history_converged')}, indent=2))
