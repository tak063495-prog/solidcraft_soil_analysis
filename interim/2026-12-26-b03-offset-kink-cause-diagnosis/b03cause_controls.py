"""Run as its own process: reuse original fixtures without leaking legacy imports."""
import argparse
import copy
import hashlib
import itertools
import json
import math
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def run(output, probe, solver):
    output.mkdir(parents=True, exist_ok=False)
    sys.path.insert(0, str(ROOT/'interim/2026-12-02-tri6-offset-kink-b03'))
    from b03_cases import make, detached_fold
    requests = []; requests_meta = []
    tau = 64*sys.float_info.epsilon*.1

    def model(name, data, expected):
        m, _, c, _ = data
        requests.append(dict(id=name, kind='model_guard', mesh=m, configuration=c['beam_system']))
        requests_meta.append(dict(id=name, expected_acceptance=expected))

    for angle in (15, 45, 90, 135, 179):
        for label, e1, e2 in (('positive', .05, .05), ('negative', -.05, -.05), ('one-zero', .05, 0.)):
            model(f'fold-{angle}-{label}', detached_fold(angle, e1, e2), False)
    for multiple in (.5, 1., 2.):
        data = detached_fold(0., 0., multiple*tau)
        model('near-step-'+str(multiple), data, multiple <= 1.)
    for sign in (-1, 1):
        for e in (math.nextafter(.05, 0.), .05, math.nextafter(.05, math.inf), .051):
            model('range-'+str(sign*e), detached_fold(0., sign*e, sign*e), abs(e) <= .05)
    # Same three panels, only declaration order changes; compare all permutations.
    for order in itertools.permutations(range(3)):
        data = detached_fold(0., .02, .02+.75*tau)
        cfg = data[2]['beam_system']; third = copy.deepcopy(cfg['plates'][0])
        # All three panels share only one edge. Reusing a complete first panel
        # would add another pairwise comparison at its private midside nodes.
        third.update(id='third', attachment_offset_m=.02+1.5*tau,
                     nodes=['v1', 'v0', 'third-v', 'm0-1', 'third-m02', 'third-m21'])
        for name, point in [('third-v', [0., -2., 1.]),
                            ('third-m02', [0., -1., 1.]), ('third-m21', [.5, -1., 1.])]:
            cfg['nodes'].append(dict(id=name, position_m=point, match_soil_node=False))
        cfg['plates'].append(third)
        cfg['plates'] = [cfg['plates'][i] for i in order]
        model('three-order-'+''.join(map(str, order)), data, None)
    for reverse in (False, True):
        data = detached_fold(0., .02, .02+1.5*tau)
        cfg = data[2]['beam_system']; other = copy.deepcopy(cfg['plate_sections'][0])
        other.update(id='thick', thickness_m=.2); cfg['plate_sections'].append(other)
        cfg['plates'][1]['section_id'] = 'thick'
        if reverse:
            cfg['plates'].reverse()
        model('thickness-order-'+str(reverse), data, not reverse)
    for kind in ('planar', 'kink-zero', 'normal-sign-pair', 'step-split', 'kink-split', 'height-step'):
        model('existing-'+kind, make(kind), True)
    write(output/'model-guard-requests.json', requests)
    p = subprocess.run([str(probe)], input=json.dumps(requests).encode(), capture_output=True)
    (output/'model-guard-responses.json').write_bytes(p.stdout)
    (output/'model-guard-stderr.log').write_bytes(p.stderr)
    p.check_returncode(); results = json.loads(p.stdout)
    for meta, result in zip(requests_meta, results):
        meta.update(result)
        expected = meta['expected_acceptance']
        meta['matches_predeclared_expectation'] = None if expected is None else expected == result['production_model_guard']['accepted']
    rows = []

    def cli(name, data, expected=None):
        folder = output/'cli'/name; folder.mkdir(parents=True)
        for key, value in zip(('mesh', 'ground', 'stages'), data[:3]):
            write(folder/(key+'.json'), value)
        cmd = [str(solver), 'solve-stages', *[str(folder/(k+'.json')) for k in ('mesh', 'ground', 'stages')], str(folder/'run.json')]
        started = time.monotonic()
        with (folder/'cli.log').open('w', encoding='utf-8') as log:
            proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT); code = proc.wait()
        text = (folder/'cli.log').read_text(); exists = (folder/'run.json').exists()
        complete = json.loads((folder/'run.json').read_text()).get('complete') if exists else None
        rows.append(dict(name=name, expected_diagnostic=expected, command=cmd, pid=proc.pid,
             exit_code=code, elapsed_seconds=time.monotonic()-started, result_exists=exists, complete=complete,
             diagnostic=text, input_sha256={k:hashlib.sha256((folder/(k+'.json')).read_bytes()).hexdigest() for k in ('mesh','ground','stages')},
             contract_confirmed=bool(code==0 and complete is True) if expected is None else bool(code==1 and not exists and expected in text)))

    for kind in ('planar', 'kink-zero', 'normal-sign-pair', 'step-split', 'kink-split', 'height-step'):
        cli(kind, make(kind))
    for name, data, diagnostic in (
        ('nonzero-fold', detached_fold(90, .05, .05), 'same signed normal offset'),
        ('one-zero-fold', detached_fold(90, .05, 0.), 'same signed normal offset'),
        ('e-step', detached_fold(0, .05, .04), 'same signed normal offset'),
        ('e-only', detached_fold(0, .05, -.05), 'same signed normal offset'),
        ('range-positive', detached_fold(0, .051, .051), 'finite |e| <= thickness/2'),
        ('range-negative', detached_fold(0, -.051, -.051), 'finite |e| <= thickness/2')):
        cli(name, data, diagnostic)
    data = detached_fold(0, .05, .05)
    plate = data[2]['beam_system']['plates'][1]; plate['nodes'] = [plate['nodes'][j] for j in (0, 2, 1, 5, 4, 3)]
    cli('normal-only', data, 'same signed normal offset')
    for axis in (0, 2):
        data = detached_fold(0, 0, 0); cfg = data[2]['beam_system']
        node = next(n for n in cfg['nodes'] if n['id'] == cfg['plates'][0]['nodes'][3]); node['position_m'][axis] += .001
        cli('curved-midpoint-'+str(axis), data, 'Mindlin Tri6 requires straight edges')
    for name, key, value in (('unknown-normal', 'normal', [0, 0, 1]), ('unknown-e', 'attachment_offset_m', .05)):
        data = make('planar'); data[2]['beam_system']['nodes'][0][key] = value
        cli(name, data, 'unknown field')
    data = make('planar'); cfg = data[2]['beam_system']; node = next(n for n in cfg['nodes'] if n['id']==cfg['plates'][0]['nodes'][3])
    node['position_m'] = node['position_m'].copy(); node['position_m'][0] += .001
    cli('attachment-mismatch', data, 'attachment coordinate mismatch')
    # Existing headless preparation refusal script also checks legacy default inputs.
    cmd = [sys.executable, str(ROOT/'interim/2026-12-02-tri6-offset-kink-b03/supplemental.py'),
           str(output/'supplemental'), '--cli', str(solver)]
    with (output/'supplemental-driver.log').open('w', encoding='utf-8') as log:
        status = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    status.check_returncode()
    supplement = json.loads((output/'supplemental/supplemental.json').read_text())
    report = dict(passed=False, b03_complete=False, diagnostic_only=True, production_scope_expanded=False,
        solver_sha256=hashlib.sha256(solver.read_bytes()).hexdigest(),
        model_guards=requests_meta, model_guard_requests=len(requests),
        cli_controls=rows, preparation_and_legacy=supplement,
        predeclared_model_expectations_confirmed=all(r['matches_predeclared_expectation'] is not False for r in requests_meta),
        all_cli_contracts_confirmed=all(r['contract_confirmed'] for r in rows),
        gui_gpu_real_cad_native_display_executed=False)
    write(output/'summary.json', report)
    print(json.dumps({k:v for k,v in report.items() if k not in ('model_guards','cli_controls','preparation_and_legacy')}))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--probe', type=Path, required=True)
    p.add_argument('--solver', type=Path, required=True)
    a = p.parse_args(); run(a.output.resolve(), a.probe.resolve(), a.solver.resolve())
