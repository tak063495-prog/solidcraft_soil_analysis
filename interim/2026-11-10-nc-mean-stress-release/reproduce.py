"""Recompute the bounded 5 kPa diagnostics, preserving each method's own history."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def commands(w):
    p = w/'prepared'
    binary = p/'observer-target/release/examples'/('material_diagnostic.exe' if sys.platform == 'win32' else 'material_diagnostic')
    common = ['--evidence', str(p/'origin'), '--binary', str(binary), '--reference', str(p/'reference400')]
    result = []
    for script, suffix in [('study.py', ''), ('regrowth_study.py', '-regrowth'), ('rate_study.py', '-rate')]:
        for kind in (['native', 'independent'] if not suffix else ['independent']):
            for order in ('mean', 'complement'):
                for steps in (128, 256):
                    label = f'{kind}-{order}{steps}{suffix}'
                    budget = (600 if steps == 128 else 900) if not suffix else 900
                    result.append((label, [sys.executable, str(HERE/script), '--map', kind, '--order', order+'-first',
                        *common, '--output', str(w/label), '--step-size', str(1/steps), '--budget-seconds', str(budget), '--count', '4096']))
    for kind,script in [('native','bounded_study.py'),('independent','bounded_rate_study.py')]:
        for order in ('mean','complement'):
            for steps in (128,256):
                label=f'bounded-{kind}-{order}{steps}'
                result.append((label,[sys.executable,str(HERE/script),'--map',kind,'--order',order+'-first',*common,
                    '--output',str(w/label),'--step-size',str(1/steps),'--budget-seconds','900','--count','4096']))
    for label, reference, material in [('independent-complement128-rate-m8', 'reference400', 8), ('independent-complement128-rate-pre800', 'reference800', 4)]:
        result.append((label, [sys.executable, str(HERE/'rate_study.py'), '--map', 'independent', '--order', 'complement-first',
            '--evidence', str(p/'origin'), '--binary', str(binary), '--reference', str(p/reference), '--output', str(w/label),
            '--step-size', str(1/128), '--material-steps', str(material), '--budget-seconds', '900', '--count', '4096']))
    for label,reference,material in [('bounded-independent-complement128-m8','reference400',8),('bounded-independent-complement128-pre800','reference800',4)]:
        result.append((label,[sys.executable,str(HERE/'bounded_rate_study.py'),'--map','independent','--order','complement-first',
            '--evidence',str(p/'origin'),'--binary',str(binary),'--reference',str(p/reference),'--output',str(w/label),
            '--step-size',str(1/128),'--material-steps',str(material),'--budget-seconds','900','--count','4096']))
    return result


def run(w, jobs):
    w.mkdir(parents=True, exist_ok=False)
    records = []
    def invoke(label, command, timeout=None):
        start = time.monotonic(); external_timeout = False
        with (w/(label+'.log')).open('x') as log:
            child = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
            try:
                status = child.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                child.terminate(); child.wait(); status = child.returncode; external_timeout = True
                folder = w/label
                if (folder/'controls.json').exists():
                    (folder/'external-watchdog.json').write_text(json.dumps(dict(passed=False, reason='external watchdog: inner material call exceeded driver wall budget; incomplete, not collapse or acceptance',
                        observed_process_elapsed_seconds=time.monotonic()-start, budget_seconds=timeout-15, signal='terminate', solver_report_present=(folder/'report.json').exists()), indent=2)+'\n')
        return dict(label=label, command=command, exit_code=status, external_timeout=external_timeout, elapsed_seconds=time.monotonic()-start)
    prep = invoke('prepare', [sys.executable, str(HERE/'prepare.py'), '--output', str(w/'prepared')])
    records.append(prep)
    if prep['exit_code']:
        raise ValueError('preparation failed; do not substitute native fields for independent history')
    reference = invoke('reference800', [sys.executable, str(HERE.parent/'2026-11-09-nc-pseudo-arclength/capture_reference_history.py'), '--output', str(w/'prepared/reference800'), '--steps', '800'])
    records.append(reference)
    if reference['exit_code'] or not (w/'prepared/reference800/unload.json').exists():
        raise ValueError('independent refined prehistory did not reach unload')
    specs = commands(w)
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = [pool.submit(invoke, label, cmd, float(cmd[cmd.index('--budget-seconds')+1])+15) for label, cmd in specs]
        for future in futures:
            row = future.result(); records.append(row)
            if row['exit_code'] and not row['external_timeout']:
                raise ValueError('diagnostic process failed: '+row['label'])
    lookup = dict(specs)
    for original, label, seed in [('native-complement128', 'native-complement-resumed', 'early.json'),
                                  ('native-complement256', 'native-complement256-resumed', 'early.json'),
                                  ('independent-complement128-rate', 'independent-complement-rate-resumed', 'early.json'),
                                  ('independent-complement128', 'independent-complement-tail', 'restart.json')]:
        checkpoint = w/original/seed
        if not checkpoint.exists() or (seed=='restart.json' and json.loads(checkpoint.read_text())['payload']['coordinate']==1.):
            continue
        cmd = list(lookup[original]); cmd[cmd.index('--output')+1] = str(w/label); cmd[cmd.index('--budget-seconds')+1] = '900'
        cmd += ['--resume', str(checkpoint)]
        records.append(invoke(label, cmd, 915))
    for kind,script in [('native','bounded_study.py'),('independent','bounded_rate_study.py')]:
        original=f'bounded-{kind}-complement128';label=f'bounded-{kind}-complement-resumed';checkpoint=w/original/'early.json'
        if checkpoint.exists():
            cmd=list(lookup[original]);cmd[cmd.index('--output')+1]=str(w/label);cmd[cmd.index('--budget-seconds')+1]='300';cmd[cmd.index('--count')+1]='24';cmd+=['--resume',str(checkpoint)]
            records.append(invoke(label,cmd,315))
    binary = w/'prepared/observer-target/release/examples'/('material_diagnostic.exe' if sys.platform == 'win32' else 'material_diagnostic')
    for label, _ in specs:
        if label not in ('native-complement128', 'native-complement256', 'independent-complement128-rate', 'independent-complement256-rate'):
            continue
        report = w/label/'report.json'
        if not report.exists() or not json.loads(report.read_text())['full_removal_reached']:
            continue
        records.append(invoke(label+'-post', [sys.executable, str(HERE/'check_post_release.py'), '--completed', str(w/label),
            '--evidence', str(w/'prepared/origin'), '--binary', str(binary), '--reference', str(w/'prepared/reference400'), '--output', str(w/(label+'-post'))]))
    solver = ROOT/'target/release'/('soil-fem-3d.exe' if sys.platform == 'win32' else 'soil-fem-3d')
    subprocess.run(['cargo', 'build', '--locked', '--release'], cwd=ROOT, check=True)
    records.append(invoke('fresh-native-prefix', [sys.executable, str(HERE/'fresh_prefix.py'), '--solver', str(solver), '--evidence', str(w/'prepared/origin'), '--output', str(w/'fresh-native-prefix')]))
    (w/'reproduction-processes.json').write_text(json.dumps(records, indent=2)+'\n')
    subprocess.run([sys.executable, str(HERE/'measure.py'), '--evidence', str(w), '--output', str(w/'measurements.json')], check=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True); p.add_argument('--jobs', type=int, choices=[1, 2, 4, 8], default=4)
    p.add_argument('--list-only', action='store_true'); a = p.parse_args()
    if a.list_only:
        print(json.dumps(commands(a.output.resolve()), indent=2))
    else:
        run(a.output.resolve(), a.jobs)
