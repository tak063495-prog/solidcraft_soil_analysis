"""Unmodified production CLI controls; identical Run, no provenance exclusions."""
import argparse
import shutil
from pathlib import Path
from b01_cases import read, write
from study import command, sha, exact


def run(study, out, cli):
    out.mkdir(parents=True, exist_ok=False)
    binary = out/'solver'/cli.name
    binary.parent.mkdir()
    shutil.copy2(cli, binary)
    records = []
    for row in read(study/'measurements.json')['cases']:
        folder = study/row['case']
        target = out/row['case']
        target.mkdir()
        base = [folder/(s+'.json') for s in ('mesh', 'ground', 'stages')]
        process = command([binary, 'solve-stages', *base, target/'run.json'], target/'solve.log', expect=None)
        item = dict(case=row['case'], process=process, solver_sha256=sha(binary),
                    complete=False, exact_except_elapsed=False)
        if (target/'run.json').exists():
            item['complete'] = read(target/'run.json')['complete']
            item['exact_except_elapsed'] = exact(read(target/'run.json'), read(folder/'native/run.json'))
        if row['case'] == 'attached-n2-s4':
            checkpoint = command([binary, 'checkpoint-stages', *base, '10', target/'checkpoint-10.json'], target/'checkpoint.log')
            resumed = command([binary, 'resume-stages', *base, target/'checkpoint-10.json', target/'resumed.json'], target/'resume.log')
            item['cli_restart'] = dict(producer=checkpoint, resume=resumed,
                checkpoint_sha256=sha(target/'checkpoint-10.json'),
                different_pid=checkpoint['pid'] != resumed['pid'],
                run_exact_except_elapsed=exact(read(target/'run.json'), read(target/'resumed.json')))
        records.append(item)
        write(out/'cli-controls.json', dict(passed=False, b01_complete=False, records=records))
    write(out/'cli-controls.json', dict(passed=False, b01_complete=False, records=records,
        production_cli_matches_observed_runs=all(r['complete'] and r['exact_except_elapsed'] for r in records),
        cli_restart_matches=all(r.get('cli_restart', {}).get('run_exact_except_elapsed', True) for r in records)))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('study', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--cli', type=Path, required=True)
    a = p.parse_args()
    run(a.study.resolve(), a.output.resolve(), a.cli.resolve())
