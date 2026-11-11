"""Fresh bounded matrix. Each child and failure is retained in a new folder."""
import argparse
import importlib.util
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('c04_orchestrator_study', HERE/'study.py')
s = importlib.util.module_from_spec(spec); spec.loader.exec_module(s)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--probe', type=Path, default=s.ROOT/'target/release/examples/wet_joint_transient_probe')
    p.add_argument('--solver', type=Path, default=s.ROOT/'target/release/soil-fem-3d')
    p.add_argument('--element', choices=['tet4', 'tet10'], default='tet4')
    p.add_argument('--mechanical-steps', type=int, default=32)
    p.add_argument('--water-steps', nargs='+', type=int, default=[8, 16, 32])
    p.add_argument('--hydraulic-steps', nargs='+', type=int, default=[4, 8, 16])
    p.add_argument('--timeout', type=int, default=900)
    a = p.parse_args(); out = a.output.resolve(); start = time.monotonic()
    assert len(a.water_steps) >= 3 and len(a.hydraulic_steps) >= 3
    assert all(n >= 2 and n & (n-1) == 0 for n in [a.mechanical_steps]+a.water_steps+a.hydraulic_steps)
    s.initialize(out, a.probe.resolve(), a.solver.resolve(), a.element, a.mechanical_steps)
    calls = []
    def call(mode, folder, extra=()):
        command = [sys.executable, HERE/'study.py', '--mode', mode, '--output', out, '--folder', folder,
                   '--mechanical-steps', a.mechanical_steps, '--timeout', a.timeout, *extra]
        proc = s.invoke(command, out/(folder+'-driver.log'), a.timeout)
        calls.append(proc)
        if proc['exit_code']: raise RuntimeError(('matrix child failed', mode, folder, proc))
        print(mode, folder, round(proc['elapsed_seconds'], 3), flush=True)
    try:
        for n in a.water_steps:
            for kind in ('native', 'independent'):
                call(kind, f'{kind}-w{n}', ['--water-steps', n])
        for n in a.hydraulic_steps:
            call('hydraulic', f'hydraulic-h{n}', ['--hydraulic-steps', n])
        call('refusals', 'scope-refusals')
        call('restart', 'mechanical-restarts', ['--source', f'native-w{a.water_steps[-1]}'])
        call('hydraulic-restart', 'hydraulic-restarts', ['--source', f'hydraulic-h{a.hydraulic_steps[-1]}'])
        call('whole', 'whole', ['--water-steps', a.water_steps[-1]])
        proc = s.invoke([sys.executable, HERE/'measure.py', '--output', out, '--report', out/'measurements.json',
                         '--water-steps', *a.water_steps, '--hydraulic-steps', *a.hydraulic_steps], out/'measurement.log', a.timeout)
        calls.append(proc)
        if proc['exit_code']: raise RuntimeError(('measurement incomplete', proc))
        s.write(out/'completion.json', dict(complete=True, fresh_mechanical_3d_paths=2*len(a.water_steps),
                fresh_hydraulic_paths=2*len(a.hydraulic_steps), restarts_are_not_fresh_independent_paths=True,
                C04='未達', passed=False, A04_complete=False, all_processes=calls, elapsed_seconds=time.monotonic()-start))
    except Exception as e:
        s.write(out/'matrix-failure.json', dict(complete=False, passed=False, error=repr(e), processes=calls, elapsed_seconds=time.monotonic()-start))
        raise


if __name__ == '__main__': main()
