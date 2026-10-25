"""Run an original-input, fully compatible independent dry MCC construction."""
import argparse
import hashlib
import json
from pathlib import Path
import mcc_compatible_sparse_reference as sparse
import mcc_rate_fem as rate
import mcc_rate_midpoint as midpoint
import fem_reference as ref
import verify_high_order
import high_order_io


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    for name in ('mesh', 'ground', 'case', 'identity', 'output'):
        parser.add_argument(name, type=Path)
    parser.add_argument('--steps', type=int, default=100)
    parser.add_argument('--dependency-lock', type=Path, required=True)
    parser.add_argument('--frame-directory', type=Path,
                        help='save accepted equilibrium frames; these are not mid-increment restart checkpoints')
    args = parser.parse_args()
    if args.steps < 2 or args.output.exists():
        parser.error('at least two steps and a new output required')
    lock = args.dependency_lock.resolve()
    runtime_lock = lock.parent/'requirements.lock.txt'
    pins = dict(line.split('==') for path in (lock, runtime_lock)
                for line in path.read_text(encoding='utf-8').splitlines() if '==' in line)
    if pins.get('scipy') != sparse.scipy.__version__ or pins.get('numpy') != sparse.np.__version__:
        parser.error('reference NumPy/SciPy versions differ from dependency locks')
    inputs = [getattr(args, key).resolve() for key in ('mesh', 'ground', 'case', 'identity')]
    sources = [Path(__file__).resolve(), lock, runtime_lock]+[Path(module.__file__).resolve()
              for module in (sparse, rate, midpoint, ref, verify_high_order, high_order_io)]
    bindings = {str(path): sha(path) for path in inputs+sources}
    frame_count = 0
    if args.frame_directory:
        args.frame_directory.mkdir(parents=True, exist_ok=False)

    def save_frame(frame):
        nonlocal frame_count
        if args.frame_directory is None:
            return
        assert bindings == {str(path): sha(path) for path in inputs+sources}, 'input or reference source changed'
        name = f'frame-{frame_count:03d}'
        temporary = args.frame_directory/(name+'.tmp')
        with temporary.open('x', encoding='utf-8') as stream:
            json.dump(dict(schema='soil-fem-compatible-reference-equilibrium-frame/1',
                           reference_steps=args.steps, frame=frame, input_sha256=bindings,
                           is_restart_checkpoint=False, goal_complete=False, p6_complete=False),
                      stream, allow_nan=False)
        temporary.rename(args.frame_directory/(name+'.json'))
        frame_count += 1

    mesh, ground, case, identity = map(read, inputs)
    result = rate.solve(mesh, ground, case, args.steps,
                        progress=lambda t, n, e: print(f'Compatible progress {t:.9f}, steps {n}, events {e}', flush=True),
                        identity=identity, system_factory=sparse.SparseCompatibleSystem, frame_callback=save_frame)
    assert bindings == {str(path): sha(path) for path in inputs+sources}, 'input or reference source changed'
    result.update(reference_steps=args.steps, input_sha256=bindings,
                  integration_implementation='cached step-start midpoint; unchanged events and yield-drift limits',
                  goal_complete=False, p6_complete=False, spatial_convergence_accepted=False)
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, allow_nan=False)
    print('Independent sparse compatible MCC construction completed', flush=True)


if __name__ == '__main__':
    main()
