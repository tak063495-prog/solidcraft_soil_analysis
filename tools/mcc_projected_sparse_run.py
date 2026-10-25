"""Explicit sparse reference entry point; dense diagnostics remain the default."""
import argparse
from pathlib import Path
import mcc_projected_sparse_reference as sparse
import mcc_projected_volume_diagnostic as initial
import mcc_projected_cycle_diagnostic as cycle


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('phase', choices=('initial', 'cycle'))
    for name in ('mesh', 'ground', 'case', 'identity', 'output'):
        parser.add_argument(name, type=Path)
    parser.add_argument('--initial', type=Path)
    parser.add_argument('--steps', type=int, default=100)
    parser.add_argument('--fraction', type=float, default=1.)
    backend = parser.add_mutually_exclusive_group()
    backend.add_argument('--factor-reuse', action='store_true',
        help='experimental corrected-factor reference; full nonlinear CAD-path equivalence is not yet accepted')
    backend.add_argument('--matrix-free-factor-reuse', action='store_true',
        help='experimental original-operator reference with corrected cached LU and assembled fallback')
    parser.add_argument('--dependency-lock', type=Path, required=True)
    args = parser.parse_args()
    paths = [getattr(args, name).resolve() for name in ('mesh', 'ground', 'case', 'identity')]
    lock = args.dependency_lock.resolve()
    runtime_lock = lock.parent/'requirements.lock.txt'
    pins = dict(line.split('==') for path in (lock, runtime_lock)
                for line in path.read_text(encoding='utf-8').splitlines() if '==' in line)
    if pins.get('scipy') != sparse.scipy.__version__ or pins.get('numpy') != sparse.np.__version__:
        parser.error('reference NumPy/SciPy versions differ from the supplied dependency locks')
    extra = [Path(__file__).resolve(), Path(sparse.__file__).resolve(), lock, runtime_lock]
    factory = sparse.SparseProjectedSystem
    if args.factor_reuse:
        import mcc_projected_factor_reference as reused
        factory = reused.ReusedFactorProjectedSystem
        extra.append(Path(reused.__file__).resolve())
    if args.matrix_free_factor_reuse:
        import mcc_projected_factor_reference as reused
        import mcc_projected_matrixfree_reference as matrixfree
        factory = matrixfree.MatrixFreeFactorProjectedSystem
        extra.extend([Path(reused.__file__).resolve(), Path(matrixfree.__file__).resolve()])
    if args.phase == 'initial':
        initial.solve(*paths, args.output.resolve(), args.fraction, args.steps, factory, extra)
    else:
        if args.initial is None or args.fraction != 1.:
            parser.error('cycle requires --initial and --fraction 1')
        cycle.solve(*paths, args.initial.resolve(), args.output.resolve(), args.steps, factory, extra)


if __name__ == '__main__':
    main()
