"""Run the independent dry projected-MCC reference from original input files.

The matrix-free backend retains the original three-field equations and their
true residual gates. A completed original-input initialization can be supplied
explicitly for a separate construction run; its bindings remain mandatory.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

import mcc_projected_matrixfree_reference as matrixfree
import mcc_projected_factor_reference as factor
import mcc_projected_sparse_reference as sparse
import mcc_projected_volume_diagnostic as initialization
import mcc_projected_cycle_diagnostic as construction
import mcc_rate_rk4
import numpy
import scipy


class Rk4MatrixFreeProjectedSystem(matrixfree.MatrixFreeFactorProjectedSystem):
    advance = mcc_rate_rk4.advance_rk4


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(mesh, ground, case, identity, output, *, steps=100, integration='midpoint',
        phase='both', initial=None):
    if type(steps) is not int or steps < 2:
        raise ValueError('at least two independent integration steps required')
    if integration not in ('midpoint', 'rk4') or phase not in ('initial', 'cycle', 'both'):
        raise ValueError('unknown reference integration or phase')
    if (phase == 'cycle') != (initial is not None):
        raise ValueError('only a separate cycle phase requires --initial')
    mesh, ground, case, identity, output = [Path(p).resolve() for p in (mesh, ground, case, identity, output)]
    if read(case).get('schema') != 'soil-fem-projected-mcc-stages/1':
        raise ValueError('dry projected-MCC construction inputs required')
    root = Path(__file__).resolve().parents[1]
    locks = [root/name for name in ('requirements.lock.txt', 'requirements-reference.lock.txt')]
    pins = dict(line.split('==') for path in locks for line in path.read_text(encoding='utf-8').splitlines() if '==' in line)
    if numpy.__version__ != pins['numpy'] or scipy.__version__ != pins['scipy']:
        raise ValueError('install the pinned NumPy and SciPy reference dependencies')
    sources = [Path(__file__).resolve(), Path(matrixfree.__file__), Path(factor.__file__), Path(sparse.__file__), *locks]
    if integration == 'rk4':
        sources.append(Path(mcc_rate_rk4.__file__))
    factory = Rk4MatrixFreeProjectedSystem if integration == 'rk4' else matrixfree.MatrixFreeFactorProjectedSystem
    inputs = [mesh, ground, case, identity]
    original = Path(initial).resolve() if initial is not None else None
    watched = inputs + sources + ([] if original is None else [original])
    bindings = {str(p): sha(p) for p in watched}
    output.mkdir(parents=True, exist_ok=False)
    # The lower-level reference binds inventory roles by canonical filenames.
    # Preserve the user's files and freeze byte-identical, named inputs here.
    (output/'inputs').mkdir()
    frozen = []
    for source, name in zip(inputs, ('mesh.json', 'ground.json', 'case.json', 'mesh-identity.json'), strict=True):
        target = output/'inputs'/name
        shutil.copyfile(source, target)
        if sha(target) != bindings[str(source)]:
            raise ValueError('reference input changed while being copied')
        frozen.append(target)
    bindings.update({str(p): sha(p) for p in frozen})
    watched += frozen
    if phase != 'cycle':
        original = output/'initial.json'
        initialization.solve(*frozen, original, 1., steps, factory, sources)
    if phase != 'initial':
        construction.solve(*frozen, original, output/'cycle.json', steps, factory, sources)
    if bindings != {str(p): sha(p) for p in watched}:
        raise ValueError('reference input or source changed during execution')
    result = dict(schema='soil-fem-projected-mcc-reference-run/1', complete=True,
                  phase=phase, integration_method=integration, steps=steps,
                  reference_backend=factory.__module__+'.'+factory.__name__,
                  versions=dict(numpy=numpy.__version__, scipy=scipy.__version__),
                  input_and_driver_sha256=bindings, initial_reference_sha256=sha(original),
                  spatial_convergence_accepted=False, p6_complete=False)
    with (output/'run.json').open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('mesh', 'ground', 'case', 'identity', 'output'):
        parser.add_argument(name, type=Path)
    parser.add_argument('--steps', type=int, default=100)
    parser.add_argument('--integration', choices=('midpoint', 'rk4'), default='midpoint')
    parser.add_argument('--phase', choices=('initial', 'cycle', 'both'), default='both')
    parser.add_argument('--initial', type=Path)
    options = parser.parse_args()
    run(**vars(options))
