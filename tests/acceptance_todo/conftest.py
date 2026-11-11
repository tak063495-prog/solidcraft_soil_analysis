"""Local-only pytest bootstrap; deliberately no package __init__.py.

unittest discover (the unchanged CI command) does not descend into this folder.
No hooks alter existing tests, collection, markers or counts.
"""
import os
from pathlib import Path
import sys

for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(key, '1')
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tools'))

import pytest


@pytest.fixture(scope='session')
def solver():
    binary = Path(os.environ.get('ACCEPTANCE_TODO_SOLVER', ROOT/'target'/'release'/('soil-fem-3d.exe' if os.name=='nt' else 'soil-fem-3d')))
    assert binary.is_file(), f'Build first: cargo build --locked --release; missing {binary}'
    return binary.resolve()
