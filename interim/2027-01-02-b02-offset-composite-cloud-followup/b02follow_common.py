"""Explicit small diagnostic; frozen helpers are reused, never imported by production."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import json
import hashlib
from contextlib import contextmanager
from types import SimpleNamespace
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OLD = ROOT/'interim/2026-12-01-tri6-offset-composite-b02'
CAUSE = ROOT/'interim/2026-12-24-b02-offset-composite-cause-diagnosis'
sys.path.insert(0, str(CAUSE))
from b02cause_binding import legacy
from b02cause_driver import command, sha
from b02cause_metrics import dot, norm, metric, common
import plate6_reference


def write(path, obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')


@contextmanager
def quadrature(order):
    if order not in (4, 5, 6):
        raise ValueError('bounded plate rule only')
    original = plate6_reference.np
    # Replace only the plate module namespace. Soil/beam numpy is untouched.
    legendre = SimpleNamespace(**{**vars(np.polynomial.legendre), 'leggauss': lambda n: np.polynomial.legendre.leggauss(order if n == 5 else n)})
    polynomial = SimpleNamespace(**{**vars(np.polynomial), 'legendre': legendre})
    plate6_reference.np = SimpleNamespace(**{**vars(np), 'polynomial': polynomial})
    try:
        yield
    finally:
        plate6_reference.np = original


def fingerprint(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def source_binding():
    import subprocess
    names = subprocess.check_output(['git','ls-files'], cwd=ROOT, text=True).splitlines()
    names += [str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.rs','.md')]
    return {n: sha(ROOT/n) for n in sorted(set(names)) if Path(n).suffix in ('.py','.rs','.toml','.ps1') or Path(n).name == 'Cargo.lock'}
