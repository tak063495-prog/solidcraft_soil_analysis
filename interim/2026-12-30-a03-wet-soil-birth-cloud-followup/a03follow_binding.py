"""Frozen helpers with distinct names and restored generic import bindings."""
import importlib.util
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OLD = ROOT/'interim/2026-11-26-wet-soil-birth-darcy-a03'
sys.path.insert(0, str(ROOT/'tools'))


def bind():
    names = ('fields', 'reference', 'study')
    saved = {n: sys.modules.get(n) for n in names}
    modules = {}
    try:
        for n in names:
            spec = importlib.util.spec_from_file_location('a03follow_legacy_'+n, OLD/(n+'.py'))
            module = importlib.util.module_from_spec(spec)
            sys.modules[n] = module
            spec.loader.exec_module(module)
            modules[n] = module
    finally:
        for n, previous in saved.items():
            if previous is None:
                sys.modules.pop(n, None)
            else:
                sys.modules[n] = previous
    return modules


legacy = bind()
fields, reference, study = (legacy[n] for n in ('fields', 'reference', 'study'))
