"""Isolate frozen B-02 helpers; no production entry point imports this module."""
import importlib.util
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
OLD = ROOT/'interim/2026-12-01-tri6-offset-composite-b02'


def bind():
    names = ('cases', 'compare', 'reference', 'study', 'work_audit')
    saved = {n: sys.modules.get(n) for n in names}
    path = sys.path[:]
    loaded = {}
    try:
        sys.path.insert(0, str(ROOT/'tools'))
        for name in names:
            spec = importlib.util.spec_from_file_location('b02cause_legacy_'+name, OLD/(name+'.py'))
            module = importlib.util.module_from_spec(spec)
            sys.modules[name] = module
            spec.loader.exec_module(module)
            loaded[name] = module
    finally:
        sys.path[:] = path
        for name, value in saved.items():
            if value is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = value
    return loaded


legacy = bind()
