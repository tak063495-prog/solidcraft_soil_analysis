"""Isolate frozen A-02 helpers; no generic module binding escapes this import."""
import importlib.util
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
OLD = ROOT/'interim/2026-11-25-joint-birth-projection-a02'
sys.path.insert(0, str(ROOT/'tools'))


def bind():
    names = ('reference', 'fields', 'study', 'projection', 'geometry_audit')
    saved = {name: sys.modules.get(name) for name in names}
    loaded = {}
    try:
        for name in names:
            spec = importlib.util.spec_from_file_location('a02cause_legacy_'+name, OLD/(name+'.py'))
            module = importlib.util.module_from_spec(spec)
            sys.modules[name] = module
            spec.loader.exec_module(module)
            loaded[name] = module
    finally:
        for name, previous in saved.items():
            if previous is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = previous
    return loaded


legacy = bind()
reference = legacy['reference']
fields = legacy['fields']
study = legacy['study']
geometry = legacy['geometry_audit']
