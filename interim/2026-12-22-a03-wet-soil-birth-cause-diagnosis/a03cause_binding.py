"""Read frozen A-03 helpers without leaking generic import names."""
import importlib.util
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
OLD = ROOT/'interim/2026-11-26-wet-soil-birth-darcy-a03'
if str(ROOT/'tools') not in sys.path:sys.path.insert(0, str(ROOT/'tools'))


def bind():
    names = ('fields', 'reference', 'study')
    saved = {name:sys.modules.get(name) for name in names}
    loaded = {}
    try:
        for name in names:
            spec = importlib.util.spec_from_file_location('a03cause_legacy_'+name, OLD/(name+'.py'))
            module = importlib.util.module_from_spec(spec)
            sys.modules[name] = module
            spec.loader.exec_module(module)
            loaded[name] = module
    finally:
        for name, previous in saved.items():
            if previous is None:sys.modules.pop(name, None)
            else:sys.modules[name] = previous
    return loaded


legacy = bind()
fields, reference, study = (legacy[name] for name in ('fields', 'reference', 'study'))
