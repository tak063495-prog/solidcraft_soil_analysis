"""Isolate frozen A-04 imports; never install a hydraulic/JOINT adapter."""
import importlib.util
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
OLD = ROOT/'interim/2026-11-29-richards-wet-joint-a04'


def bind():
    names = ('soil_reference', 'study')
    saved = {name: sys.modules.get(name) for name in names}
    path = sys.path[:]
    loaded = {}
    try:
        sys.path.insert(0, str(ROOT/'tools'))
        for name in names:
            spec = importlib.util.spec_from_file_location('a04cause_legacy_'+name, OLD/(name+'.py'))
            module = importlib.util.module_from_spec(spec)
            sys.modules[name] = module
            spec.loader.exec_module(module)
            loaded[name] = module
    finally:
        sys.path[:] = path
        for name, previous in saved.items():
            if previous is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = previous
    return loaded
