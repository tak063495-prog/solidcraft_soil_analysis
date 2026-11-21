"""Bind frozen A-01 helpers without leaking generic import names into tests."""
import importlib.util
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
OLD=ROOT/'interim/2026-11-23-joint-combo-extension'
sys.path.insert(0,str(ROOT/'tools'))

def bind():
    saved={name:sys.modules.get(name) for name in ('fixtures','reference','fields','study')}
    loaded={}
    try:
        for name in saved:
            spec=importlib.util.spec_from_file_location('a01cause_legacy_'+name,OLD/(name+'.py'))
            module=importlib.util.module_from_spec(spec)
            sys.modules[name]=module
            spec.loader.exec_module(module)
            loaded[name]=module
    finally:
        for name,value in saved.items():
            if value is None:sys.modules.pop(name,None)
            else:sys.modules[name]=value
    return loaded
legacy=bind()
fixtures=legacy['fixtures'];fields=legacy['fields'];study=legacy['study'];original_reference=legacy['reference']
fem=fixtures.fem
