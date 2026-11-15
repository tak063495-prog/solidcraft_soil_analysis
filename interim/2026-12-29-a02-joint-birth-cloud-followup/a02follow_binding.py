"""Reuse frozen A-02 helpers without leaking generic module aliases."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
CAUSE=ROOT/'interim/2026-12-21-a02-joint-birth-cause-diagnosis'
saved=sys.path[:]
try:
    sys.path.insert(0,str(CAUSE))
    from a02cause_binding import reference,fields,study,geometry,OLD
    from a02cause_metrics import metric,norm,dot,compare,common
finally:sys.path[:]=saved
