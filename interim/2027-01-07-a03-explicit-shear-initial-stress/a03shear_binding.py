"""Reuse frozen input/observation helpers, with restored search path."""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FOLLOW = ROOT/'interim/2026-12-30-a03-wet-soil-birth-cloud-followup'
_saved = sys.path[:]
try:
    sys.path.insert(0, str(FOLLOW))
    from a03follow_binding import study, fields, reference
    from a03follow_controls import controls, Sink, layout, canonical, compare, common, metric, norm, dot
    from a03follow_report import raw_observations, load, values, order
finally:
    sys.path[:] = _saved
