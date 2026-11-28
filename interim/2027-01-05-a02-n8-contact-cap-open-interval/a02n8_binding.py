"""Bind frozen A-02 helpers under their existing unique namespaces."""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FOLLOW = ROOT/'interim/2026-12-29-a02-joint-birth-cloud-followup'
_saved = sys.path[:]
try:
    sys.path.insert(0, str(FOLLOW))
    from a02follow_binding import study, fields, reference, metric, norm, dot, common, compare
    from a02follow_policies import Predictor, serial, required
    from a02follow_report import switches, switchdiff
finally:
    sys.path[:] = _saved
