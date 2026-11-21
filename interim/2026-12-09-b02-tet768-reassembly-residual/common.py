"""Diagnostic-only imports. Existing B-02 inputs, gates and reference stay intact."""
from pathlib import Path
import os
import sys
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE = ROOT / 'interim/2026-12-01-tri6-offset-composite-b02'
sys.path.insert(0, str(BASE))
sys.path.insert(0, str(ROOT/'tools'))
from cases import read, write, inputs
from study import command, sha
