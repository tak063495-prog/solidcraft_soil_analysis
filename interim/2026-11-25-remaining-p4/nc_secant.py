"""Recover the exact incoming secant by replaying the P3 run with a recorder."""
import gzip
import inspect
import json
from pathlib import Path
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
P3=HERE.parent/'2026-11-25-remaining-p3'
sys.path.insert(0,str(P3))
import nc_branch

inputs,binary,out,source=map(Path,sys.argv[1:])
code=inspect.getsource(nc_branch.run)
needle='    checkpoint("last")'
assert code.count(needle)==1
code=code.replace(needle,needle+'\n    np.save(out/"secant.npy",t)')
namespace=dict(vars(nc_branch))
exec(compile(code,nc_branch.__file__,'exec'),namespace)
namespace['run'](inputs,binary.resolve(),out,0.02 if source.name=='0.02' else 0.04)
saved=json.load(gzip.open(P3/'nc'/source/'last.json.gz'))
actual=json.load(gzip.open(out/'last.json.gz'))
assert saved==actual
