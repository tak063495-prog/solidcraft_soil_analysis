"""Explicit observation helpers. No production import or dispatch."""
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time
for k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[k] = '1'
import numpy as np
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE = '71655b7b09e1895e364ba476b9089075ba7130c0'

def read(p):
    return json.loads(Path(p).read_text(encoding='utf-8-sig'))

def convert(v):
    if isinstance(v, np.ndarray):
        return v.astype(float).tolist() if v.dtype.kind == 'f' else v.tolist()
    if isinstance(v, np.generic):
        return float(v) if v.dtype.kind == 'f' else v.item()
    raise TypeError(type(v).__name__)

def write(p, v):
    with Path(p).open('x', encoding='utf-8', newline='\n') as f:
        json.dump(v, f, default=convert, ensure_ascii=False, indent=2, allow_nan=False)
        f.write('\n')

def sha(p):
    with Path(p).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def load(name, p):
    spec = importlib.util.spec_from_file_location(name, p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m

def invoke(args, folder, name, switches=None, expected=0):
    env = {k:v for k,v in os.environ.items() if not k.startswith('C04FOLLOW_')}
    env.update(switches or {})
    started = time.monotonic()
    with (folder/(name+'.log')).open('x', encoding='utf-8') as f:
        p = subprocess.Popen(list(map(str,args)), cwd=ROOT, env=env, stdout=f, stderr=subprocess.STDOUT)
        timed = False
        try:
            code = p.wait(timeout=int(os.environ.get('C04FOLLOW_PROCESS_SECONDS','3600')))
        except subprocess.TimeoutExpired:
            timed = True
            p.kill()
            code = p.wait()
    result = dict(command=list(map(str,args)), pid=p.pid, exit_code=code, timed_out=timed,
                  elapsed_seconds=time.monotonic()-started, switches=switches or {})
    write(folder/(name+'-process.json'), result)
    if code != expected or timed:
        raise RuntimeError(('numerically incomplete; retain prefix and logs', result))
    return result

def rows(p):
    with Path(p).open(encoding='utf-8') as f:
        for line in f:
            yield json.loads(line)

def norm(v):
    return math.sqrt(math.fsum(float(x)*float(x) for x in np.asarray(v).ravel(order='C')))

def compare(a,b):
    x=np.asarray(a,dtype=float); y=np.asarray(b,dtype=float)
    if x.shape != y.shape or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError(('field shape/nonfinite',x.shape,y.shape))
    d=x-y; dn=norm(d); rn=norm(y)
    return dict(absolute_l2=dn, reference_l2=rn, relative_l2=dn/max(rn,1e-12),
                maximum_absolute=float(np.max(np.abs(d),initial=0)), components=x.size)

def sink(p):
    f=Path(p).open('x',encoding='utf-8',newline='\n')
    def emit(row,fields):
        f.write(json.dumps(dict(row=row,fields=fields),default=convert,allow_nan=False,separators=(',',':'))+'\n')
        f.flush()
    return f,emit
