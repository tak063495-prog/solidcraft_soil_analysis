"""Restart from a committed P3 point, with the original bordered corrector."""
import gzip
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np

HERE = Path(__file__).resolve().parent
P3 = HERE.parent/'2026-11-25-remaining-p3'
sys.path.insert(0,str(P3))
from nc_branch import Section, correct


def read(path):
    with gzip.open(path,'rt') as stream:
        return json.load(stream)


def run(inputs, binary, out, source):
    start = time.monotonic()
    out.mkdir(parents=True)
    m = Section(inputs,binary)
    src = P3/'nc'/source
    path = src/'last.json.gz'
    parent = read(path)
    old = parent['active_material_states']
    lam, arc = parent['load'], parent['arc_m']
    u = np.ravel(parent['displacement_m'])
    # Evaluate the local Jacobian with zero displacement increment. Orient its
    # null direction toward decreasing lambda, as observed on the incoming path.
    _, k, _, _ = m.evaluate(old,np.zeros(len(m.free)),np.ones(len(old),dtype=int))
    t = np.r_[np.linalg.solve(k,m.h),-1.]
    t /= np.linalg.norm(t)
    ds = .02
    rows = []
    inherited = json.loads((src/'summary.json').read_bytes())['attempts']+1
    stop = 'original attempt cap'
    with (out/'trials.jsonl').open('x',buffering=1) as log:
        for attempt in range(1,2001-inherited):
            result, its, reason = correct(m,old,lam,t,ds,min(1e-12,m.settings['relative_force_tolerance']))
            row = dict(attempt=attempt,parent_lambda=lam,ds_m=ds,iterations=its,accepted=result is not None)
            if result is None:
                row['reason'] = reason
                ds /= 2
                if ds < m.settings['minimum_increment']:
                    stop = 'failed corrector at minimum step'
            else:
                dz,old,diag = result
                u[m.free] += dz[:-1]
                lam += dz[-1]
                arc += ds
                t = dz/np.linalg.norm(dz)
                row.update(diag,load=lam,delta_lambda=dz[-1],arc_m=arc)
                if dz[-1]>=0 or diag['lambda_tangent']>=0 or lam>=1:
                    stop = 'load event bracket'
                ds = min(.02,ds*(1.3 if len(its)<=8 else 1.))
            log.write(json.dumps(row,allow_nan=False)+'\n')
            rows.append(row)
            if stop!='original attempt cap':
                break
    with gzip.open(out/'last.json.gz','wt') as f:
        json.dump(dict(load=lam,arc_m=arc,displacement_m=u.reshape(-1,3).tolist(),active_material_states=old),f,allow_nan=False)
    accepted = [r for r in rows if r['accepted']]
    summary = dict(source=str(path.relative_to(HERE.parents[1])),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                   start_lambda=parent['load'],last_lambda=lam,arc_m=arc,stop=stop,inherited_attempts=inherited,
                   attempts=len(rows),accepted=len(accepted),elapsed_s=time.monotonic()-start)
    for key in ('sigma_next','bordered_min','lambda_tangent','delta_lambda','relative_force'):
        if accepted:
            values=[r[key] for r in accepted]
            summary[key]=[min(values),max(values)]
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':
    run(Path(sys.argv[1]),Path(sys.argv[2]).resolve(),Path(sys.argv[3]),sys.argv[4])
