"""Compare independently marched branches on the shared accumulated arclength."""
import json
from pathlib import Path
import sys
import numpy as np

folder=Path(sys.argv[1])
routes=[]
for name in ('fine-extended-retry','extended-004'):
    route=folder/name
    meta=json.loads((route/'summary.json').read_bytes())
    rows=[json.loads(line) for line in (route/'trials.jsonl').open()]
    accepted=[r for r in rows if r['accepted']]
    assert len(rows)==meta['attempts'] and len(accepted)==meta['accepted']
    events=dict(nonnegative_increments=sum(r['delta_lambda']>=0 for r in accepted),
                nonnegative_tangents=sum(r['lambda_tangent']>=0 for r in accepted),
                reached_one=any(r['load']>=1 for r in accepted))
    error_counts={}
    for row in rows:
        if not row['accepted']:
            key=row['reason'].strip()
            error_counts[key]=error_counts.get(key,0)+1
    routes.append(dict(name=name,events=events,arc=[r['arc_m'] for r in accepted],
                       load=[r['load'] for r in accepted],rejections=error_counts))
a,b=routes
common=np.array(a['arc'])
mask=(common>=b['arc'][0])&(common<=b['arc'][-1])
difference=np.interp(common[mask],b['arc'],b['load'])-np.array(a['load'])[mask]
print(json.dumps(dict(overlapping_points=int(mask.sum()),max_abs_lambda_interpolation_difference=float(np.max(abs(difference))),
                      final_common_arc=float(common[mask][-1]),events={r['name']:r['events'] for r in routes},rejections={r['name']:r['rejections'] for r in routes}),indent=2))
