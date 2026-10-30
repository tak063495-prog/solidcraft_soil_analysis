"""Recheck the pinned1200-element static witness with unchanged independent budgets."""
from pathlib import Path
import argparse,json
from independent_srm_slope_work import audit
a=argparse.ArgumentParser();a.add_argument('mesh',type=Path);a.add_argument('problem',type=Path);a.add_argument('result',type=Path);a.add_argument('output',type=Path)
args=a.parse_args();assert not args.output.exists()
proof=audit(args.mesh,args.problem,args.result,1.45);assert proof['passed']
with args.output.open('x',encoding='utf-8') as f:json.dump(proof,f,indent=2,allow_nan=False)
print('PASS pinned1200-element F1.45 full32400GP static exclusion; not exactFs or continuum safety')
