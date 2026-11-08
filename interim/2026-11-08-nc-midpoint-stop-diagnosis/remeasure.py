"""Verify archived bytes and remeasure stored evidence without running any solver."""
import argparse
import hashlib
import json
from pathlib import Path
from measure import measure
HERE=Path(__file__).resolve().parent

def remeasure(evidence,out):
    manifest=json.loads((HERE/'archive-sha256.json').read_text())
    for n,h in manifest['entries'].items():
        assert hashlib.sha256((evidence/n).read_bytes()).hexdigest()==h,('evidence changed',n)
    out.mkdir(parents=True,exist_ok=False)
    r=measure(evidence/'prior',evidence/'diagnostic',None,out/'recomputed',evidence/'measurements')
    matches={n:json.loads((evidence/'measurements'/n).read_text())==json.loads((out/'recomputed'/n).read_text())
             for n in ('measurements.json','iterate-measurements.json','material-history.json','contact-measurements.json')}
    report=dict(passed=False,history_converged=False,solver_executed=False,material_map_executed=False,
                archive_entry_hashes_verified=len(manifest['entries']),exact_saved_json_matches=matches,
                physical_acceptance_granted=False,windows_ci='unverified')
    (out/'remeasurement.json').write_text(json.dumps(report,indent=2))
    return report
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--evidence',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    print(json.dumps(remeasure(a.evidence.resolve(),a.output.resolve()),indent=2))
