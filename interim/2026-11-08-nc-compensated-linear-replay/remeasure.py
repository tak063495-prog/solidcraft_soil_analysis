"""Hash-check and remeasure saved evidence; executes no Rust/constitutive solver."""
import argparse
import hashlib
import json
from pathlib import Path
from audit_linear import audit
from measure_trace import measure
HERE=Path(__file__).resolve().parent


def remeasure(evidence,output):
    manifest=json.loads((HERE/'archive-sha256.json').read_text())
    for n,h in manifest['entries'].items():
        assert hashlib.sha256((evidence/n).read_bytes()).hexdigest()==h,('changed evidence',n)
    output.mkdir(parents=True,exist_ok=False)
    matches={}
    for request,response,name in [('matrices.json','linear-results.json','linear-audit.json'),
                                  ('candidate-matrices.json','candidate-linear-results.json','candidate-linear-audit.json')]:
        r=audit(json.loads((evidence/request).read_text()),json.loads((evidence/response).read_text()))
        matches[name]=r==json.loads((evidence/name).read_text())
        (output/name).write_text(json.dumps(r,indent=2)+'\n')
    measure(evidence/'evidence',evidence/'candidate-trial',output/'measurements')
    for name in ('trial-measurements.json','material-history.json'):
        matches[name]=json.loads((output/'measurements'/name).read_text())==json.loads((evidence/'measurements'/name).read_text())
    result=dict(passed=False,history_converged=False,solver_executed=False,material_map_executed=False,
                archive_entry_hashes_verified=len(manifest['entries']),exact_saved_json_matches=matches,
                windows_ci='unverified')
    (output/'remeasurement.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--evidence',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    print(json.dumps(remeasure(a.evidence.resolve(),a.output.resolve()),indent=2))
