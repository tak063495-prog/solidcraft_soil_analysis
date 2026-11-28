"""Record an incomplete sparse evidence dependency before invoking saved controls."""
import argparse
import json
from pathlib import Path
import sys
from pcfail_evidence import manifest, verify


def run(output,solver):
    output.mkdir(parents=True,exist_ok=False)
    folder,index=manifest('b01-cli');result=verify(folder,index)
    (output/'evidence-status.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    if result['status']!='verified':return 2 if result['status']=='incomplete' else 1
    sys.path.insert(0,str(folder.parent/'2026-12-25-b01-curved-tri6-cause-diagnosis'))
    from b01cause_cli_controls import run as controls
    controls(output/'controls',solver)
    value=json.loads((output/'controls/cli-controls-observations.json').read_text())
    return 0 if value['all_original_contracts_confirmed'] else 1

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--solver',type=Path,required=True);a=p.parse_args()
    raise SystemExit(run(a.output.resolve(),a.solver.resolve()))
