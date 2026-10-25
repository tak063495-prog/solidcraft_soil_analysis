"""Compare unchanged MCC physics after typed MC/MCC dispatch integration."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def verify(baseline,current,output):
    if output.exists():raise ValueError('comparison exists')
    a=json.loads(baseline.read_text());b=json.loads(current.read_text())
    assert a['complete'] and b['complete'] and a['ground']==b['ground'] and a['case']==b['case']
    assert a['last_committed']==b['last_committed'],'MCC final displacement/material/load history changed'
    assert a['frames']==b['frames'],'MCC stage equilibrium/monitor results changed'
    report={'schema':'soil-fem-mcc-dispatch-regression/1','exact_committed_state':True,'exact_frames':True,
        'baseline_solver_sha256':a['provenance']['solver_binary_sha256'],'current_solver_sha256':b['provenance']['solver_binary_sha256'],
        'files_sha256':{str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in (baseline,current,Path(__file__))}}
    output.write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('baseline',type=Path);p.add_argument('current',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();verify(a.baseline,a.current,a.output)
