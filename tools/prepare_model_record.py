"""Record geometry-only preparation; never mark it an equilibrated ground."""
import argparse
import hashlib
import json
from pathlib import Path


def record(folder,accepted=False):
    if accepted:
        result=json.loads((folder/'initial-equilibrium.json').read_text())
        assert result['complete'] and not result['case']['stages'] and result['geostatic']['complete']
        verification=json.loads((folder/'accepted-stage-run.json').read_text());assert verification['solver_binary_sha256']==result['provenance']['solver_binary_sha256']
    files=[p for p in folder.iterdir() if p.is_file() and p.name not in ('prepared-model.json','accepted-run.json')]
    result={'schema':'soil-fem-accepted-ground-run/1' if accepted else 'soil-fem-prepared-model/1',
        'equilibrium_verified':accepted,'files_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
    with (folder/('accepted-run.json' if accepted else 'prepared-model.json')).open('x',encoding='utf-8') as f:json.dump(result,f,indent=2)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('folder',type=Path);p.add_argument('--accepted',action='store_true');a=p.parse_args();record(a.folder,a.accepted)
