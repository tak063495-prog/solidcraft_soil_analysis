"""Prepare documented MCC column series only; do not run a solver."""
import argparse
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
from prepare_mixed_mcc_reference import column, prepare

def generate(output):
    output=Path(output);prepare(output)
    for n in (8,16,32):
        p=output/f'n{n}';p.mkdir()
        (p/'mesh.json').write_text(json.dumps(column(n)),encoding='utf-8')
        (p/'stages.json').write_bytes((output/'n4/stages.json').read_bytes())
    for label,n,load,material in (
        ('space8',8,1/32,1e-5),('space16',16,1/32,1e-5),('space32',32,1/32,1e-5),
        ('load16',16,1/16,1e-5),('load32',16,1/32,1e-5),('load64',16,1/64,1e-5),
        ('mat1e4',16,1/32,1e-4),('mat1e5',16,1/32,1e-5),('mat1e6',16,1/32,1e-6)):
        p=output/label;p.mkdir()
        (p/'mesh.json').write_bytes((output/f'n{n}/mesh.json').read_bytes())
        case=json.loads((output/f'n{n}/stages.json').read_text())
        case['settings'].update(initial_increment=load,material_maximum_substep=material)
        (p/'stages.json').write_text(json.dumps(case,indent=2)+'\n',encoding='utf-8')
    print('Prepared 9 cases, no solve; same documented MCC physics; separate space/load/material axes')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);generate(p.parse_args().output)
