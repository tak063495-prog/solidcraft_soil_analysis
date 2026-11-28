"""Prepare the OLD rainfall RT0 space/time columns into a fresh directory.

Input definition transcribed from the unchanged rain-original/run.py. No solve.
"""
import argparse
from collections import Counter
import copy
import hashlib
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
from richards_dam_mesh import make
SOURCE=ROOT/'interim/2026-10-29-native-kernel-and-rain-gates/rain-original'

def prepare(output):
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    base=json.loads((SOURCE/'reference-case.json').read_text(encoding='utf-8-sig'))
    for n,dt in ((4,100),(8,100),(16,100),(8,50),(8,25)):
        folder=output/f'rt0-{n}-dt{dt}';folder.mkdir()
        mesh=make(n)
        faces=Counter(tuple(sorted(t[i] for i in face)) for t in mesh['tetrahedra']
                      for face in ((1,2,3),(0,3,2),(0,1,3),(0,2,1)))
        mesh['face_sets']=dict(bottom=[],rain=[],dry=[])
        for face,count in faces.items():
            if count!=1:continue
            xyz=[mesh['nodes'][i] for i in face]
            if all(q[2]==0. for q in xyz):mesh['face_sets']['bottom'].append(list(face))
            elif all(q[2]==1. for q in xyz):
                patch=all(.25<=q[0]<=.75 and .25<=q[1]<=.75 for q in xyz)
                mesh['face_sets']['rain' if patch else 'dry'].append(list(face))
        case=copy.deepcopy(base)
        case['settings']=dict(initial_time_step_s=float(dt),maximum_time_step_s=float(dt),
                              absolute_flow_tolerance_m3_s=1e-17,store_integration_points=False)
        for period in case['analysis']['periods']:
            top=period['boundaries'].pop()
            assert top['face_set']=='top'
            top['face_set']='rain'
            period['boundaries'].extend([top,dict(face_set='dry',condition=dict(kind='atmospheric_no_ponding',rainfall_vertical_m_s=0.))])
        for name,value in (('mesh.json',mesh),('case.json',case)):
            (folder/name).write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8')
    binding=dict(source_paths=['rain-original/reference-case.json','rain-original/run.py'],
                 source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (SOURCE/'reference-case.json',SOURCE/'run.py')},
                 numerical_solve_performed=False,p6_03_complete=False)
    (output/'preparation.json').write_text(json.dumps(binding,indent=2)+'\n',encoding='utf-8')
    print('Prepared old RT0 4/8/16 space and 8-fixed 100/50/25s time cases; no solve')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);prepare(p.parse_args().output)
