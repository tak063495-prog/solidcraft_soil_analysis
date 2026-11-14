"""Repeat the SAME 21 saved jump paths at x64/x128, then compare vectors.

No global solve and no extrapolation. Last finite-reference differences are
truncation indicators, never rigorous upper bounds or new acceptance gates.
"""
import argparse
import json
from pathlib import Path
import shutil
import sys
import tempfile
import zipfile
from math import sqrt
import numpy as np

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];PRIOR=ROOT/'interim/2026-12-15-c01-local-path-decomposition'
sys.path.insert(0,str(PRIOR))
import c01_decomposition as d
import c01_refinement_audit as audit
from c01_local_study import write


def extract_old(out):
    import c01_local_evidence
    c01_local_evidence.verify()
    idx=d.read(PRIOR/'archive-index.json');out.mkdir(parents=True,exist_ok=False)
    with tempfile.TemporaryDirectory(dir=out) as tmp:
        arc=Path(tmp)/'old.zip'
        with arc.open('xb') as f:
            for part in idx['parts']:
                with (PRIOR/part['file']).open('rb') as src:shutil.copyfileobj(src,f)
        with zipfile.ZipFile(arc) as z:
            for name in ('decomposition/decomposition.json','decomposition/point-increments.jsonl'):
                target=out/Path(name).name
                with z.open(name) as f,target.open('xb') as dest:shutil.copyfileobj(f,dest)
                entry=next(e for e in idx['entries'] if e['path']==name)
                if d.sha(target)!=entry['sha256']:raise ValueError('old selected trace SHA')


def compact(trace,tags):
    result={}
    with trace.open(encoding='utf-8') as f:
        for line in f:
            r=json.loads(line);key=(r['control'],r['kind'],r['steps']);coord=(r['stage'],r['progress'])
            row={field:{tag:np.asarray(r['fields'][field][tag]) for tag in tags} for field in d.FIELDS}
            rows=result.setdefault(key,{})
            if coord in rows:raise ValueError('duplicate coordinate')
            rows[coord]=row
    return result


def components(values,coordinates,intervals,field,tag):
    aa,bb=values;stored_a=np.stack([aa[c][field]['stored'] for c in coordinates]);stored_b=np.stack([bb[c][field]['stored'] for c in coordinates]);delta=stored_b-stored_a
    path=np.stack([bb[c][field][tag]-aa[c][field][tag] for c in coordinates]);sw=np.zeros_like(path[0]);before=np.zeros_like(sw);rows=[]
    for i,interval in enumerate(intervals):
        mask=np.zeros(len(sw),dtype=bool);mask[interval['switch_timing_mismatch_gp']]=True
        if sw.ndim==2:mask=mask[:,None]
        sw+=np.where(mask,path[i]-before,0.);rows.append(sw.copy());before=path[i]
    switch=np.stack(rows)
    return delta,dict(zip(d.COMPONENTS,(delta-path,path-switch,switch)))


def compare(out):
    old=out/'old32';extract_old(old)
    report=d.read(out/'x128/decomposition.json')
    traces={32:compact(old/'point-increments.jsonl',('stored','coarser_replay','refined_replay')),64:compact(out/'x64/point-increments.jsonl',('stored','refined_replay')),128:compact(out/'x128/point-increments.jsonl',('stored','refined_replay'))}
    result=dict(schema='c01-regen-finite-reference/1',c01_passed=False,strict_upper_bound=False,extrapolation_used=False,reference_multipliers=[16,32,64,128],pairs={},other_iv=dict(vector_identifiable=False,reason='nonlinear interaction, unresolved true between-endpoint path, finite contact-event history; no unique separate vector',finite_reference_indicator='norm of x128 minus x64 component, without Richardson correction'))
    for label,pair in report['pairs'].items():
        control,kind,levels=label.split('/');a,b=map(int,levels.split('__'));coords=[(s,i/a) for s in d.STAGES for i in range(1,a+1)];fields={}
        for field in d.FIELDS:
            versions={};delta=None
            for n,src,tag in ((16,32,'coarser_replay'),(32,32,'refined_replay'),(64,64,'refined_replay'),(128,128,'refined_replay')):
                delta,versions[n]=components((traces[src][(control,kind,a)],traces[src][(control,kind,b)]),coords,pair['intervals'],field,tag)
            total=d.dot(delta,delta);projections={str(n):{c:d.dot(v,delta)/total if total else None for c,v in cs.items()} for n,cs in versions.items()}
            changes={}
            for lo,hi in ((16,32),(32,64),(64,128),(16,128),(32,128)):
                changes[f'{lo}__{hi}']={c:dict(squared_norm=d.dot(versions[hi][c]-versions[lo][c],versions[hi][c]-versions[lo][c]),norm_over_original_difference=sqrt(d.dot(versions[hi][c]-versions[lo][c],versions[hi][c]-versions[lo][c])/total) if total else None) for c in d.COMPONENTS}
            ratios={}
            for c in d.COMPONENTS:
                x=changes['32__64'][c]['squared_norm'];y=changes['64__128'][c]['squared_norm'];ratios[c]=sqrt(y/x) if x else None
            fields[field]=dict(original=pair['fields'][field]['original'],signed_projections=projections,component_vector_changes=changes,last_change_ratio=ratios,finite_truncation_indicator=changes['64__128'],maximum_identity_residual=max(float(np.max(abs(delta-sum(cs.values())))) for cs in versions.values()))
        result['pairs'][label]=fields
    write(out/'finite-summary.json',result)
    return result


def main():
    import subprocess
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--inputs',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--prepare',action='store_true');p.add_argument('--compare-only',action='store_true');a=p.parse_args()
    if not a.compare_only:
        a.output.mkdir(parents=True,exist_ok=False)
        if a.prepare:d.prepare(a.inputs)
        for n in (64,128):
            subprocess.run([sys.executable,str(PRIOR/'c01_decomposition.py'),'--inputs',str(a.inputs),'--output',str(a.output/f'x{n}'),'--refinement',str(n)],check=True)
            write(a.output/f'x{n}/finite-reference.json',audit.build(d.read(a.output/f'x{n}/decomposition.json'),a.output/f'x{n}/point-increments.jsonl'))
    compare(a.output)

if __name__=='__main__':main()
