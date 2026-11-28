"""Observe localization and the attachment rotation in the unchanged load source."""
import argparse
import json
import sys
from pathlib import Path
D=Path(__file__).resolve().parent
if D.name=='work':sys.path.insert(0,str(D.parent/'interim/2026-12-09-b02-tet768-reassembly-residual'))
from common import read,write
from plate6_offset_reference import System
import numpy as np

def observe(formal,analysis,output):
    results=[]
    for case in read(formal/'fresh-solves.json')['cases']:
        folder=formal/case['name'];m,c=read(folder/'mesh.json'),read(folder/'stages.json');model=System(m,c['beam_system'])
        attachment=model.data[0][0][0];point=model.positions[attachment];tip=model.node_ids['tip'];arm=model.positions[tip]-point
        fs=read(analysis/(case['name']+'-fields.json'))
        for kind in ('native','independent'):
            endpoint=next(r for r in fs[kind] if r['stage_index']==2 and r['progress']==1.)
            q=np.asarray(endpoint['shear_direct_global_kn_m']).reshape(-1,3);queries=np.asarray(endpoint['queries_m'])
            spatial=[dict(probe=i,position_m=p.tolist(),in_plane_distance_from_beam_m=float(np.linalg.norm((p-point)[:2])),direct_shear_global_kn_m=v.tolist(),norm_kn_m=float(np.linalg.norm(v))) for i,(p,v) in enumerate(zip(queries,q,strict=True))]
            row=dict(case=case['name'],level=case['level'],steps=case['steps'],kind=kind,common_probes=spatial,loaded_work_components_kn_m=endpoint['work_components_kn_m'])
            if 'load_source_probes' in endpoint:
                probes=endpoint['load_source_probes'];a=np.asarray(probes['nodal_attachment_response']);t=np.asarray(probes['nodal_tip_response']);force=np.zeros(3)
                for load in c['beam_system']['nodal_loads']:
                    if load['node_id']=='tip':force+=load.get('boundary_force_kn',[0.,0.,0.])
                rotational=np.cross(a[3:],arm);flex=t[:3]-a[:3]-rotational
                row['nodal_source_tip_work_kn_m']=dict(attachment_translation=float(force@a[:3]),attachment_rotation=float(force@rotational),beam_flexure=float(force@flex),sum=float(force@t[:3]))
                row['nodal_source_shear_global_kn_m']=endpoint['nodal_source_shear_global_kn_m'];row['other_source_shear_global_kn_m']=endpoint['other_source_shear_global_kn_m']
            source=folder/('native/accepted.jsonl' if kind=='native' else 'independent.jsonl')
            with source.open() as f:
                for line in f:
                    value=json.loads(line)
                    if value['stage_index']!=2 or value['progress']!=1.:continue
                    report=value['beam_results'];plates=report['plates'];plates={p['id']:p for p in plates} if isinstance(plates,list) else plates
                    raw=[];energies={k:0. for k in ('membrane_energy_kn_m','bending_energy_kn_m','transverse_shear_energy_kn_m','drilling_penalty_energy_kn_m')}
                    for pi,(nodes,_,kernel) in enumerate(model.plate_data):
                        p=plates[model.plates[pi]['id']]
                        for key in energies:energies[key]+=p[key]
                        for qi,s in enumerate(p['samples']):
                            pos=np.asarray(s['barycentric'])@model.positions[nodes][:3]
                            raw.append(dict(plate=model.plates[pi]['id'],sample=qi,attachment_position_m=pos.tolist(),in_plane_distance_from_beam_m=float(np.linalg.norm((pos-point)[:2])),norm_kn_m=float(np.linalg.norm(s['transverse_shear_force_kn_m']))))
                    row['plate_energy_kn_m']=energies;row['largest_raw_shear_samples']=sorted(raw,key=lambda v:v['norm_kn_m'],reverse=True)[:12]
                    break
            results.append(row)
    write(output,dict(passed=False,b02_complete=False,observations=results,note='Finite three-level localization evidence; no proof of continuum divergence, no new tolerance or recovery replacement.'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('formal',type=Path);p.add_argument('analysis',type=Path);p.add_argument('output',type=Path);a=p.parse_args();observe(a.formal.resolve(),a.analysis.resolve(),a.output.resolve())
