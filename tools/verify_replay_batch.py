"""Audit unit replay, adjusted inputs, per-trial exports and independent solves."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from verify_construction_linear import verify


def audit(folder,output):
    if output.exists():raise ValueError('verification output exists')
    native=json.loads((folder/'replay-results.json').read_text());package=json.loads((folder/'trials/verification.json').read_text())
    baseline=json.loads((folder/'source-stage-results.json').read_text())
    assert package['source_report_sha256']==hashlib.sha256((folder/'replay-results.json').read_bytes()).hexdigest()
    assert len(native['trials'])==len(package['trials'])==len(native['specification']['factors'])
    assert native['source_ground']==baseline['ground'] and native['case']==baseline['case']
    assert native['specification']==json.loads((folder/'specification.json').read_text())
    assert native['stability_classification']==package['stability_classification']=='not_evaluated'
    assert native['replay_origin']=='original_inputs_before_initial_equilibration'
    digest=hashlib.sha256((folder/'solver/soil-fem-3d.exe').read_bytes()).hexdigest();assert native['source_identity']['solver_binary_sha256']==digest==package['solver_sha256']
    unit=next(t for t in native['trials'] if t['factor']==1.)
    assert unit['run']['last_committed']==baseline['last_committed'] and unit['run']['frames']==baseline['frames']
    original={m['id']:m for m in native['source_ground']['materials']};selected=native['specification']['material_ids'];records=[]
    files=[folder/'replay-results.json',folder/'source-stage-results.json',folder/'trials/verification.json',folder/'specification.json']
    for i,trial in enumerate(native['trials']):
        assert trial['outcome']=='completed' and trial['run']['complete'] and package['trials'][i]['accepted_export']
        assert trial['run']['provenance']['ground_sha256']==trial['adjusted_ground_sha256']
        for m in trial['adjusted_ground']['materials']:
            before=original[m['id']]
            if m['id'] not in selected:assert m==before;continue
            a=before['mohr_coulomb'];b=m['mohr_coulomb'];factor=trial['factor']
            assert np.isclose(b['cohesion_kpa'],a['cohesion_kpa']/factor,atol=1e-12,rtol=1e-12)
            assert np.isclose(np.tan(np.deg2rad(b['friction_deg'])),np.tan(np.deg2rad(a['friction_deg']))/factor,atol=1e-12,rtol=1e-12)
            assert b['young_kpa']==a['young_kpa'] and b['poisson']==a['poisson'] and m['unit_weight_kn_m3']==before['unit_weight_kn_m3']
        p=folder/f'trials/trial-{i:03}';vp=p/'independent.json';verify(p,vp)
        independent=json.loads(vp.read_text());maximum=max(f['relative_3d_displacement_error'] for f in independent['frames'])
        assert maximum<1e-9
        records.append({'factor':trial['factor'],'frames':len(independent['frames']),'maximum_independent_3d_displacement_error':maximum})
        files.extend([vp,p/'accepted-stage-run.json',p/'stage-results.json',p/'ground/ground.json'])
    report={'schema':'soil-fem-replay-verification/1','solver_sha256':digest,'unit_factor_exact':True,'stability_classification':'not_evaluated','trials':records,
        'scope':'All 3D elastic-history reference paths; plastic order independence and failure handling additionally tested in Rust',
        'sha256':{str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
    with output.open('x') as f:json.dump(report,f,indent=2)
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('folder',type=Path);p.add_argument('output',type=Path);a=p.parse_args();audit(a.folder,a.output)
