"""Export only completed replays; retain incomplete and rejected diagnostics."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
from export_stages import export
from verify_stage_package import verify


def capture(mesh,result,output,solver):
    report=json.loads(result.read_text());assert report['schema']=='soil-fem-replay-result/1'
    assert report['stability_classification']=='not_evaluated'
    digest=hashlib.sha256(solver.read_bytes()).hexdigest();assert report['source_identity']['solver_binary_sha256']==digest
    output.mkdir(parents=True,exist_ok=False);records=[]
    for i,trial in enumerate(report['trials']):
        folder=output/f'trial-{i:03}';folder.mkdir()
        with (folder/'trial.json').open('x') as f:json.dump({k:v for k,v in trial.items() if k not in ('run','adjusted_ground')},f,indent=2)
        row={'factor':trial['factor'],'outcome':trial['outcome'],'folder':folder.name,'accepted_export':False}
        if trial['run'] is not None:
            run=trial['run'];assert run['provenance']['solver_binary_sha256']==digest and run['ground']==trial['adjusted_ground']
            (folder/'ground').mkdir();shutil.copyfile(mesh,folder/'ground/mesh.json')
            for p,v in ((folder/'ground/ground.json',run['ground']),(folder/'stages.json',run['case']),(folder/'stage-results.json',run)):
                with p.open('x') as f:json.dump(v,f,indent=2)
            if run['complete']:
                assert trial['outcome']=='completed'
                export(folder/'ground/mesh.json',folder/'stage-results.json',folder/'results')
                verify(folder/'ground/mesh.json',folder/'ground/ground.json',folder/'stage-results.json',folder/'results',solver=solver)
                row['accepted_export']=True
            else:assert trial['outcome']=='incomplete' and not (folder/'results').exists()
        records.append(row)
    package={'schema':'soil-fem-replay-export/1','solver_sha256':digest,'source_report_sha256':hashlib.sha256(result.read_bytes()).hexdigest(),
        'source_mesh_sha256':hashlib.sha256(mesh.read_bytes()).hexdigest(),'stability_classification':'not_evaluated','trials':records}
    with (output/'verification.json').open('x') as f:json.dump(package,f,indent=2)
    print(json.dumps(package,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mesh',type=Path);p.add_argument('result',type=Path);p.add_argument('output',type=Path);p.add_argument('--solver',type=Path,required=True)
    a=p.parse_args();capture(a.mesh,a.result,a.output,a.solver)
