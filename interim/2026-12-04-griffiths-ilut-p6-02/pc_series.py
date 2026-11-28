"""ILUT-only wrapper around the unchanged pinned 661-node comparator gates."""
import argparse
import copy
from pathlib import Path
import tempfile
from prepare import ROOT, read, write
from study import sha
from prepare_slope_h_series import prepare as original_prepare
from griffiths_2d import (base_report, displacement_metrics, independent_audit,
    load_reference, sample_nodes, validate_run)


def prepare(out, reference):
    manifest=original_prepare(out,1,reference)
    (out/'preparation.json').rename(out/'original-preparation.json')
    for record in manifest['levels']:
        sub=out/record['folder']
        for name in ('case.json','stages.json'):
            case=read(sub/name)
            case['settings']['linear_solver']['options']['ilut']=dict(drop_tolerance=1e-4,fill_limit=200)
            # Only newly generated files are changed. The source manifest is retained.
            import json
            (sub/name).write_text(json.dumps(case,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
        record['files_sha256']={name:sha(sub/name) for name in record['files_sha256']}
    manifest.update(schema='soil-fem-griffiths-ilut-preparation/1',
        original_preparation_sha256=sha(out/'original-preparation.json'),
        ilut=dict(drop_tolerance=1e-4,fill_limit=200),
        note='Only explicit ILUT options differ from regenerated pinned original inputs; not a solve or Fs.')
    write(out/'preparation-ilut.json',manifest)
    return manifest


def preflight(folder,reference):
    manifest=read(folder/'preparation-ilut.json')
    with tempfile.TemporaryDirectory(prefix='ilut-pinned-inputs-') as temp:
        expected=prepare(Path(temp)/'expected',reference)
        if manifest!=expected: raise ValueError('ILUT preparation manifest changed')
        if sha(folder/'original-preparation.json')!=expected['original_preparation_sha256']:
            raise ValueError('original preparation changed')
        for record in expected['levels']:
            for name,h in record['files_sha256'].items():
                if sha(folder/record['folder']/name)!=h: raise ValueError('prepared input changed: '+name)
    return manifest


def compare(folder,reference_path,solver,output):
    reference=load_reference(reference_path);manifest=preflight(folder,reference_path)
    binary_sha=sha(solver);levels=[];files={}
    # All levels must be valid before any independent audit is run.
    for record in manifest['levels']:
        sub=folder/record['folder']
        mesh,ground,case,run=[read(sub/(name+'.json')) for name in ('mesh','ground','case','result')]
        u=validate_run(run,mesh,ground,case,binary_sha)
        mapping=read(sub/'reference-map.json')
        if mapping!=sample_nodes(mesh,reference): raise ValueError('original node map changed')
        samples={depth:u[indices][:,[0,2]] for depth,indices in mapping.items()}
        levels.append(dict(folder=record['folder'],elements=record['elements'],
                           **displacement_metrics(reference['displacement_xz_m'],samples)))
        files.update({str(sub/name):sha(sub/name) for name in (*record['files_sha256'],'result.json')})
    for record in manifest['levels']: independent_audit(folder/record['folder'],solver)
    if sha(solver)!=binary_sha or any(sha(Path(p))!=h for p,h in files.items()):
        raise ValueError('binary/input/result changed during audit')
    result=dict(schema='soil-fem-griffiths-ilut-661-comparison/1',**base_report(reference_path),levels=levels,
        new_slope_solve_performed=False,measurement_source='separately completed user-PC production solves',
        independent_audit_passed=True,all_levels_five_percent_passed=all(x['five_percent_verdict']=='pass' for x in levels),
        continuum_verdict='indeterminate',solver_binary_sha256=binary_sha,files_sha256=files)
    write(output,result)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['prepare','preflight','compare'])
    p.add_argument('folder',type=Path);p.add_argument('--reference',type=Path,default=ROOT/'tests/data/griffiths_2d_reference.json')
    p.add_argument('--solver',type=Path);p.add_argument('--output',type=Path);a=p.parse_args()
    if a.mode=='prepare': prepare(a.folder,a.reference)
    elif a.mode=='preflight': preflight(a.folder,a.reference)
    else:
        if not a.solver or not a.output: p.error('compare requires --solver and --output')
        report=compare(a.folder,a.reference,a.solver,a.output)
        raise SystemExit(0 if report['all_levels_five_percent_passed'] else 1)
