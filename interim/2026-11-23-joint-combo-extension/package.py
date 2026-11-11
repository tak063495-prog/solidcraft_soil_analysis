"""Seal new A-01 evidence; never overwrite historical or failed evidence."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import zipfile
from fixtures import ROOT,read,write


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def seal(matrix,plate,refusals,stacked,diagnostics,retained,destination):
    destination.mkdir(parents=True,exist_ok=True)
    selected={p.name:p for p in matrix.iterdir() if p.is_dir() and (p/'comparison.json').is_file()}
    selected['joint-plate']=plate/'joint-plate';assert len(selected)==7
    measurements={};restarts=[]
    for name,folder in sorted(selected.items()):
        assert read(folder/'native/run.json')['complete'] and read(folder/'independent-result.json')['complete']
        native=read(folder/'native-verification.json');comparison=read(folder/'comparison.json')
        assert native['accepted_states']==comparison['accepted_pairs']==96
        result=read(folder/'restarts/verification.json');assert result['passed'] and len(result['checks'])==3
        restarts.extend(dict(case=name,**row) for row in result['checks'])
        measurements[name]=dict(verdict='未達',scope='dry non-tip elastic MC bulk; fixed input increment 1/8; contact BIRTH/death/rebirth and complete boundary cycle',
            fresh_native_solve=True,fresh_independent_solve=True,**comparison)
    diagnostic=read(diagnostics);assert set(diagnostic)==set(measurements)
    refusals_report=read(refusals/'refusals/verification.json');assert refusals_report['passed'] and len(refusals_report['checks'])==18
    report=dict(schema='soil-fem-a01-combination-extension/1',base='c70fb37ea025b6e40dda4db4adf67379ea2a6ee3',
        A01='未達',passed=False,requested_all_combinations_passed=False,new_limited_acceptance=False,
        label='2026-11-23',executed_jst='2026-11-12',selected_fresh_native_paths=7,selected_fresh_independent_paths=7,
        compared_accepted_pairs=sum(x['accepted_pairs'] for x in measurements.values()),
        comparison_method='all accepted input grid points; all nodes/material GP/contact GP; absolute maxima and relative pooled L2. Fields derived from actual saved values are distinguished from actual internal/target diagnostics. Saved-state/local audits are not fresh 3D solves.',
        cases=measurements,actual_checkpoint_restarts=dict(passed=True,count=len(restarts),
            suffix_states_exact=sum(x['suffix_states_exact'] for x in restarts),whole_run_exact_except_elapsed=True,checks=restarts),
        scope_and_mutation_refusals=refusals_report,interacting_stacked_attempt=read(stacked/'verification.json'),
        remaining=dict(NC01='未達・対象外',P6='未達',A02_A04='未達・対象外',B01_B03='未達・対象外',C01_C04='未達・対象外'),
        superseded_plate_observation='Main matrix plate comparison detected absent Kirchhoff-zero transverse shear observations; retained as failed. Corrected plate case was freshly solved on both sides, including three fresh-process resumes. Not relabelled saved re-audit as new solving.')
    write(destination/'measurements.json',report);write(destination/'cause-diagnostics.json',diagnostic)
    write(destination/'refusals.json',refusals_report)
    # Keep every measured/failed preliminary folder, including partial output.
    sources={}
    for name,folder in selected.items():sources['selected/'+name]=folder
    sources['main-matrix-including-failed-plate']=matrix
    sources['corrected-plate-matrix']=plate
    sources['scope-and-mutation-refusals']=refusals
    sources['interacting-stacked-attempt']=stacked
    for folder in retained:sources['retained-failures/'+folder.name]=folder
    archive=ROOT/'work/a01-new-evidence.zip';assert not archive.exists()
    manifest={}
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for prefix,folder in sorted(sources.items()):
            for path in sorted(folder.rglob('*')):
                if path.is_file() and '__pycache__' not in path.parts and '_bin' not in path.parts:
                    name=prefix+'/'+path.relative_to(folder).as_posix()
                    z.write(path,name);manifest[name]=dict(size=path.stat().st_size,sha256=sha(path))
    parts=[]
    with archive.open('rb') as stream:
        while chunk:=stream.read(24*1024*1024):
            path=destination/f'evidence.zip.part{len(parts):03}'
            with path.open('xb') as part:part.write(chunk)
            parts.append(dict(name=path.name,size=len(chunk),sha256=sha(path)))
    write(destination/'evidence-manifest.json',dict(archive_sha256=sha(archive),archive_bytes=archive.stat().st_size,
        parts=parts,entries=manifest,solver_binaries_archived=False,all_raw_states_inputs_checkpoints_and_failed_runs_retained=True))
    print('sealed',len(manifest),'entries',archive.stat().st_size,'bytes',len(parts),'parts')


def main():
    p=argparse.ArgumentParser();p.add_argument('--matrix',type=Path,required=True);p.add_argument('--plate',type=Path,required=True)
    p.add_argument('--refusals',type=Path,required=True);p.add_argument('--stacked',type=Path,required=True)
    p.add_argument('--diagnostics',type=Path,required=True);p.add_argument('--retained',type=Path,action='append',default=[])
    p.add_argument('--destination',type=Path,required=True);a=p.parse_args()
    seal(a.matrix,a.plate,a.refusals,a.stacked,a.diagnostics,a.retained,a.destination)


if __name__=='__main__':main()
