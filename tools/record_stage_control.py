"""Seal the external staged-mechanics control foundation, not native UI/P6."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import zipfile
from freeze_sources import freeze


def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def receipt(path):
    r=read(path);assert r['passed'],path
    for name,digest in r.get('input_and_source_sha256',{}).items():assert sha(Path(name))==digest,name
    for name,digest in r.get('evidence_sha256',{}).items():assert sha(path.parent/name)==digest,name
    return r


def record(out):
    root=Path(__file__).resolve().parents[1];runs=root/'runs';development=runs/'stage-control-development-20261025-v1'
    binary=development/'solver/soil-fem-3d.exe';digest=sha(binary)
    assert digest=='59207c2cd184562ab6a4419054cea00d59ba67be3c468bdcc3597a93a86d3530'
    cad=runs/'stage-job-cad-verification-20261025-v2';cr=receipt(cad/'verification.json')
    exported=runs/'stage-control-mixed-export-20261025-v1';er=receipt(exported/'verification.json')
    legacy=runs/'stage-control-legacy-20261025-v1';lr=receipt(legacy/'verification.json')
    mixed=receipt(development/'mixed-compatibility.json')
    assert cr['solver_sha256']==lr['solver_sha256']==digest
    assert cr['live_newton_observed_before_stop_request'] and cr['plastic_history_present'] and cr['all_nonruntime_fields_exact']
    assert len(er['rejected_inputs'])==15 and len(mixed['cases'])==2
    paused=development/'terrain-job';status=read(paused/'status.json');checkpoint=read(paused/'checkpoint.json')
    audit=read(paused/'material-strain-audit.json')
    assert status['state']=='checkpointed' and status['artifact_sha256']==sha(paused/'checkpoint.json')
    assert checkpoint['run']['failure_reason']=='checkpoint_requested' and checkpoint['run']['geostatic']['progress']==0.140625
    assert audit['audit_completed'] and audit['initial_state_within_configured_strain_domain'] and not audit['native_run_complete']
    assert abs(audit['maximum_material_engineering_strain_from_seed']-checkpoint['run']['geostatic']['maximum_absolute_engineering_strain_from_seed'])<1e-10
    for name,value in audit['input_and_source_sha256'].items():assert sha(Path(name))==value,name
    logs=['stage-control-job-material-tests-20261025-v1.log','mixed-material-domain-unit-20261025-v1.log',
          'stage-control-clippy-20261025-v1.log','stage-control-fmt-20261025-v1.log','stage-control-release-tests-20261025-v1.log',
          'stage-job-cad-verification-20261025-v2.log','stage-control-mixed-export-20261025-v1.log',
          'stage-control-legacy-20261025-v1.log','stage-control-mixed-compatibility-20261025-v1.log']
    for name in logs:assert (runs/name).read_text(encoding='utf-8').rstrip().endswith('Exit code: 0'),name
    counts=re.findall(r'test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored',(runs/logs[4]).read_text())
    assert sum(int(n) for n,_,_ in counts)==238 and all(f==i=='0' for _,f,i in counts)
    build=development/'build-sources';sources=read(build/'source-snapshot.json');assert len(sources)==478
    with zipfile.ZipFile(build/'source-snapshot.zip') as z:
        assert set(z.namelist())==set(sources)
        assert all(hashlib.sha256(z.read(k)).hexdigest()==v for k,v in sources.items())
    compiled={k:v for k,v in sources.items() if k.startswith(('src/','tests/','examples/')) or k in ('Cargo.toml','Cargo.lock','rust-toolchain.toml')}
    assert all(sha(root/k)==v for k,v in compiled.items()),'current compiled sources differ from this build snapshot'
    out.mkdir(exist_ok=False);(out/'solver').mkdir();shutil.copyfile(binary,out/'solver/soil-fem-3d.exe')
    for name in ('source-snapshot.json','source-snapshot.zip'):shutil.copyfile(build/name,out/name)
    for tag,source in [('actual-cad-job',cad),('export',exported),('legacy',legacy)]:shutil.copytree(source,out/tag)
    columns=out/'mixed-cases';columns.mkdir();shutil.copytree(development/'n4',columns/'n4')
    for name in ('ground.json','two-material-example.json','mixed-compatibility.json'):shutil.copyfile(development/name,columns/name)
    # Only the terminal paused job is frozen. The resumed terrain job is still
    # live and is explicitly outside this acceptance evidence inventory.
    shutil.copytree(paused,out/'bounded-terrain-checkpoint')
    failed=out/'not-accepted-previous-terrain';failed.mkdir()
    old=runs/'mixed-mcc-ordering-resolution-20261025-v1/terrain'
    for name in ('material-strain-audit.json','case.json','ground.json','mesh.json','result.json'):shutil.copyfile(old/name,failed/name)
    (out/'logs').mkdir()
    for name in logs:shutil.copyfile(runs/name,out/'logs'/name)
    v=out/'verification-sources';v.mkdir();freeze(v)
    report=dict(schema='soil-fem-stage-control-foundation-final-validation/1',passed=True,goal_complete=False,p6_complete=False,
        scope='External staged-mechanics CLI: immutable job, observed Newton/linear/attempt data, accepted-boundary pause and separate-process resume; initial mixed-material strain-domain calculation/guard',
        solver_sha256=digest,source_snapshot_sha256=sha(out/'source-snapshot.zip'),build_source_files=478,
        verification_source_snapshot_sha256=sha(v/'source-snapshot.zip'),release_rust_tests_passed=238,
        observer_and_cli_tests=9,mixed_strain_formula_and_cumulative_domain_unit_test=True,
        real_cad_live_plastic_construction_pause_resume_exact=True,legacy_mcc_mc_results_exact=True,
        mixed_columns_primary_fields_and_histories_exact=True,new_material_strain_cache_recomputed=True,
        invalid_export_refusals=15,terminal_terrain_checkpoint_within_initial_strain_domain=True,
        bounded_checkpoint_progress=0.140625,maximum_material_engineering_strain_at_checkpoint=audit['maximum_material_engineering_strain_from_seed'],
        native_ui_stop_resume_connected=False,general_mixed_terrain_accepted=False,
        limitations=['Native host/worker Stop still uses its previous process-stop path until connected to this protocol',
                     'ONE_WAY and Richards-specific job controllers are not wrapped by these new StageCase CLI commands',
                     'The complete general mixed-MCC terrain run remains unresolved; a bounded checkpoint is not a complete equilibrium',
                     'All other unfinished P0-P6 requirements remain required'],
        sha256={p.relative_to(root).as_posix():sha(p) for p in out.rglob('*') if p.is_file()})
    with (out/'final-validation.json').open('x',encoding='utf-8') as stream:json.dump(report,stream,indent=2,allow_nan=False)
    print(json.dumps({k:report[k] for k in ('passed','goal_complete','scope','release_rust_tests_passed','native_ui_stop_resume_connected','general_mixed_terrain_accepted')}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);record(p.parse_args().output.resolve())
