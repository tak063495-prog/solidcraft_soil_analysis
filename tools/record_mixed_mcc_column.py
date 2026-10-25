"""Seal the dry mixed-MCC column scope, keeping terrain/P6 explicitly open.

The solver's historical build sources are the primary archive. Current sources
are archived separately for verification tools; they are not mislabeled as the
source of a previously built executable.
"""
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


def verify_receipt(path):
    report=read(path);assert report['passed'],path
    for name,digest in report.get('input_and_source_sha256',{}).items():assert sha(Path(name))==digest,name
    for name,digest in report.get('evidence_sha256',{}).items():assert sha(path.parent/name)==digest,name
    return report


def record(out):
    root=Path(__file__).resolve().parents[1];runs=root/'runs'
    study=runs/'mixed-mcc-ordering-resolution-20261025-v1'
    origin=runs/'mixed-mcc-implicit-resolution-20261025-v1'
    solver=study/'solver/soil-fem-3d.exe';digest=sha(solver)
    assert digest=='99d154374e9b83f9dbf8971f49089faf3b079ed6dfe2a86863f86b982e770f08'
    metrics=verify_receipt(study/'resolution-assessment.json')
    reference=verify_receipt(study/'reference-refinement-v2.json')
    assert metrics['candidate_resolution']==32 and not metrics['all_coarse_and_fine_field_targets_passed']
    assert metrics['passed'] and reference['passed']
    fine=next(c for c in metrics['comparisons'] if c['vertical_subdivisions']==32)
    assert len(fine['frames'])==4
    assert all(c['integration_layout_checked_independently'] for c in metrics['comparisons'])
    coarse_failures=[g for g in metrics['gates'] if not g['passed']]
    assert len(coarse_failures)==3 and all(g['name'].startswith('n16:') and g['name'].endswith(':pc') for g in coarse_failures)
    receipts={}
    for tag,folder in [('lifecycle','mixed-mcc-ordering-lifecycle-20261025-v1'),
                       ('export','mixed-mcc-ordering-export-verification-20261025-v1'),
                       ('legacy','mixed-mcc-ordering-legacy-20261025-v1')]:
        receipts[tag]=(runs/folder,verify_receipt(runs/folder/'verification.json'))
    assert len(receipts['lifecycle'][1]['restarts'])==3 and receipts['lifecycle'][1]['rejected_converged_trials']==3
    assert len(receipts['export'][1]['rejected_inputs'])==14
    assert receipts['legacy'][1]['solver_sha256']==digest
    verify_receipt(study/'first-increment-comparison.json')
    native={n:read(study/f'n{n}/native-ultrafine.json') for n in (16,32)}
    for n,d in native.items():
        assert d['complete'] and len(d['frames'])==4 and d['provenance']['solver_binary_sha256']==digest
        assert sum(a['accepted'] for a in d['attempts'])==384
        assert all(a['accepted'] for a in d['attempts'])
    failed=read(study/'terrain/result.json')
    assert not failed['complete'] and failed['geostatic']['progress']<1.
    assert 'pivot failure' in failed['failure_reason']
    logs=['mixed-mcc-ordering-tests-20261025-v2.log','mixed-mcc-ordering-material-global-tests-20261025-v1.log',
          'mixed-mcc-ordering-clippy-20261025-v1.log','mixed-mcc-ordering-fmt-20261025-v1.log',
          'mixed-mcc-ordering-release-tests-20261025-v1.log','mixed-mcc-ordering-lifecycle-20261025-v1.log',
          'mixed-mcc-ordering-export-verification-20261025-v1.log','mixed-mcc-ordering-legacy-20261025-v1.log',
          'mixed-mcc-ordering-reference-refinement-20261025-v2.log','mixed-mcc-ordering-resolution-assessment-20261025-v1.log',
          'mixed-mcc-ordering-native-n16-resumed-20261025-v1.log','mixed-mcc-ordering-native-n32-resumed-20261025-v1.log']
    for name in logs:assert (runs/name).read_text(encoding='utf-8').rstrip().endswith('Exit code: 0'),name
    testlog=(runs/'mixed-mcc-ordering-release-tests-20261025-v1.log').read_text()
    counts=re.findall(r'test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored',testlog)
    assert sum(int(n) for n,_,_ in counts)==228 and all(f==i=='0' for _,f,i in counts)
    build=study/'build-sources';sources=read(build/'source-snapshot.json');assert len(sources)==470
    with zipfile.ZipFile(build/'source-snapshot.zip') as z:
        assert set(z.namelist())==set(sources)
        assert all(hashlib.sha256(z.read(k)).hexdigest()==v for k,v in sources.items())
    out.mkdir(exist_ok=False);(out/'solver').mkdir();shutil.copyfile(solver,out/'solver/soil-fem-3d.exe')
    for name in ('source-snapshot.json','source-snapshot.zip'):shutil.copyfile(build/name,out/name)
    data=out/'column-data';data.mkdir()
    shutil.copyfile(study/'ground.json',data/'ground.json')
    for n in (4,16,32):shutil.copytree(study/f'n{n}',data/f'n{n}')
    for name in ('resolution-assessment.json','assessment-n16.json','reference-refinement-v2.json','first-increment-comparison.json','reference-reuse.json'):
        shutil.copyfile(study/name,data/name)
    shutil.copytree(origin/'reference-sources',out/'reference-sources')
    shutil.copyfile(origin/'preparation.json',out/'reference-sources/original-preparation.json')
    for tag,(folder,_) in receipts.items():shutil.copytree(folder,out/tag)
    diagnostic=out/'not-accepted-terrain';diagnostic.mkdir()
    for name in ('mesh.json','ground.json','case.json','result.json'):shutil.copyfile(study/'terrain'/name,diagnostic/name)
    shutil.copyfile(runs/'mixed-mcc-ordering-terrain-resumed-20261025-v1.log',diagnostic/'solver.log')
    (out/'logs').mkdir()
    for name in logs:shutil.copyfile(runs/name,out/'logs'/name)
    verification=out/'verification-sources';verification.mkdir();freeze(verification)
    peak=lambda field,key:max(r['fields'][field][key] for r in fine['frames'])
    report=dict(schema='soil-fem-mixed-mcc-column-final-validation/1',passed=True,goal_complete=False,p6_complete=False,
        scope='Dry small-strain NC two-layer confined Tet10/P1/GP-volume MCC column: K0=1, pressure 100->130->100 kPa, upper-half excavation with retained 100 kPa',
        solver_sha256=digest,source_snapshot_sha256=sha(out/'source-snapshot.zip'),build_source_files=470,
        primary_source_snapshot_meaning='Exact historical sources of solver 99d15437; current-tree control/strain-guard additions are outside this accepted executable',
        verification_source_snapshot_sha256=sha(verification/'source-snapshot.zip'),
        release_rust_tests_passed=228,finest_vertical_subdivisions=32,finest_tet10_elements=192,finest_material_points=5184,
        independent_reference_steps=[400,800],reference_initialized_from_rust_results=False,
        original_field_targets_unchanged=True,coarse_n16_pc_failures_retained=coarse_failures,
        maximum_displacement_relative_error=max(r['displacement_relative_all_nodes'] for r in fine['frames']),
        maximum_displacement_component_error_m=max(r['displacement_maximum_component_error_m'] for r in fine['frames']),
        maximum_stress_relative_error=peak('stress','relative_all_points'),maximum_pc_error_kpa=peak('pc','maximum_absolute_all_points'),
        maximum_void_ratio_error=peak('void_ratio','maximum_absolute_all_points'),maximum_plastic_strain_error=peak('plastic_strain','maximum_absolute_all_points'),
        floor_excavation_increment=fine['frames'][-1]['monitors']['floor'],
        lifecycle_scope='n4: 3 separate-process midpoint resumes and 3 converged trial rejections; n16/n32: actual separate-process resume from the first accepted load increment',
        primary_pressure_and_all_gp_history_exports_verified=True,legacy_mcc_mc_nonruntime_results_exact=True,
        general_terrain_accepted=False,terrain_failure_progress=failed['geostatic']['progress'],
        limitations=['No general terrain/material-layout stability acceptance; the stored slope stops during geostatic initialization',
                     'Later small-strain-domain checks and cooperative job control need their own solver/source verification',
                     'No wet, structural, JOINT, birth/material-transition or native UI acceptance for mixed MCC',
                     'All other unfinished P0-P6 requirements remain required'],
        sha256={p.relative_to(root).as_posix():sha(p) for p in out.rglob('*') if p.is_file()})
    with (out/'final-validation.json').open('x',encoding='utf-8') as f:json.dump(report,f,indent=2,allow_nan=False)
    print(json.dumps({k:report[k] for k in ('passed','goal_complete','scope','release_rust_tests_passed','maximum_pc_error_kpa','general_terrain_accepted')}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);record(p.parse_args().output.resolve())
