"""Seal the finite-bulk mixed-elastic scope after actual regression gates pass."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import math
import numpy as np
from freeze_sources import freeze
from verify_mcc_cycle_lifecycle import no_elapsed


def read(path):return json.loads(path.read_text(encoding='utf-8'))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def record(out):
    root=Path(__file__).resolve().parents[1]
    mixed=root/'runs/mixed-incompressibility-assessment-20261025-v1'
    column=root/'runs/mixed-column-verification-20261025-v1'
    curved=root/'runs/mixed-curved-cad-20261025-v1'
    exports=root/'runs/mixed-export-verification-20261025-v1'
    current=root/'runs/mixed-mechanics-current-mcc-20261025-v1'
    cycle=root/'runs/mcc-rate-general-development-20261024-v1'
    solver=current/'solver/soil-fem-3d.exe';binary_hash=sha(solver)
    metrics=read(mixed/'mixed-assessment-v2.json');assessment=read(mixed/'assessment.json');columns=read(column/'verification.json')
    curved_report=read(curved/'verification.json')
    exported=read(exports/'verification.json')
    assert exported['passed'] and len(exported['cases'])==3 and len(exported['rejected_inputs'])==5
    for name,digest in exported['evidence_sha256'].items():assert sha(exports/name)==digest,name
    for name,digest in exported['input_and_source_sha256'].items():assert sha(Path(name))==digest,name
    assert all(c.get('material_pressures_kept_separate',True) for c in exported['cases'])
    assert curved_report['passed'] and curved_report['all_native_integration_owners_positions_weights_verified']
    assert curved_report['pressure_coupling_nullity']==0
    for name,digest in curved_report['evidence_sha256'].items():assert sha(curved/name)==digest,name
    for name,digest in curved_report['input_and_source_sha256'].items():assert sha(Path(name))==digest,name
    curve_results=[read(curved/f'result-q{n}.json') for n in (3,4,5)]
    assert all(c['solver_binary_sha256']==binary_hash for c in curve_results)
    u=[np.asarray(c['displacement_m']) for c in curve_results]
    curve_changes=[float(np.linalg.norm(a-b)/np.linalg.norm(b)) for a,b in zip(u,u[1:])]
    assert curve_changes[1]<.6*curve_changes[0] and curve_changes[1]<1e-5
    for order in (8,9):
        independent=read(curved/f'independent-gl{order}.json')
        assert independent['passed'] and all(v['relative_l2']<1e-6 for v in independent['metrics'].values())
    assert metrics['passed'] and assessment['all_cases_solved'] and columns['passed']
    assert len(assessment['cases'])==9
    primary_pressure=[]
    for nu in (.3,.49,.4999):
        rows=sorted([r for r in metrics['cases'] if r['poisson']==nu],key=lambda r:r['subdivisions'])
        values=[r['primary_pressure_l2_over_sqrt_bulk']*math.sqrt(2*1000*(1+nu)/(3*(1-2*nu))) for r in rows]
        assert len(values)==3 and all(a>b for a,b in zip(values,values[1:])),(nu,values)
        primary_pressure.append(dict(poisson=nu,subdivisions=[r['subdivisions'] for r in rows],
            pressure_error_l2_kpa_m_to_3_over_2=values,exact_pressure_kpa=0.))
    for name,digest in columns['evidence_sha256'].items():assert sha(column/name)==digest,name
    for name,digest in assessment['evidence_sha256'].items():assert sha(mixed/name)==digest,name
    for path,digest in metrics['sha256'].items():assert sha(Path(path))==digest,path
    for name,digest in columns['input_and_source_sha256'].items():assert sha(Path(name))==digest,name
    for row in assessment['cases']:
        tag=f"tet10-n{row['subdivisions']}-nu{str(row['poisson']).replace('.','_')}"
        actual=read(mixed/(tag+'-result.json'));assert actual['solver_binary_sha256']==binary_hash
        assert actual['pressure_equation_relative_residual']<1e-9
    for row in columns['cases']:
        assert row['independent']['passed']
        assert row['maximum_displacement_error_m']<1e-10 and row['maximum_pressure_error_kpa']<1e-7
    # Bind the comparison to exactly the same original mesh and loads, including
    # the previously failed sparse solve and its bounded dense diagnostic.
    original=root/'runs/incompressibility-assessment-20261025-v2'
    for name in ('tet10-n4-mesh.json','tet10-n4-nu0_4999-case.json'):
        assert sha(original/name)==sha(mixed/name),name
    assert read(original/'dense-supplement-metrics.json')['displacement_relative_l2']>.39
    legacy=read(cycle/'slope-cycle-rust.json');native=read(current/'result.json')
    assert native['complete'] and native['provenance']['solver_binary_sha256']==binary_hash
    for key in ('frames','last_committed','reference_displacement_m'):assert native[key]==legacy[key],key
    for key in ('geostatic','attempts'):assert no_elapsed(native[key])==no_elapsed(legacy[key]),key
    for key in ('mesh_sha256','ground_sha256','case_sha256'):
        assert native['provenance'][key]==legacy['provenance'][key]
    frozen=read(current/'source-snapshot.json')
    for name,digest in frozen.items():
        if name.startswith('src/') or name in ('Cargo.toml','Cargo.lock','rust-toolchain.toml'):
            assert sha(root/name)==digest,name
    logs=['mixed-mechanics-tests-20261025-v1.log','mixed-mechanics-build-20261025-v1.log',
          'mixed-mechanics-release-tests-20261025-v1.log','mixed-mechanics-clippy-20261025-v1.log',
          'mixed-mechanics-fmt-20261025-v1.log','mixed-column-verification-20261025-v1.log',
          'mixed-incompressibility-assessment-20261025-v1.log','mixed-curved-cad-20261025-v1.log',
          'mixed-export-verification-20261025-v1.log']
    for name in logs:assert (root/'runs'/name).read_text(encoding='utf-8').rstrip().endswith('Exit code: 0'),name
    test_log=(root/'runs/mixed-mechanics-release-tests-20261025-v1.log').read_text(encoding='utf-8')
    counts=re.findall(r'test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored;',test_log)
    passed=sum(int(n) for n,_,_ in counts)
    assert passed>=210 and all(failed=='0' and ignored=='0' for _,failed,ignored in counts),counts
    out.mkdir(exist_ok=False);(out/'solver').mkdir();shutil.copyfile(solver,out/'solver/soil-fem-3d.exe')
    for name,source in (('mixed-data',mixed),('column-data',column),('current-mcc',current),('curved-cad-data',curved),('export-data',exports)):
        destination=out/name;destination.mkdir()
        for p in source.iterdir():
            if p.is_file():shutil.copyfile(p,destination/p.name)
            elif name in ('curved-cad-data','export-data'):shutil.copytree(p,destination/p.name)
    old=out/'original-displacement-diagnostic';old.mkdir()
    for name in ('tet10-n4-mesh.json','tet10-n4-nu0_4999-case.json','tet10-n4-nu0_4999.log',
                 'tet10-n4-nu0_4999-dense-case.json','tet10-n4-nu0_4999-dense-result.json',
                 'tet10-n4-nu0_4999-dense.log','dense-supplement-metrics.json'):
        shutil.copyfile(original/name,old/name)
    shutil.copyfile(cycle/'slope-cycle-rust.json',out/'current-mcc/legacy-result.json')
    (out/'logs').mkdir()
    for name in logs:shutil.copyfile(root/'runs'/name,out/'logs'/name)
    freeze(out)
    report=dict(schema='soil-fem-mixed-mechanics-final-validation/1',passed=True,goal_complete=False,p6_complete=False,
        scope='Small-strain finite-bulk Tet10/material-split P1 elasticity: manufactured cube, two-layer/rotated column and actual curved two-layer CAD mesh',
        solver_sha256=binary_hash,source_snapshot_sha256=sha(out/'source-snapshot.zip'),
        release_rust_tests_passed=passed,manufactured_cases=9,pressure_rank_checks=[26,124],
        primary_pressure_refinement=primary_pressure,
        analytic_column_and_independent_all_field_comparison=True,
        curved_cad_all_field_comparison_and_quadrature_refinement=True,curved_cad_displacement_refinements=curve_changes,
        material_split_primary_pressure_vtu_csv_export_verified=True,invalid_export_refusals=5,
        current_binary_legacy_mcc_all_history_and_nontiming_diagnostics_exact=True,
        accepted_mixed_assessment='mixed-data/mixed-assessment-v2.json',
        historical_metric_note='mixed-assessment.json is retained history; v2 corrects the dimensional label of the pointwise L2 volume-relation diagnostic',
        limitations=['No general inf-sup proof for arbitrary terrain/material configurations or exact incompressibility',
                     'Mixed MCC/MC, construction, history restoration and native SolidCraft UI integration remain required',
                     'The separately archived 5 kPa MCC refinement uses its retained production binary; its acceptance is distinct from mixed elasticity',
                     'All other unfinished P0-P6 requirements remain unchanged'],
        sha256={p.relative_to(root).as_posix():sha(p) for p in out.rglob('*') if p.is_file()})
    with (out/'final-validation.json').open('x') as f:json.dump(report,f,indent=2)
    print(json.dumps({k:report[k] for k in ('passed','goal_complete','scope','release_rust_tests_passed','manufactured_cases')}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);a=p.parse_args();record(a.output.resolve())
