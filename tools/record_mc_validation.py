"""Bind independently run MC integration checks and source/binary inventories."""
import argparse
import hashlib
import json
from pathlib import Path
import re


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def record(run,rust_log,lint_log,format_log):
    root=Path(__file__).resolve().parents[1];suite=json.loads((run/'verification.json').read_text());solver=run/'solver/soil-fem-3d.exe'
    digest=sha(solver);assert suite['passed'] and suite['solver_sha256']==digest==sha(root/'target/release/soil-fem-3d.exe')
    text=rust_log.read_text();assert text.rstrip().endswith('Exit code: 0')
    count=sum(map(int,re.findall(r'test result: ok\. (\d+) passed;',text)));assert count==63
    for p in (lint_log,format_log):assert p.read_text().rstrip().endswith('Exit code: 0')
    cad=run/'cad-tests.log';text=cad.read_text();assert 'Ran 5 tests' in text and 'OK' in text and text.rstrip().endswith('Exit code: 0')
    paths=json.loads((run/'material-paths/verification.json').read_text());assert paths['solver_sha256']==digest and len(paths['paths'])==8 and len(paths['observed_return_modes'])==9
    column=json.loads((run/'columns/verification.json').read_text());assert column['solver_sha256']==digest and column['separate_process_tet10_plastic_restart_exact']
    for name in ('mc-cad','mcc-cad'):
        package=json.loads((run/name/'accepted-stage-run.json').read_text());assert package['solver_binary_sha256']==digest
        assert sha(run/name/'solver/soil-fem-3d.exe')==digest
        for path,expected in package['inputs_and_results_sha256'].items():assert sha(Path(path))==expected
    regression=json.loads((run/'mcc-regression.json').read_text());assert regression['exact_committed_state'] and regression['exact_frames']
    failed=root/'runs/mc-local-first-20261023';failure=json.loads((failed/'stage-results.json').read_text())
    assert not failure['complete'] and not (failed/'accepted-stage-run.json').exists() and not (failed/'results').exists()
    assert all(a['committed_state_unchanged_on_rejection'] for a in failure['attempts'] if not a['accepted'])
    assert failure['provenance']['solver_binary_sha256']==sha(failed/'solver/soil-fem-3d.exe')
    files=[solver,rust_log,lint_log,format_log,cad,run/'verification.json',run/'mcc-regression.json',run/'mcc-theory/verification.json',
        run/'material-paths/verification.json',run/'columns/verification.json',failed/'stage-results.json',failed/'solver/soil-fem-3d.exe',
        root/'Cargo.lock',root/'sources.lock.json',root/'requirements.lock.txt',root/'rust-toolchain.toml']
    files.extend(p for folder in ('src','tools','scripts','tests','examples') for p in (root/folder).rglob('*')
        if p.is_file() and p.suffix in ('.rs','.py','.ps1','.json'))
    report={'schema':'soil-fem-mc-final-validation/1','passed':True,'rust_tests':count,'cad_tests':5,
        'independent_material_paths':8,'return_modes_observed':9,'typed_tet10_plastic_restart_exact':True,
        'mcc_compatible_exact':True,'mc_local_excavation_converged':True,
        'old_weaker_mc_case':'Numerical nonconvergence retained without accepted export; no physical collapse/Fs inference',
        'scope':'MC/MCC constitutive and staged FEM milestone; general 3D localization/uniqueness, 2D execution comparison and P6 overall remain incomplete',
        'solver_sha256':digest,'sha256':{str(p.relative_to(root)):sha(p) for p in files}}
    with (run/'final-validation.json').open('x',encoding='utf-8') as f:json.dump(report,f,indent=2)
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--rust-log',type=Path,required=True)
    p.add_argument('--lint-log',type=Path,required=True);p.add_argument('--format-log',type=Path,required=True)
    a=p.parse_args();record(*[v.resolve() for v in (a.run,a.rust_log,a.lint_log,a.format_log)])
