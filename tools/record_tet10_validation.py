"""Bind separately executed checks and integration results to current sources."""
import argparse
import hashlib
import json
from pathlib import Path
import re


def record(root,run,cad,rust_log,lint_log,format_log):
    suite=json.loads((run/'verification.json').read_text());spatial=json.loads((run/'spatial-sensitivity.json').read_text())
    solver=root/'target/release/soil-fem-3d.exe';digest=hashlib.sha256(solver.read_bytes()).hexdigest()
    assert suite['passed'] and suite['solver_sha256']==digest
    rust=rust_log.read_text();assert 'Exit code: 0' in rust and '0 failed' in rust
    count=sum(map(int,re.findall(r'test result: ok\. (\d+) passed;',rust)));assert count==54
    for log in (lint_log,format_log):assert log.read_text().rstrip().endswith('Exit code: 0')
    cad_log=cad/'cad-tests.log';text=cad_log.read_text();assert 'Ran 5 tests' in text and 'OK' in text and text.rstrip().endswith('Exit code: 0')
    checkpoint=json.loads((run/'checkpoint/checkpoint-verification.json').read_text());assert checkpoint['solver_sha256']==digest
    assert checkpoint['exact_final_committed_state'] and checkpoint['exact_frames'] and checkpoint['exact_increment_path']
    local=json.loads((run/'local-reference/independent-local-verification.json').read_text());assert local['solver_sha256']==digest
    for name in ('flat','curved','curved-fine'):
        package=json.loads((cad/name/'accepted-stage-run.json').read_text());assert package['solver_binary_sha256']==digest
        for path,expected in package['inputs_and_results_sha256'].items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==expected
    files=[solver,root/'Cargo.lock',root/'rust-toolchain.toml',root/'sources.lock.json',rust_log,lint_log,format_log,cad_log,
        run/'verification.json',run/'spatial-sensitivity.json',run/'spatial.log',run/'checkpoint/checkpoint-verification.json',
        run/'local-reference/independent-local-verification.json',run/'high-order/verification.json',run/'independent-plastic/verification.json']
    files.extend(p for folder in ('src','tools','scripts','tests','examples') for p in (root/folder).rglob('*')
        if p.is_file() and p.suffix in ('.rs','.py','.ps1','.json'))
    inventory={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    report={'schema':'soil-fem-tet10-final-validation/1','passed':True,'rust_tests':count,'cad_geometry_tests':5,
        'solver_sha256':digest,'cad_directory':str(cad),'suite_directory':str(run),
        'scope':'Local Tet10 excavation, general boundary and MCC integration. P6 remains incomplete.',
        'independent_3d_mcc_comparison':'elastic branch; plastic hardening has separate 1D theory and full-history recovery tests',
        'spatial_convergence_proven':False,'floor_area_mean_relative_change':spatial['floor_area_mean_relative_change'],
        'sha256':inventory}
    with (run/'final-validation.json').open('x',encoding='utf-8') as f:json.dump(report,f,indent=2)
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--cad',type=Path,required=True)
    p.add_argument('--rust-log',type=Path,required=True);p.add_argument('--lint-log',type=Path,required=True);p.add_argument('--format-log',type=Path,required=True)
    a=p.parse_args();record(Path(__file__).resolve().parents[1],*[v.resolve() for v in (a.run,a.cad,a.rust_log,a.lint_log,a.format_log)])
