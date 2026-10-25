"""Bind full new-engine core and finite-width steady validation to one freeze."""
import argparse
import json
from pathlib import Path
from freeze_sources import manifest
from verify_acceptance_packet import verify as packet,sha


def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def record(run):
    root=Path(__file__).resolve().parents[1];assert not (run/'final-validation.json').exists()
    assert read(run/'source-snapshot.json')==manifest(),'source mutation during acceptance'
    core=packet(run/'core');digest=sha(run/'solver/soil-fem-3d.exe');assert digest==core['solver_sha256']==sha(root/'target/release/soil-fem-3d.exe')
    for n in ('core','reference-tests','dam','counterexample'):assert (run/f'{n}.log').read_text(encoding='utf-8').rstrip().endswith('Exit code: 0'),n
    assert 'Ran 2 tests' in (run/'reference-tests.log').read_text(encoding='utf-8')
    core_report=read(run/'core/final-validation.json');assert core_report['rust_tests']>=165 and core_report['separate_process_restarts']==8 and core_report['exact_legacy_cases']==6
    assert 'finite_width_reservoir_converges_without_changing_local_mass_or_atmospheric_gates ... ok' in (run/'core/rust-tests.log').read_text(encoding='utf-8')
    dam=read(run/'dam/verification.json');counter=read(run/'counterexample/verification.json')
    assert dam['passed'] and counter['passed'] and dam['solver_sha256']==counter['solver_sha256']==digest
    assert len(dam['comparisons'])==8 and len(dam['independent_reference_refinements'])==2
    for label in ('full','strip'):
        for n in (4,8,16,24):
            result=read(run/f'dam/rt0-{label}-{n}/result.json');assert result['complete'] and result['solver_binary_sha256']==digest
            assert read(run/f'dam/rt0-{label}-{n}/verification.json')['passed']
        for n in (16,32,64):assert read(run/f'dam/fv-{label}-{n}.json')['maximum_cell_mass_residual_m3_s']<1e-17
    old_manifest=root/counter['old_parent_manifest'];assert sha(old_manifest)==counter['old_parent']['manifest_sha256']
    files=[p for p in run.rglob('*') if p.is_file()]
    hashes={p.relative_to(root).as_posix():sha(p) for p in sorted(files)};hashes[old_manifest.relative_to(root).as_posix()]=sha(old_manifest)
    report=dict(schema='soil-fem-richards-dam-final-validation/1',passed=True,solver_sha256=digest,source_snapshot_sha256=sha(run/'source-snapshot.zip'),rust_tests=core_report['rust_tests'],cad_tests=5,joint_split_tests=2,independent_reference_tests=2,separate_process_restarts=8,exact_legacy_cases=6,
        core_manifest_sha256=core['manifest_sha256'],comparisons=dam['comparisons'],independent_reference_refinements=dam['independent_reference_refinements'],native_refinements=dam['native_refinements'],exact_counterexample_fixed=True,
        scope='Fresh complete hydraulic-core suite plus two steady homogeneous cube benchmarks, including genuinely finite-width reservoir flow. Local residual certification unchanged; fixed-transmissivity L2 merit globalizes Newton. No suction mechanics, domain-change water transfer, pond dynamics or general geology certification.',sha256=hashes)
    with (run/'final-validation.json').open('x',encoding='utf-8') as f:json.dump(report,f,indent=2)
    print(json.dumps({k:v for k,v in report.items() if k not in ('sha256','comparisons','native_refinements','independent_reference_refinements')},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);record(p.parse_args().run.resolve())
