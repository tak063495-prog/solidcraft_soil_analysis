"""Resume frozen650 QA after identifying independent quadrature family mismatch."""
from pathlib import Path
p=Path(__file__).resolve().parent
original=(p/'check-virtual-work-quality.py').read_text('utf-8')
header='''from pathlib import Path
import datetime,hashlib,json,os,re,shutil,subprocess,sys,time
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['MKL_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from qa_log_counts import unittest_summary
import virtual_work_references as reference
p=Path(__file__).resolve().parent;root=p.parents[1]
candidate=p/'srm-general-virtual-work-v1';runtime=candidate/'runtime'
old=candidate/'quality-v1';out=candidate/'quality-v1-continuation-v1';out.mkdir()
read=lambda q:json.loads(Path(q).read_text('utf-8-sig'))
def sha(q):
    with Path(q).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
    with Path(q).open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
def active():
    if (p/'pause-requested').exists():raise SystemExit(2)
active()
sys.path.insert(0,str(runtime/'tools'))
from freeze_sources import manifest,freeze
sources=read(old/'source-snapshot.json');base=read(candidate/'baseline645.json')
assert len(sources)==650 and len(base)==645
prior=p/'srm-adaptive-search-v1/runtime'
def unchanged():return manifest()==sources and all(sha(prior/n)==h for n,h in base.items())
assert unchanged();freeze(out);assert read(out/'source-snapshot.json')==sources
root635=read(p/'Root-source635-promotion-v1/quality-v1/source-snapshot.json')
solver=old/'soil-fem-3d.exe'
assert solver.exists() and sha(solver)==sha(candidate/'target/release/soil-fem-3d.exe')
assert sha(solver)==sha(runtime/'target/release/soil-fem-3d.exe')
dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,
    script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
steps=['format','clippy','new-twelve-tests','prior-search-eleven-tests','prior-mean-nine-tests','CLI-build','CLI-usage']
records=[read(old/(n+'-terminal.json')) for n in steps]
assert all(r['exit_code']==r['expected_exit_code'] and r['both_private650_and645_sources_unchanged'] for r in records)
for n,count in [('new-twelve-tests',12),('prior-search-eleven-tests',11),('prior-mean-nine-tests',9)]:
    assert f'{count} passed; 0 failed; 0 ignored' in (old/(n+'.log')).read_text('utf-8')
dump(out/'previous-harness-terminal.json',dict(passed=False,actual_tool_session=55398,actual_driver_exit_code=1,
    actual_Tet10_CLI_exit_code0=True,actual_new12_and_prior11_and9_passed=True,
    failed_comparison='Tet10 point0 weight:0.008770474929651058 versus0.014972747367083927',
    cause='Independent v1 Legendre-Duffy points differ from Root shifted Gauss-Jacobi points; both global integrals valid',
    original_failure_logs_and_reference_v1_retained=True,Rust_source_or_tolerance_changes=0,
    independent_reference_v2_all56_exact_moments_checked=True,p6_complete=False))
refdir=candidate/'independent-FEM-v2';assert read(refdir/'terminal.json')['passed']
cases=reference.prepare(runtime,out,refdir)
oldcases=read(old/'independent-reference-binding.json')['cases']
reuse=[]
for case,oldcase in zip(cases,oldcases):
    assert case['index']==oldcase['index']
    assert sha(case['mesh'])==sha(oldcase['mesh']) and sha(case['problem'])==sha(oldcase['problem'])
    if case['index']<4:
        receipt=read(old/('CLI-independent-'+str(case['index'])+'-terminal.json'))
        assert receipt['exit_code']==0 and receipt['both_private650_and645_sources_unchanged']
        case['output']=oldcase['output']
        reuse.append(dict(index=case['index'],actual_result=case['output'],sha256=sha(case['output']),
            inputs_byte_identical=True,actual_CLI_receipt=receipt))
    else:assert not Path(oldcase['output']).exists()
dump(out/'reused-actual-CLI-results.json',dict(passed=True,binary_sha256=sha(solver),
    no_passing_Rust_tests_or_build_or_CLI_repeated_for_reference_fix=True,
    all6_inputs_byte_identical_to_previous_protocol=True,all_original2e_minus12_limits_unchanged=True,
    quadrature_reference_correction=read(candidate/'reference-family-correction-v1/preparation.json'),rows=reuse))
env=dict(os.environ,CARGO_TARGET_DIR=str(candidate/'target'),CARGO_BUILD_JOBS='1',RUST_TEST_THREADS='2',
    PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
for name in ['PYTHONPATH','PYTHONHOME']:env.pop(name,None)
'''
call=original[original.index('def call(label,'):original.index("cargo='C:/Users/link_/.cargo/bin/cargo.exe'")]
compare='''cargo='C:/Users/link_/.cargo/bin/cargo.exe'
verified=[]
for case in cases:
    if case['index']>=4:
        call('CLI-independent-'+str(case['index']),[solver,'mc-virtual-work-bound',case['mesh'],case['problem'],case['output']])
    verified.append(reference.verify(case))
assert sum(r['all_GP_checked'] for r in verified)==504
dump(out/'full-independent-reference-comparison.json',dict(passed=True,all6_actual_CLI_cases=True,
    all504_GP_records_checked=True,all_framed_point_hashes_from_raw_CLI_lexemes_verified=True,
    all_original2e_minus12_numeric_comparison_limits_unchanged=True,
    v2_reference_SciPy_Jacobi_family_and56_exact_moments_verified=True,
    four_prior_actual_CLI_results_reused_byte_exactly=True,cases=verified,p6_complete=False))
'''
tail=original[original.index("before=sha(Path(cases[0]['output']))"):]
target=p/'resume-virtual-work-quality-v1.py'
with target.open('x',encoding='utf-8') as f:f.write(header+call+compare+tail)
print(target)
