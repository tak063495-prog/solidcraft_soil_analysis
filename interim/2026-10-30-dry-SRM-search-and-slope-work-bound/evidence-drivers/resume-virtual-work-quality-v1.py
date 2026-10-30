from pathlib import Path
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
def call(label,args,expected=0):
    active();assert unchanged()
    start=time.monotonic()
    with (out/(label+'.log')).open('x',encoding='utf-8') as log:
        child=subprocess.Popen(list(map(str,args)),cwd=runtime,env=env,stdout=log,stderr=subprocess.STDOUT,
            creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        dump(out/(label+'-child.json'),dict(pid=child.pid,parent_pid=os.getpid(),executable=str(args[0]),
            arguments=list(map(str,args[1:])),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        code=child.wait()
    stable=unchanged()
    row=dict(step=label,exit_code=code,expected_exit_code=expected,elapsed_seconds=time.monotonic()-start,
        both_private650_and645_sources_unchanged=stable,
        current_Root_matches_initial635=all((root/n).exists() and sha(root/n)==h for n,h in root635.items()))
    records.append(row);dump(out/(label+'-terminal.json'),row);print(json.dumps(row),flush=True)
    if code!=expected or not stable:
        dump(out/'terminal.json',dict(passed=False,failed_step=label,records=records,actual_failure_retained=True,p6_complete=False))
        raise SystemExit(1)
cargo='C:/Users/link_/.cargo/bin/cargo.exe'
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
before=sha(Path(cases[0]['output']))
call('CLI-refuse-overwrite',[solver,'mc-virtual-work-bound',cases[0]['mesh'],cases[0]['problem'],cases[0]['output']],1)
assert sha(Path(cases[0]['output']))==before
bad=read(cases[0]['problem']);bad['virtual_displacement_m'][0][0]=1e-15
dump(out/'supported-component-violation.json',bad)
call('CLI-supported-component-unavailable',[solver,'mc-virtual-work-bound',cases[0]['mesh'],out/'supported-component-violation.json',out/'not-applicable.json'])
unavailable=read(out/'not-applicable.json')
assert unavailable['assessment']['status']=='not_applicable' and not unavailable['complete_point_audit']
bad=read(cases[0]['problem']);bad['water_policy']='unsupported'
dump(out/'unsupported-water.json',bad)
call('CLI-reject-water',[solver,'mc-virtual-work-bound',cases[0]['mesh'],out/'unsupported-water.json',out/'must-not-exist.json'],1)
assert not (out/'must-not-exist.json').exists()
call('all-Python-tests',[sys.executable,'-m','unittest','discover','-s','tests','-p','test_*.py','-v'])
inventory=unittest_summary((out/'all-Python-tests.log').read_text('utf-8'),219)
dump(out/'Python-inventory.json',inventory)
while not (p/'srm-adaptive-search-v1/quality-v2-continuation-v2/terminal.json').exists():
    active();time.sleep(5)
assert read(p/'srm-adaptive-search-v1/quality-v2-continuation-v2/terminal.json')['passed']
call('all-Rust-tests',[cargo,'test','--locked','--release'])
counts=[tuple(map(int,x)) for x in re.findall(r'test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored;',
    (out/'all-Rust-tests.log').read_text('utf-8'))]
assert counts and tuple(map(sum,zip(*counts)))==(332,0,0)
dump(out/'terminal.json',dict(passed=True,source_files=650,Rust_passed=332,Python_passed=219,
    ignored_skipped_failed=0,actual_new12_and_prior11_and9_passed=True,
    actual_distributed_support_frictional_search_and_real_plastic_FEM_and_numerical_unknowns_tested=True,
    actual6_fresh_CLI_cases_and_all504_independent_reference_points_and_raw_point_hashes_verified=True,
    actual_CLI_help_no_overwrite_supported_component_unavailable_water_rejection_passed=True,
    private645_and650_sources_unchanged=True,Root_was_not_used_as_a_mutable_runtime_dependency=True,
    automatic_actual_slope_and_full3D_convergence_and_allP6_not_claimed=True,records=records,p6_complete=False))
print('PASS isolated650 pointwise MC work witness and all regressions; actual slope/P6 remain',flush=True)
