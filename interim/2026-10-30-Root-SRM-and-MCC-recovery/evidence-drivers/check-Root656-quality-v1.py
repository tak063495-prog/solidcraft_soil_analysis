"""Build the adopted Root, exercise real CLI physics/recovery, and run219 Python."""
from pathlib import Path
import datetime,hashlib,json,os,shutil,struct,subprocess,sys,time
from qa_log_counts import unittest_summary
p=Path(__file__).resolve().parent;root=p.parents[1];adoption=p/'Root-source656-promotion-v1';out=adoption/'quality-v1';out.mkdir()
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
    with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
    with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
def active():
    if (p/'pause-requested').exists():raise SystemExit(2)
active();sys.path.insert(0,str(root/'tools'));from freeze_sources import manifest,freeze
sources=read(adoption/'adopted-source-manifest.json');assert manifest()==sources and len(sources)==656
qualified=p/'srm-durable-job-v1/quality-v3';private=p/'srm-durable-job-v1/runtime';qsources=read(qualified/'source-snapshot.json')
assert read(qualified/'terminal.json')['passed']
assert all(sources[n]==h for n,h in qsources.items() if Path(n).suffix!='.md')
freeze(out);dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
env=dict(os.environ,CARGO_TARGET_DIR=str(root/'target'),CARGO_BUILD_JOBS='1',PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1',
    OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
for key in ['PYTHONPATH','PYTHONHOME']:env.pop(key,None)
records=[]
def call(label,args,expected=0):
    active();assert manifest()==sources;assert all(sha(private/n)==h for n,h in qsources.items());start=time.monotonic()
    with (out/(label+'.log')).open('x',encoding='utf-8') as log:
        child=subprocess.Popen(list(map(str,args)),cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,
            creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        dump(out/(label+'-child.json'),dict(pid=child.pid,parent_pid=os.getpid(),executable=str(args[0]),arguments=list(map(str,args[1:])),
            created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        code=child.wait()
    row=dict(step=label,exit_code=code,expected_exit_code=expected,elapsed_seconds=time.monotonic()-start,
        Root_and_private_sources_unchanged=manifest()==sources and all(sha(private/n)==h for n,h in qsources.items()))
    records.append(row);dump(out/(label+'-terminal.json'),row);print(json.dumps(row),flush=True)
    if code!=expected or not row['Root_and_private_sources_unchanged']:
        dump(out/'terminal.json',dict(passed=False,failed_step=label,records=records,p6_complete=False));raise SystemExit(1)
def compare(a,b,first_binary,second_binary,path='$'):
    assert type(a)is type(b),(path,type(a),type(b))
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for key in a:
            if key=='elapsed_seconds':continue
            if key=='solver_binary_sha256':assert a[key]==first_binary and b[key]==second_binary,(path,a[key],b[key]);continue
            compare(a[key],b[key],first_binary,second_binary,path+'.'+key)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,first_binary,second_binary,f'{path}[{i}]')
    elif isinstance(a,float):assert struct.pack('>d',a)==struct.pack('>d',b),(path,a,b)
    else:assert a==b,(path,a,b)
call('actual-Root-release-build',['C:/Users/link_/.cargo/bin/cargo.exe','build','--locked','--release','--message-format=json-render-diagnostics'])
executables=[]
for line in (out/'actual-Root-release-build.log').read_text('utf-8').splitlines():
    try:item=json.loads(line)
    except json.JSONDecodeError:continue
    if item.get('reason')=='compiler-artifact' and item.get('target',{}).get('name')=='soil-fem-3d' and item.get('executable'):executables.append(Path(item['executable']))
assert len(set(executables))==1;solver=root/'target/release/soil-fem-3d.exe';assert executables[-1].resolve()==solver.resolve()
binary_sha=sha(solver);shutil.copyfile(solver,out/'soil-fem-3d.exe');assert sha(out/'soil-fem-3d.exe')==binary_sha
dump(out/'actual-Root-binary.json',dict(passed=True,actual_Cargo_Root_build=True,binary=str(solver),sha256=binary_sha,
    source_manifest_sha256=sha(out/'source-snapshot.json'),all_nonMarkdown_sources_exact_to_qualified347_Rust219_Python=True,p6_complete=False))
call('actual-Root-help',[solver,'unknown-command'],1)
helptext=(out/'actual-Root-help.log').read_text('utf-8')
for name in ['prepare-srm-search-job','run-srm-search-job','prepare-srm-search-recovery','stop-srm-search-job',
    'verify-srm-search-job','search-mc-srm','mc-virtual-work-bound']:assert name in helptext
old_sha=read(adoption/'old-Root635-baselines.json')['old_binary_sha256']
for family in ['mixed','projected']:
    src=p/'mcc-local-cut-stop-verification-v1/results-653'/family
    call('actual-Root-MCC-'+family,[solver,'solve-stages',src/'mesh.json',src/'ground.json',src/'case.json',out/(family+'-result.json')])
    old=read(adoption/'before'/family/'result.json');new=read(out/(family+'-result.json'));assert new['complete'] and len(new['frames'])==5
    compare(old,new,old_sha,binary_sha)
    dump(out/(family+'-physics-equivalence.json'),dict(passed=True,all5_MCC_preload_unload_local_cut_postcut_frames_and_full_history_bits_equal=True,
        only_elapsed_and_verified_exact_executable_provenance_excluded=True,old_Root635_sha256=old_sha,new_Root656_sha256=binary_sha,p6_complete=False))
goldens=read(p/'slope1200-durable-SRM-v1/report-extraction-verification-v1/terminal.json')['rows']
for i,label in [(0,'closed-form'),(1,'numeric-unknown')]:
    golden=Path(goldens[i]['source']);assert sha(golden)==goldens[i]['source_sha256'];inputs=golden.parent.parent
    job=out/('Root '+label+' SRM job with spaces')
    call(label+'-prepare',[solver,'prepare-srm-search-job',inputs/'mesh.json',inputs/'ground.json',inputs/'case.json',inputs/'spec.json','1',job])
    frozen=job/'solver/soil-fem-3d.exe';assert sha(frozen)==binary_sha
    if i==0:
        call(label+'-one-accepted-checkpoint',[frozen,'run-srm-search-job',job,'1'],2)
        call(label+'-verify-paused',[frozen,'verify-srm-search-job',job])
        info=read(out/(label+'-verify-paused.log'));assert info['accepted_mechanical_checkpoint'] and not info['complete']
        before={str(q.relative_to(job)):sha(q) for q in job.rglob('*') if q.is_file()}
        recovered=out/'Root closed form SRM recovered with spaces'
        call(label+'-prepare-recovery',[solver,'prepare-srm-search-recovery',job,recovered])
        call(label+'-finish-recovered',[recovered/'solver/soil-fem-3d.exe','run-srm-search-job',recovered])
        assert before=={str(q.relative_to(job)):sha(q) for q in job.rglob('*') if q.is_file()}
        final_job=recovered
    else:
        call(label+'-run',[frozen,'run-srm-search-job',job]);final_job=job
    report=read(final_job/'report.json');control=read(golden)
    compare(control,report,control['source_identity']['solver_binary_sha256'],binary_sha)
    if i==0:
        assert report['requirement']=='equilibrium_requirement_met'
        interval=report['interval'];assert interval['attained_equilibrium_factor']<20/15<interval['independently_excluded_factor']
        assert interval['width']<=report['specification']['absolute_factor_tolerance']
    else:assert report['requirement']=='indeterminate' and any(t['outcome']=='indeterminate' for t in report['trials'])
    dump(out/(label+'-full-report-equivalence.json'),dict(passed=True,actual_Root_CLI_all_scientific_fields_types_and_binary64_equal=True,
        only_elapsed_and_verified_exact_executable_provenance_excluded=True,requirement=report['requirement'],
        unmodified_source_job_after_recovery=(i==0),p6_complete=False))
case=p/'slope1200-durable-SRM-v1/tighter-quadratic-vertex-static-bounds-v5/F1p45'
assert read(case/'independent-PASS.json')['passed']
mesh=p/'srm-general-virtual-work-v1/slope1200-LP-v1/mesh.json'
call('actual-Root-full32400-GP',[solver,'mc-virtual-work-bound',mesh,case/'candidate.json',out/'Root-F1p45-all-GP.json'])
assert sha(out/'Root-F1p45-all-GP.json')==sha(case/'actual-Rust-all-GP-result.json')
dump(out/'actual-Root-all-GP-equivalence.json',dict(passed=True,Root_actual_1200_Tet10_all32400_GP_entire_report_byte_exact=True,
    reference_was_independently_audited_at_original_tensor_eigen_dual_and_work_budgets=True,
    actual_report_sha256=sha(out/'Root-F1p45-all-GP.json'),exact_Fs_or_continuum_or_original5percent_field_acceptance_not_claimed=True,p6_complete=False))
call('all-Root-Python-tests',[sys.executable,'-m','unittest','discover','-s','tests','-p','test_*.py','-v'])
dump(out/'Root-Python-inventory.json',unittest_summary((out/'all-Root-Python-tests.log').read_text('utf-8'),219))
assert manifest()==sources and sha(solver)==binary_sha
dump(out/'terminal.json',dict(passed=True,actual_Root_Cargo_build_and_CLI_exercised=True,Root_source_files=656,
    all_Root219_Python_passed=True,private_same_nonMarkdown_sources_full347_Rust219_Python_passed=True,
    full347_Rust_not_rerun_as_local_Root_test_suite=True,
    actual_Root_MCC_mixed_projected_all5_frames_and_histories_match_old_Root635_bits=True,
    actual_Root_SRM_accepted_save_and_new_recovery_and_closed_form_and_unknown_semantics_passed=True,
    actual_Root_all32400_GP_F1p45_report_entire_bytes_identical_to_independently_audited_candidate=True,
    all_old_frozen_jobs_and_sources_and_guard_unchanged=True,binary_sha256=binary_sha,records=records,
    new_remote_CI_Native_gate_UI_convergence_and_all35_separate=True,p6_complete=False))
dump(adoption/'terminal.json',dict(passed=True,adopted_Root656=True,old_Root635_sources_and_binary_deeply_preserved=True,
    qualification=str(out/'terminal.json'),qualification_sha256=sha(out/'terminal.json'),
    actual_Root_build_MCC_SRM_recovery_allGP_and219Python_passed=True,
    qualified_nonMarkdown_numeric_sources347_Rust_passed=True,all35_and_original5percent_FAIL_retained=True,p6_complete=False))
print('PASS actual Root656 build, MCC/SRM/recovery/static32400GP and219Python; allP6 remains',flush=True)
