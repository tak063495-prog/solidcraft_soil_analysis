"""Inspect, or after weekly remaining falls below20%, preserve and stop owned jobs at a safe milestone."""
from pathlib import Path
import argparse, ctypes, ctypes.wintypes as w, datetime, hashlib, json, subprocess, sys, time

p = Path(__file__).resolve().parent
root = p.parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--execute', action='store_true')
parser.add_argument('--usage-percent', type=int)
args = parser.parse_args()
if args.execute and (args.usage_percent is None or args.usage_percent <= 80):
    parser.error('execution requires a fresh weekly usage observation above80% (remaining below20%)')
folders = [root/'runs'/n for n in [
    'native-mcc-finite-pit-reference-resume-20261028-v1',
    'native-mcc-compact-pit-spatial-fourth-20261028-v1',
    'native-mcc-compact-integration-20261028-v5',
    'mcc-acoustic-tool-prototype-20261028-v1',
    'native-mcc-finite-pit-temporal-20261028-v1',
    'clean-runtime-20261028-v2',
    'mcc-nonlinear-performance-20261028-v1',
    'native-mcc-finite-pit-regional-overlap-20261028-v1',
    'native-mcc-compact-pit-local-refinement-20261028-v1',
    'native-mcc-compact-pit-local-refinement-20261028-v2',
    'native-mcc-compact-pit-local-refinement-20261028-v3',
    'native-mcc-compact-pit-local-refinement-20261028-v4',
    'native-mcc-compact-pit-fixed-surface-refinement-20261028-v1',
    'native-mcc-compact-pit-fixed-surface-refinement-20261028-v2',
    'mesh-diagnostics-figure-20261028-v1',
    'native-mcc-compact-pit-1324-temporal-20261028-v1',
    'mcc-midpoint-quality-20261028-v10',
    'native-mcc-finite-pit-reference-rust-pair-20261028-v1',
    'mcc-midpoint-quality-20261028-v11',
    'native-mcc-fixed-wall-comparison-20261028-v1',
    'native-mcc-finite-pit-increment-residual-sensitivity-20261028-v1',
    'native-mcc-boundary-depth-20261028-v1',
    'mcc-large-structures-development-20261028-v1',
    'completed-studies-portable-solve-20261028-v1',
    'mcc-midpoint-quality-20261029-v13',
    'native-large-source-stream-20261029-v1',
    'published-large-case-reproduction-20261029-v1',
    'published-large-case-durable-reproduction-20261029-v2',
    'native-fem-smoke-20261029-v1',
    'bridge-install-20261029-v1',
    'native64-fresh-install-20261029-v1',
    'native65-fresh-install-20261029-v1',
    'native65-periodic-gui-recovery-20261029-v1',
    'fem-ci-262c8ae-success-20261029-v1',
    'fem-ci-3870c2f-success-20261029-v1',
    'native-coordinate-index-20261029-v1',
    'durable-stage-jobs-development-20261029-v1',
    'durable-root-promotion-20261029-v1',
    'native-durable-worker-development-20261029-v1',
    'native-large-worker-combined-20261029-v1',
    'native-large-worker-combined-20261029-v2',
    'native-durable-owned-ui-20261029-v1',
    'native-recovery-worker-development-20261029-v1',
    'native-snapshot-recovery-development-20261029-v1',
    'mcc-kernel-profiling-development-20261029-v1',
    'mcc-mixed-kernel-profiling-development-20261029-v1',
    'mcc-mixed-kernel-profiling-development-20261029-v2',
    'mcc-mixed-assembly-development-20261029-v1',
    'mc-stage-gate-development-20261029-v1',
    'mcc-parallel-material-development-20261029-v1',
    'mcc-fixed-wall-durable-pair-20261029-v1',
    'mcc-next-fixed-wall-durable-pair-20261029-v1',
    'fixed-wall-public-review-20261029-v1',
    'mcc-hydrostatic-scaling-20261029-v1',
    'mcc-hydrostatic-durable-restart-20261029-v1',
    'mcc-hydrostatic-100k-durable-20261029-v1',
    'mcc-parallel-100k-durable-20261029-v1',
    'native-mcc-compact-pit-fixed-surface-refinement-20261029-v3',
    'native-mcc-compact-pit-fixed-surface-refinement-20261029-v4',
    'richards-local-rain-3d-20261029-v1',
    'mc-stage-gate-root-promotion-20261029-v1',
    'mcc-105k-public-review-20261029-v1',
    'interim-upload-20261029-v32',
    'fem-ci-94e70ec-20261029-v1',
    'richards-local-rain-3d-refinement-20261029-v1',
    'native-boolean-kernel-regression-20261029-v1',
    'interim-upload-20261029-v33',
    'mcc-assembly-and-endpoints-combined-20261029-v1',
    'mc-parameter-transition-development-20261029-v1',
    'interim-upload-20261029-v34',
    'mcc-csr-root-promotion-20261029-v1',
    'mcc-csr-root-promotion-20261029-v2',
    'interim-upload-20261029-v35',
    'fem-ci-e00378f-20261029-v1',
    'mc-epoch-checkpoint-development-20261029-v1',
    'interim-upload-20261029-v36',
    'mc-epoch-cli-development-20261029-v1',
    'mc-epoch-root-promotion-20261029-v1',
    'interim-upload-20261029-v37',
    'fem-ci-a18d81b-20261029-v1',
    'native-kernel-cap-repair-20261029-v1',
    'mc-epoch-export-development-20261029-v1',
    'rain-refined-public-review-20261029-v1',
    'interim-upload-20261029-v38',
    'fem-ci-052a208-20261029-v1',
    'native-mc-epoch-view-development-20261029-v1',
    'interim-upload-20261029-v39',
    'native-mcc-compact-pit-fixed-surface-refinement-20261029-v5',
    'mcc-fixed-wall-finer14909-20261029-v1',
    'interim-upload-20261029-v40',
    'native-typed-source-development-20261029-v1',
    'native-typed-source-integration-20261029-v1',
    'interim-upload-20261029-v41',
    'native-source-point-view-development-20261030-v1',
]]
norm = lambda s: str(s).lower().replace('\\', '/')
prefixes = [norm(x)+'/' for x in folders]
ps = "Get-CimInstance Win32_Process | Select-Object ProcessId,ParentProcessId,ExecutablePath,CommandLine,@{n='Created';e={$_.CreationDate.ToUniversalTime().ToString('o')}} | ConvertTo-Json -Depth 2 -Compress"
records = json.loads(subprocess.check_output(['pwsh','-NoProfile','-Command',ps], text=True, encoding='utf-8-sig'))
records = records if isinstance(records, list) else [records]
records = {r['ProcessId']: r for r in records}
excluded={__import__('os').getpid()}
ancestor=records.get(next(iter(excluded)))
while ancestor and ancestor['ParentProcessId'] not in excluded:
    excluded.add(ancestor['ParentProcessId']);ancestor=records.get(ancestor['ParentProcessId'])
def scoped(r):
    command=norm(r.get('CommandLine') or '')
    return any(prefix in command or 'runs/'+prefix.split('/runs/',1)[1] in command for prefix in prefixes)
owned=set();seed_evidence=[]
# A path mentioned by an unrelated shell is not ownership evidence. Seed only
# processes recorded by our drivers, with matching image and creation time.
for folder in folders:
    receipts=set(folder.rglob('*-child.json'))|set(folder.rglob('child.json'))
    for path in receipts:
        try:
            evidence=json.loads(path.read_text(encoding='utf-8-sig'));pid=evidence['pid'];r=records[pid]
            if pid in excluded:continue
            observed=datetime.datetime.fromisoformat(r['Created']).timestamp()
            expected=datetime.datetime.fromisoformat(evidence['created_utc']).timestamp()
            if abs(observed-expected)>5:continue
            executable=evidence.get('executable')
            if executable:
                actual=norm(r.get('ExecutablePath') or '')
                expected_image=norm(executable)
                if Path(executable).is_absolute():
                    if actual!=expected_image:continue
                elif Path(actual).stem!=Path(expected_image).stem:continue
            elif not (evidence.get('script') and norm(evidence['script']) in norm(r.get('CommandLine') or '')):continue
            owned.add(pid);seed_evidence.append(dict(pid=pid,receipt=str(path)))
            child=r
            while child['ParentProcessId'] in records:
                parent=records[child['ParentProcessId']]
                if parent['ProcessId'] in excluded or not scoped(parent):break
                if datetime.datetime.fromisoformat(parent['Created'])>datetime.datetime.fromisoformat(child['Created']):break
                owned.add(parent['ProcessId']);child=parent
        except (OSError,ValueError,KeyError,TypeError):continue
# Only descendants of an explicitly scoped process are added. Never traverse upward.
while True:
    expanded = (owned | {pid for pid,r in records.items() if r['ParentProcessId'] in owned})-excluded
    if expanded == owned:
        break
    owned = expanded
plan = [records[pid] for pid in sorted(owned)]
if not args.execute:
    print(json.dumps(dict(inspect_only=True, usage_threshold_percent=80, threshold_comparison='strictly_greater_than', pause_at_next_safe_milestone=True, ownership_seed_receipts=seed_evidence, own_process_and_ancestors_excluded=True, owned_processes=plan), indent=2))
    raise SystemExit(0)

out = p/('pause-remaining20-'+datetime.datetime.now().strftime('%Y%m%d-%H%M%S'))
out.mkdir()
kernel = ctypes.WinDLL('kernel32', use_last_error=True)
kernel.OpenProcess.argtypes=[w.DWORD,w.BOOL,w.DWORD];kernel.OpenProcess.restype=w.HANDLE
kernel.CloseHandle.argtypes=[w.HANDLE]
kernel.WaitForSingleObject.argtypes=[w.HANDLE,w.DWORD]
kernel.QueryFullProcessImageNameW.argtypes=[w.HANDLE,w.DWORD,w.LPWSTR,ctypes.POINTER(w.DWORD)]
kernel.GetProcessTimes.argtypes=[w.HANDLE,*([ctypes.POINTER(w.FILETIME)]*4)]
kernel.TerminateProcess.argtypes=[w.HANDLE,w.UINT]
handles=[];saved={};stopped=[];cooperative={}
sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
try:
    # Hold all identities before snapshotting or requesting cooperative cancellation.
    for r in plan:
        h=kernel.OpenProcess(0x1000|0x100000|1,False,r['ProcessId'])
        if not h:
            continue
        handles.append((h,r))
        image=ctypes.create_unicode_buffer(32768);size=w.DWORD(len(image))
        assert kernel.QueryFullProcessImageNameW(h,0,image,ctypes.byref(size))
        assert norm(image.value)==norm(r['ExecutablePath'])
        times=[w.FILETIME() for _ in range(4)]
        assert kernel.GetProcessTimes(h,*[ctypes.byref(x) for x in times])
        actual=(times[0].dwHighDateTime<<32)|times[0].dwLowDateTime
        expected=(datetime.datetime.fromisoformat(r['Created']).timestamp()+11644473600)*10000000
        assert abs(actual-expected)<10000,'PID was reused'
    fourth=folders[1];temporal=folders[4];clean=folders[5]/'Compact MCC portable workspace'
    local=root/'runs/native-mcc-compact-pit-local-refinement-20261028-v1'
    local2=root/'runs/native-mcc-compact-pit-local-refinement-20261028-v2'
    local3=root/'runs/native-mcc-compact-pit-local-refinement-20261028-v3'
    local4=root/'runs/native-mcc-compact-pit-local-refinement-20261028-v4'
    fixed=root/'runs/native-mcc-compact-pit-fixed-surface-refinement-20261028-v2'
    temporal1324=root/'runs/native-mcc-compact-pit-1324-temporal-20261028-v1'
    reference_pair=root/'runs/native-mcc-finite-pit-reference-rust-pair-20261028-v1'
    fixed_comparison=root/'runs/native-mcc-fixed-wall-comparison-20261028-v1'
    increment_residual=root/'runs/native-mcc-finite-pit-increment-residual-sensitivity-20261028-v1'
    boundary_depth=root/'runs/native-mcc-boundary-depth-20261028-v1'
    large_structure=root/'runs/mcc-large-structures-development-20261028-v1'
    for driver in [folders[0],fourth,temporal,local,local2,local3,local4,fixed,temporal1324,reference_pair,fixed_comparison,increment_residual,boundary_depth,large_structure]:
        marker=driver/'pause-requested'
        if not marker.exists():marker.touch(exist_ok=False)
    jobs=[fourth/'workflow/analysis/solver-job',local/'workflow/analysis/solver-job',local2/'workflow/analysis/solver-job',local3/'workflow/analysis/solver-job',local4/'workflow/analysis/solver-job',fixed/'workflow/analysis/solver-job',*[temporal/(str(n)+'-v2')/'analysis/solver-job' for n in [32,64]]]
    jobs.append(temporal1324/'32-v2/analysis/solver-job')
    jobs.extend(reference_pair/str(n)/'analysis/solver-job' for n in [16,32])
    jobs.extend(increment_residual/n/'analysis/solver-job' for n in ['default64','tight16','tight32'])
    jobs.extend(boundary_depth/('depth'+str(d))/('mesh'+str(s))/'workflow/analysis/solver-job' for d in [6,8] for s in [4,2])
    jobs.extend(large_structure/'controls'/n/'analysis/solver-job' for n in ['default16','bounded16','bounded32'])
    jobs.extend(large_structure/'large-case'/n for n in ['solver-job-continuous','checkpoint-job','resume-job'])
    jobs.append(root/'runs/completed-studies-portable-solve-20261028-v1/Portable verification with spaces/solver-job')
    jobs.append(root/'runs/mcc-midpoint-quality-20261029-v13/default-control/analysis/solver-job')
    jobs.append(root/'runs/published-large-case-reproduction-20261029-v1/Portable full recomputation/workflow/analysis/solver-job')
    jobs.append(root/'runs/native-fem-smoke-20261029-v1/Native63 cooperative smoke/workflow/analysis/solver-job')
    jobs.extend(root/'runs/mcc-hydrostatic-scaling-20261029-v1'/('n'+str(n))/'solver-job' for n in [2,6,12,26])
    (root/'runs/mcc-hydrostatic-scaling-20261029-v1/pause-requested').touch(exist_ok=True)
    jobs.append(root/'runs/mcc-hydrostatic-durable-restart-20261029-v1/n12/solver-job')
    jobs.append(root/'runs/mcc-hydrostatic-durable-restart-20261029-v1/recovery-fork/n12/solver-job')
    jobs.append(root/'runs/mcc-hydrostatic-durable-restart-20261029-v1/loaded-recovery/n12/solver-job')
    (root/'runs/mcc-hydrostatic-durable-restart-20261029-v1/pause-requested').touch(exist_ok=True)
    jobs.append(root/'runs/mcc-hydrostatic-100k-durable-20261029-v1/n26/solver-job')
    jobs.append(root/'runs/mcc-kernel-profiling-development-20261029-v1/n6/solver-job')
    (root/'runs/mcc-kernel-profiling-development-20261029-v1/pause-requested').touch(exist_ok=True)
    jobs.append(root/'runs/mcc-mixed-kernel-profiling-development-20261029-v1/workflow/analysis/solver-job')
    (root/'runs/mcc-mixed-kernel-profiling-development-20261029-v1/pause-requested').touch(exist_ok=True)
    jobs.append(root/'runs/mcc-mixed-kernel-profiling-development-20261029-v2/workflow/analysis/solver-job')
    (root/'runs/mcc-mixed-kernel-profiling-development-20261029-v2/pause-requested').touch(exist_ok=True)
    (root/'runs/mcc-parallel-material-development-20261029-v1/pause-requested').touch(exist_ok=True)
    jobs.extend(root/'runs/mcc-parallel-material-development-20261029-v1'/('n'+str(n))/'solver-job' for n in [6,12])
    jobs.append(root/'runs/mcc-parallel-material-development-20261029-v1/controls/native274/workflow/analysis/solver-job')
    jobs.append(root/'runs/mcc-parallel-material-development-20261029-v1/loaded-recovery/n12/solver-job')
    jobs.append(root/'runs/mcc-parallel-material-development-20261029-v1/root-control/n6/solver-job')
    (root/'runs/mcc-hydrostatic-100k-durable-20261029-v1/pause-requested').touch(exist_ok=True)
    jobs.append(root/'runs/native-mcc-compact-pit-fixed-surface-refinement-20261029-v3/workflow/analysis/solver-job')
    (root/'runs/native-mcc-compact-pit-fixed-surface-refinement-20261029-v3/pause-requested').touch(exist_ok=True)
    (root/'runs/native-coordinate-index-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/durable-stage-jobs-development-20261029-v1/pause-requested').touch(exist_ok=True)
    jobs.extend(root/'runs/durable-stage-jobs-development-20261029-v1/real274-control'/n for n in ['legacy','durable'])
    jobs.extend(root/'runs/durable-root-promotion-20261029-v1/control'/n for n in ['legacy','durable'])
    jobs.extend(root/'runs/native-durable-worker-development-20261029-v1/controls'/n/'workflow/analysis/solver-job' for n in ['legacy','durable'])
    (root/'runs/native-durable-worker-development-20261029-v1/pause-requested').touch(exist_ok=True)
    jobs.append(root/'runs/native-large-worker-combined-20261029-v1/controls/durable/workflow/analysis/solver-job')
    (root/'runs/native-large-worker-combined-20261029-v1/pause-requested').touch(exist_ok=True)
    jobs.append(root/'runs/native-large-worker-combined-20261029-v2/controls/durable/workflow/analysis/solver-job')
    (root/'runs/native-large-worker-combined-20261029-v2/pause-requested').touch(exist_ok=True)
    jobs.append(root/'runs/native-durable-owned-ui-20261029-v1/desktop/owned-job/job/workflow/analysis/solver-job')
    jobs.extend(root/'runs/native-recovery-worker-development-20261029-v1/controls'/n/'workflow/analysis/solver-job' for n in ['continuous','abrupt-loss','recovered','root-mid-excavation'])
    (root/'runs/native-recovery-worker-development-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/native-snapshot-recovery-development-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/native-durable-owned-ui-20261029-v1/pause-requested').touch(exist_ok=True)
    jobs.extend(root/'runs/mcc-fixed-wall-durable-pair-20261029-v1'/('mesh'+str(n))/'workflow/analysis/solver-job' for n in [2865,5144])
    (root/'runs/mcc-fixed-wall-durable-pair-20261029-v1/pause-requested').touch(exist_ok=True)
    jobs.append(root/'runs/mcc-mixed-assembly-development-20261029-v1/workflow/analysis/solver-job')
    jobs.append(root/'runs/mcc-mixed-assembly-development-20261029-v1/large2865/workflow/analysis/solver-job')
    jobs.append(root/'runs/mcc-mixed-assembly-development-20261029-v1/recovery/workflow/analysis/solver-job')
    jobs.append(root/'runs/mcc-assembly-and-endpoints-combined-20261029-v1/workflow/analysis/solver-job')
    jobs.append(root/'runs/mcc-assembly-and-endpoints-combined-20261029-v1/recovery/workflow/analysis/solver-job')
    jobs.append(root/'runs/mcc-parallel-100k-durable-20261029-v1/n26/solver-job')
    jobs.append(root/'runs/mcc-parallel-100k-durable-20261029-v1/loaded-recovery/n26/solver-job')
    jobs.append(root/'runs/native65-periodic-gui-recovery-20261029-v1/desktop/recovered-owned/job/workflow/analysis/solver-job')
    (root/'runs/fem-ci-3870c2f-success-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/fem-ci-e00378f-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/fem-ci-a18d81b-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/native-kernel-cap-repair-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/mc-epoch-export-development-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/rain-refined-public-review-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/interim-upload-20261029-v38/pause-requested').touch(exist_ok=True)
    (root/'runs/fem-ci-052a208-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/native-mc-epoch-view-development-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/native-mcc-compact-pit-fixed-surface-refinement-20261029-v5/pause-requested').touch(exist_ok=True)
    (root/'runs/mcc-fixed-wall-finer14909-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/native-typed-source-development-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/native-typed-source-integration-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/native-source-point-view-development-20261030-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/mc-epoch-checkpoint-development-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/mc-epoch-cli-development-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/native65-periodic-gui-recovery-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/mcc-mixed-assembly-development-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/mcc-parallel-100k-durable-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/fixed-wall-public-review-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/native-mcc-compact-pit-fixed-surface-refinement-20261029-v4/pause-requested').touch(exist_ok=True)
    jobs.extend(root/'runs/mcc-next-fixed-wall-durable-pair-20261029-v1'/('mesh'+str(n))/'workflow/analysis/solver-job' for n in [5144,8061])
    jobs.append(root/'runs/mcc-fixed-wall-finer14909-20261029-v1/mesh14909/workflow/analysis/solver-job')
    jobs.append(root/'runs/native-source-point-view-development-20261030-v1/mesh14909/workflow/analysis/solver-job')
    jobs.append(root/'runs/native-source-point-view-development-20261030-v1/global-volume-sensitivity-study-v1/mesh16088/workflow/analysis/solver-job')
    jobs.append(root/'runs/native-source-point-view-development-20261030-v1/native-owned-v4-real-job-v1/run/owned-case/job/workflow/analysis/solver-job')
    # Registered original1200 origin: same immutable new executable and original inputs.
    registered_origin = root/'runs/native-source-point-view-development-20261030-v1/slope1200-epoch-SRM-origin-v1/origin-job'
    if (registered_origin/'job.json').exists():
        assert sha(registered_origin/'job.json') == '316ef07bfd6708cc1b112d757c7b1fef981b2fc183bdd886574aaa0750847ad6', 'registered epoch origin manifest changed'
        registered_manifest = json.loads((registered_origin/'job.json').read_text(encoding='utf-8-sig'))
        assert registered_manifest['schema'] == 'soil-fem-stages-job/1'
        assert registered_manifest['job_id'] == '27d31292fdac2853dc5589716bcd5530bf7ff8bba529cb07a4b3cc70cc68568a'
        assert registered_manifest['binding']['checkpoint_interval'] == 1
        assert sha(registered_origin/'solver/soil-fem-3d.exe') == 'f205ef2345354fb889d114c2c2508574819f8028fea7ceec3fc8b37b0d865ee7'
        assert registered_manifest['binding']['files_sha256']['solver/soil-fem-3d.exe'] == 'f205ef2345354fb889d114c2c2508574819f8028fea7ceec3fc8b37b0d865ee7'
        jobs.append(registered_origin)
    jobs.append(root/'runs/native-source-point-view-development-20261030-v1/plane-strain-slope-comparison-v1/xy1-ny1/solver-job')
    jobs.append(root/'runs/native-source-point-view-development-20261030-v1/plane-strain-slope-comparison-v1/xy1-ny2/solver-job')
    jobs.append(root/'runs/native-source-point-view-development-20261030-v1/plane-strain-slope-comparison-v1/xy2-ny2/solver-job')
    jobs.append(root/'runs/native-source-point-view-development-20261030-v1/plane-strain-slope-comparison-v1/xy1-ny1/recovery-probe-v1')
    jobs.append(root/'runs/native-source-point-view-development-20261030-v1/plane-strain-slope-comparison-v1/xy1-ny1-global64/solver-job')
    jobs.append(root/'runs/native-source-point-view-development-20261030-v1/plane-strain-slope-comparison-v1/xy1-ny1-local-half/solver-job')
    jobs.append(root/'runs/native-source-point-view-development-20261030-v1/Root635-default-owned-workflow-v1/run/owned-case/job/workflow/analysis/solver-job')
    jobs.append(root/'runs/native-source-point-view-development-20261030-v1/fixed-wall-finer-study-v1/mesh34350/workflow/analysis/solver-job')
    (root/'runs/mcc-next-fixed-wall-durable-pair-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/mc-stage-gate-development-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/richards-local-rain-3d-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/mc-stage-gate-root-promotion-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/mcc-105k-public-review-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/interim-upload-20261029-v32/pause-requested').touch(exist_ok=True)
    (root/'runs/fem-ci-94e70ec-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/richards-local-rain-3d-refinement-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/native-boolean-kernel-regression-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/interim-upload-20261029-v33/pause-requested').touch(exist_ok=True)
    (root/'runs/mcc-assembly-and-endpoints-combined-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/mc-parameter-transition-development-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/interim-upload-20261029-v34/pause-requested').touch(exist_ok=True)
    jobs.append(root/'runs/mc-stage-gate-root-promotion-20261029-v1/control/workflow/analysis/solver-job')
    portable_marker=root/'runs/completed-studies-portable-solve-20261028-v1/pause-requested'
    portable_marker.touch(exist_ok=True)
    (root/'runs/mcc-midpoint-quality-20261029-v13/pause-requested').touch(exist_ok=True)
    (root/'runs/native-large-source-stream-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/published-large-case-reproduction-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/mcc-csr-root-promotion-20261029-v1/pause-requested').touch(exist_ok=True)
    jobs.append(root/'runs/mcc-csr-root-promotion-20261029-v1/control/workflow/analysis/solver-job')
    (root/'runs/mcc-csr-root-promotion-20261029-v2/pause-requested').touch(exist_ok=True)
    jobs.append(root/'runs/mcc-csr-root-promotion-20261029-v2/control/workflow/analysis/solver-job')
    (root/'runs/mc-epoch-root-promotion-20261029-v1/pause-requested').touch(exist_ok=True)
    jobs.append(root/'runs/mc-epoch-root-promotion-20261029-v1/control/workflow/analysis/solver-job')
    jobs.append(root/'runs/published-large-case-durable-reproduction-20261029-v2/workflow/analysis/solver-job')
    (root/'runs/published-large-case-durable-reproduction-20261029-v2/pause-requested').touch(exist_ok=True)
    (root/'runs/native-fem-smoke-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/bridge-install-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/native64-fresh-install-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/native65-fresh-install-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/fem-ci-262c8ae-success-20261029-v1/pause-requested').touch(exist_ok=True)
    (root/'runs/published-large-case-reproduction-20261029-v1/published-bundle/pause-requested').touch(exist_ok=True)
    gui_folders=[root/'runs/native-large-source-stream-20261029-v1'/n for n in ['desktop','desktop-ci']]
    gui_folders.extend(root/'runs/native-coordinate-index-20261029-v1'/n for n in ['desktop','desktop-release'])
    gui_folders.append(root/'runs/native-durable-owned-ui-20261029-v1/desktop')
    gui_folders.append(root/'runs/native-snapshot-recovery-development-20261029-v1/desktop-release')
    gui_folders.append(root/'runs/native65-periodic-gui-recovery-20261029-v1/desktop')
    for gui in gui_folders:
        if (gui/'launch.json').exists() and not (gui/'owned-app-exited.json').exists():
            launch=json.loads((gui/'launch.json').read_text(encoding='utf-8'))
            if launch['process_id'] in owned:
                tag='safe-quit-'+datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
                result=subprocess.run([sys.executable,str(gui/'gui-control.py'),tag,str(gui/'request-safe-quit.json')],capture_output=True,text=True,timeout=60)
                assert result.returncode==0 and (gui/'owned-app-exited.json').exists(),'owned app graceful quit/preservation still pending'
    native_job=clean/'job';native_status=native_job/'status.json'
    if native_status.exists() and json.loads(native_status.read_text(encoding='utf-8-sig'))['state'] not in ['succeeded','failed','cancelled','checkpointed']:
        request=subprocess.run([str(clean/'.venv/Scripts/python.exe'),str(clean/'tools/native_job.py'),'cancel',str(native_job)],capture_output=True,text=True,timeout=30)
        assert request.returncode==0
        for _ in range(36):
            time.sleep(5)
            if json.loads(native_status.read_text(encoding='utf-8-sig'))['state'] in ['succeeded','failed','cancelled','checkpointed']:
                break
        assert json.loads(native_status.read_text(encoding='utf-8-sig'))['state'] in ['succeeded','failed','cancelled','checkpointed'],'native cooperative save still pending'
    jobs.append(native_job/'workflow/analysis/solver-job')
    from match_stage_job_process import matches_stage_job_process
    def windows_argv(command):
        count=ctypes.c_int()
        parse=ctypes.WinDLL('shell32',use_last_error=True).CommandLineToArgvW
        parse.argtypes=[w.LPCWSTR,ctypes.POINTER(ctypes.c_int)];parse.restype=ctypes.POINTER(w.LPWSTR)
        pointer=parse(command,ctypes.byref(count))
        if not pointer:raise OSError('Cannot read owned MC process arguments')
        try:return [pointer[i] for i in range(count.value)]
        finally:
            free=ctypes.WinDLL('kernel32',use_last_error=True).LocalFree
            free.argtypes=[ctypes.c_void_p];free.restype=ctypes.c_void_p
            free(ctypes.cast(pointer,ctypes.c_void_p))
    for job in jobs:
        status=job/'status.json'
        if not status.exists():continue
        value=json.loads(status.read_text(encoding='utf-8-sig'))
        if value['state']!='running':continue
        if value.get('process_id') not in owned:continue
        if not matches_stage_job_process(job,value,owned,records,windows_argv):continue
        solver=job/'solver/soil-fem-3d.exe'
        manifest=json.loads((job/'job.json').read_text(encoding='utf-8-sig'))
        assert solver.is_file() and sha(solver)==manifest['binding']['files_sha256']['solver/soil-fem-3d.exe'],'frozen job solver identity differs'
        result=subprocess.run([str(solver),'stop-stages-job',str(job)],capture_output=True,text=True,timeout=30)
        cooperative[str(job)]=dict(exit_code=result.returncode,stdout=result.stdout,stderr=result.stderr)
        assert result.returncode==0
    for job in jobs:
        status=job/'status.json'
        if not status.exists():continue
        if str(job) not in cooperative:continue
        for _ in range(36):
            value=json.loads(status.read_text(encoding='utf-8-sig'))
            if value['state']!='running':break
            time.sleep(5)
        if value['state']=='running':
            # A large Newton trial may outlast the cooperative wait. Discard
            # only its unaccepted trial after a copied accepted snapshot has
            # passed full Rust history verification in a new immutable job.
            from preserve_durable_snapshot import preserve_snapshot
            key=hashlib.sha256(str(job).encode('utf-8')).hexdigest()[:16]
            attempts=[];preserved=None
            for attempt in range(3):
                try:
                    frozen=job/'solver/soil-fem-3d.exe'
                    binding=json.loads((job/'job.json').read_text(encoding='utf-8-sig'))['binding']
                    expected=binding['files_sha256']['solver/soil-fem-3d.exe']
                    preserved=preserve_snapshot(job,out/'accepted-recovery-forks'/key/str(attempt),frozen,expected)
                    break
                except (OSError,ValueError,subprocess.TimeoutExpired) as error:
                    attempts.append(dict(attempt=attempt,error=str(error)))
                    if not any(word in str(error) for word in ['changed','advanced']):break
            cooperative[str(job)]['snapshot_preservation_attempt_errors']=attempts
            assert preserved and preserved['passed'],'no verified accepted snapshot; do not kill the pending solver'
            cooperative[str(job)]['accepted_snapshot_preserved_before_owned_stop']=preserved
            cooperative[str(job)]['discarded_unaccepted_trial_only']=True
        cooperative.setdefault(str(job),{})['final_state']=value['state']
        cooperative[str(job)]['final_state_is_last_report_not_process_liveness']=True
    # MC parameter jobs use their own sealed origin/material/path schema.
    # Select only currently owned process IDs established above by held handles
    # and creation/image identity; a stale running status never grants ownership.
    mc_scopes=[root/'runs/mc-epoch-cli-development-20261029-v1',root/'runs/mc-epoch-root-promotion-20261029-v1',root/'runs/native-mc-epoch-view-development-20261029-v1']
    mc_jobs=[]
    for manifest_path in (q for scope in mc_scopes for q in scope.rglob('job.json')):
        try:
            manifest=json.loads(manifest_path.read_text(encoding='utf-8-sig'))
            if manifest.get('schema')!='soil-fem-mc-epoch-job/1':continue
            job=manifest_path.parent
            value=json.loads((job/'status.json').read_text(encoding='utf-8-sig'))
            if value.get('state')!='running' or value.get('process_id') not in owned:continue
            record=records[value['process_id']]
            command=windows_argv(record.get('CommandLine') or '')
            if len(command)<3 or command[1]!='run-mc-epoch-job' or not Path(command[2]).is_absolute():continue
            if norm(Path(command[2]).resolve())!=norm(job.resolve()):continue
            if not record.get('ExecutablePath') or sha(Path(record['ExecutablePath']))!=manifest['binding']['files_sha256']['solver/soil-fem-3d.exe']:continue
            mc_jobs.append((job,manifest))
        except (OSError,ValueError):continue
    for job,manifest in mc_jobs:
        solver=job/'solver/soil-fem-3d.exe'
        expected=manifest['binding']['files_sha256']['solver/soil-fem-3d.exe']
        assert sha(solver)==expected,'MC frozen executable differs'
        process=subprocess.run([str(solver),'stop-mc-epoch-job',str(job)],capture_output=True,text=True,timeout=30)
        assert process.returncode==0,'MC cooperative stop request failed'
        for _ in range(36):
            value=json.loads((job/'status.json').read_text(encoding='utf-8-sig'))
            if value['state']!='running':break
            time.sleep(5)
        from preserve_mc_epoch_job import preserve_epoch_job
        key=hashlib.sha256(str(job).encode('utf-8')).hexdigest()[:16]
        fork=preserve_epoch_job(job,out/'mc-accepted-recovery-forks'/key,solver,expected)
        assert fork['passed'],'MC last accepted state was not deeply preserved; do not stop owned processes'
        cooperative[str(job)]=dict(material_epoch_job=True,stop_request_exit_code=process.returncode,
            final_state=value['state'],final_state_is_last_report_not_process_liveness=True,
            accepted_snapshot_preserved_before_owned_stop=fork,discarded_unaccepted_trial_only=True)
    # registered_slope1200_srm_job: exact path/manifest/executable only.
    # Ownership still requires held PID creation/image identity and exact argv.
    srm_job = root / 'runs/native-source-point-view-development-20261030-v1/slope1200-durable-SRM-v1/job'
    srm_manifest_path = srm_job / 'job.json'
    if srm_manifest_path.exists():
        srm_manifest = json.loads(srm_manifest_path.read_text(encoding='utf-8-sig'))
        assert sha(srm_manifest_path) == 'c5ef8c9c91d4ae2e3259f713c65e2173faa744073885b10024771182abac1bd8', 'registered SRM manifest changed'
        assert srm_manifest.get('schema') == 'soil-fem-srm-search-job/1'
        assert srm_manifest['job_id'] == '2075eb773e18cfabb79a2a5477828c8825971041d46d1dda511a223e01e6b799'
        assert srm_manifest['binding']['files_sha256']['solver/soil-fem-3d.exe'] == '5786e5885ab903d42a7702610257fa6ef1879230fe7c3493e9b3c9fbb0e027bd'
        srm_status = json.loads((srm_job / 'status.json').read_text(encoding='utf-8-sig'))
        if srm_status.get('state') == 'running' and srm_status.get('process_id') in owned:
            record = records[srm_status['process_id']]
            command = windows_argv(record.get('CommandLine') or '')
            assert len(command) >= 3 and command[1] == 'run-srm-search-job'
            assert Path(command[2]).is_absolute() and norm(Path(command[2]).resolve()) == norm(srm_job.resolve())
            assert record.get('ExecutablePath') and norm(Path(record['ExecutablePath']).resolve()) == norm((srm_job / 'solver/soil-fem-3d.exe').resolve())
            assert sha(Path(record['ExecutablePath'])) == '5786e5885ab903d42a7702610257fa6ef1879230fe7c3493e9b3c9fbb0e027bd'
            solver = srm_job / 'solver/soil-fem-3d.exe'
            process = subprocess.run([str(solver),'stop-srm-search-job',str(srm_job)],capture_output=True,text=True,timeout=30)
            assert process.returncode == 0, 'SRM cooperative stop request failed'
            for _ in range(36):
                value = json.loads((srm_job / 'status.json').read_text(encoding='utf-8-sig'))
                if value['state'] != 'running':break
                time.sleep(5)
            from preserve_srm_search_job import preserve_search_job
            preserved = None;preservation_errors = []
            for attempt in range(3):
                try:
                    preserved = preserve_search_job(srm_job,out/'srm-accepted-recovery-forks'/str(attempt),solver,'5786e5885ab903d42a7702610257fa6ef1879230fe7c3493e9b3c9fbb0e027bd')
                    break
                except (OSError,ValueError,subprocess.TimeoutExpired) as error:
                    preservation_errors.append(dict(attempt=attempt,error=str(error)))
                    if not any(word in str(error).lower() for word in ['changed','advanced']):break
            assert preserved and preserved['passed'], 'SRM full accepted snapshot not preserved; do not stop owned solver'
            cooperative[str(srm_job)] = dict(srm_search_job=True,stop_request_exit_code=process.returncode,
                final_state=value['state'],final_state_is_last_report_not_process_liveness=True,
                accepted_snapshot_preserved_before_owned_stop=preserved,preservation_attempt_errors=preservation_errors,
                discarded_unaccepted_trial_only=True)
    # registered_slope1200_srm_v2_job: exact path/manifest/executable only.
    # Ownership still requires held PID creation/image identity and exact argv.
    srm_job = root / 'runs/native-source-point-view-development-20261030-v1/slope1200-durable-SRM-v2/job'
    srm_manifest_path = srm_job / 'job.json'
    if srm_manifest_path.exists():
        srm_manifest = json.loads(srm_manifest_path.read_text(encoding='utf-8-sig'))
        assert sha(srm_manifest_path) == 'f5b210c081ebca389c9977162d7e3fa110f452408e91949ff2349fd84d54f6de', 'registered SRM manifest changed'
        assert srm_manifest.get('schema') == 'soil-fem-srm-search-job/1'
        assert srm_manifest['job_id'] == '3b933c98b62c8521665caa3579527733cde1bed0e69e7801ee70b3d4470f3ed8'
        assert srm_manifest['binding']['files_sha256']['solver/soil-fem-3d.exe'] == 'd8756fbb0d6db6f764746b15abbc0c3fc19bbcf1411a25c483a3e8369bcdf3f8'
        srm_status = json.loads((srm_job / 'status.json').read_text(encoding='utf-8-sig'))
        if srm_status.get('state') == 'running' and srm_status.get('process_id') in owned:
            record = records[srm_status['process_id']]
            command = windows_argv(record.get('CommandLine') or '')
            assert len(command) >= 3 and command[1] == 'run-srm-search-job'
            assert Path(command[2]).is_absolute() and norm(Path(command[2]).resolve()) == norm(srm_job.resolve())
            assert record.get('ExecutablePath') and norm(Path(record['ExecutablePath']).resolve()) == norm((srm_job / 'solver/soil-fem-3d.exe').resolve())
            assert sha(Path(record['ExecutablePath'])) == 'd8756fbb0d6db6f764746b15abbc0c3fc19bbcf1411a25c483a3e8369bcdf3f8'
            solver = srm_job / 'solver/soil-fem-3d.exe'
            process = subprocess.run([str(solver),'stop-srm-search-job',str(srm_job)],capture_output=True,text=True,timeout=30)
            assert process.returncode == 0, 'SRM cooperative stop request failed'
            for _ in range(36):
                value = json.loads((srm_job / 'status.json').read_text(encoding='utf-8-sig'))
                if value['state'] != 'running':break
                time.sleep(5)
            from preserve_srm_search_job import preserve_search_job
            preserved = None;preservation_errors = []
            for attempt in range(3):
                try:
                    preserved = preserve_search_job(srm_job,out/'srm-accepted-recovery-forks'/'slope1200-v2'/str(attempt),solver,'d8756fbb0d6db6f764746b15abbc0c3fc19bbcf1411a25c483a3e8369bcdf3f8')
                    break
                except (OSError,ValueError,subprocess.TimeoutExpired) as error:
                    preservation_errors.append(dict(attempt=attempt,error=str(error)))
                    if not any(word in str(error).lower() for word in ['changed','advanced']):break
            assert preserved and preserved['passed'], 'SRM full accepted snapshot not preserved; do not stop owned solver'
            cooperative[str(srm_job)] = dict(srm_search_job=True,stop_request_exit_code=process.returncode,
                final_state=value['state'],final_state_is_last_report_not_process_liveness=True,
                accepted_snapshot_preserved_before_owned_stop=preserved,preservation_attempt_errors=preservation_errors,
                discarded_unaccepted_trial_only=True)
    # Python references save completed phases only; transient RK4 state is not a checkpoint.
    candidates=[]
    for mode in ['cached100','cached200','original100']:
        candidates.extend((folders[0]/mode/'checkpoints').glob('phase-*.json'))
    for job in jobs:candidates.extend(job.glob('checkpoint*.json'))
    for job in jobs:candidates.extend(job.glob('recovery*.json'))
    for source in candidates:
        raw=source.read_bytes();json.loads(raw)
        target=out/'json-readable-saves'/source.relative_to(root/'runs')
        target.parent.mkdir(parents=True,exist_ok=True)
        with target.open('xb') as stream:stream.write(raw)
        saved[str(source)]=dict(copy=str(target),sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw))
    sys.path.insert(0,str(root/'tools'))
    from freeze_sources import freeze
    freeze(out)
    # Stop waiting drivers before their children so no subsequent phase can start.
    def depth(r):
        n=0;seen=set()
        while r['ParentProcessId'] in owned and r['ParentProcessId'] not in seen:
            seen.add(r['ParentProcessId']);r=records[r['ParentProcessId']];n+=1
        return n
    for h,r in sorted(handles,key=lambda pair:depth(pair[1])):
        if kernel.WaitForSingleObject(h,0)==0:
            continue
        terminated=kernel.TerminateProcess(h,1223)
        assert terminated or kernel.WaitForSingleObject(h,10000)==0
        assert kernel.WaitForSingleObject(h,10000)==0
        stopped.append(dict(pid=r['ProcessId'],image=r['ExecutablePath'],created=r['Created'],held_identity_verified=True))
    assert all(kernel.WaitForSingleObject(h,0)==0 for h,r in handles)
    assert all(sha(Path(row['copy']))==row['sha256'] for row in saved.values())
    receipt=dict(usage_percent=args.usage_percent,owned_processes_stopped=stopped,cooperative_rust_stop=cooperative,
                 completed_phase_copies=saved,current_partial_python_phase_not_claimed_resumable=True,
                 unrelated_processes_touched=False,p6_complete=False)
    with (out/'receipt.json').open('x') as stream:json.dump(receipt,stream,indent=2)
    print(json.dumps(dict(receipt=str(out/'receipt.json'),stopped=len(stopped),saved=len(saved))))
finally:
    for h,r in handles:kernel.CloseHandle(h)
