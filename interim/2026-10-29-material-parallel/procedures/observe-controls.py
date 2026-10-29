"""Read CPU/peak memory through verified process handles; never signal a solver."""
from pathlib import Path
import ctypes,ctypes.wintypes as w,datetime,hashlib,json,os,sys,time
p=Path(__file__).resolve().parent
with (p/'telemetry-driver-child.json').open('x') as f:json.dump(dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),f)
k=ctypes.WinDLL('kernel32',use_last_error=True)
k.OpenProcess.argtypes=[w.DWORD,w.BOOL,w.DWORD];k.OpenProcess.restype=w.HANDLE
k.CloseHandle.argtypes=[w.HANDLE];k.WaitForSingleObject.argtypes=[w.HANDLE,w.DWORD]
k.QueryFullProcessImageNameW.argtypes=[w.HANDLE,w.DWORD,w.LPWSTR,ctypes.POINTER(w.DWORD)]
k.GetProcessTimes.argtypes=[w.HANDLE,*([ctypes.POINTER(w.FILETIME)]*4)]
class Counters(ctypes.Structure):
    _fields_=[('cb',w.DWORD),('faults',w.DWORD),*[(n,ctypes.c_size_t) for n in ['peak_ws','ws','peak_paged_pool','paged_pool','peak_nonpaged_pool','nonpaged_pool','commit','peak_commit','private']]]
k.K32GetProcessMemoryInfo.argtypes=[w.HANDLE,ctypes.POINTER(Counters),w.DWORD]
read=lambda q:json.loads(q.read_text('utf-8-sig'));sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
ticks=lambda t:(t.dwHighDateTime<<32)|t.dwLowDateTime
reports=[]
for n in [6,12]:
    out=p/('n'+str(n));receipt=out/'solve-child.json'
    while not receipt.exists():
        if (p/'pause-requested').exists():raise SystemExit(2)
        for q in (p/'full-quality').glob('*-terminal.json'):
            if read(q).get('exit_code',0)!=0:raise RuntimeError('candidate QA failed before controls; no telemetry accepted')
        time.sleep(1)
    child=read(receipt);job=out/'solver-job';solver=job/'solver/soil-fem-3d.exe'
    assert Path(child['executable']).resolve()==solver.resolve()
    expected=read(job/'job.json')['binding']['files_sha256']['solver/soil-fem-3d.exe'];assert sha(solver)==expected
    h=k.OpenProcess(0x1000|0x0400|0x0010|0x100000,False,child['pid']);assert h,'control process already gone; do not invent resource measurements'
    samples=0;peak_ws=0;peak_commit=0
    try:
        image=ctypes.create_unicode_buffer(32768);size=w.DWORD(len(image))
        assert k.QueryFullProcessImageNameW(h,0,image,ctypes.byref(size)) and Path(image.value).resolve()==solver.resolve()
        times=[w.FILETIME() for _ in range(4)];assert k.GetProcessTimes(h,*[ctypes.byref(t) for t in times])
        created=ticks(times[0]);expected_time=datetime.datetime.fromisoformat(child['created_utc']).timestamp()
        assert abs(created/1e7-11644473600-expected_time)<5
        with (out/'process-observations.jsonl').open('x') as f:
            while k.WaitForSingleObject(h,0)==258:
                if (p/'pause-requested').exists():raise SystemExit(2)
                memory=Counters();memory.cb=ctypes.sizeof(memory)
                if not k.K32GetProcessMemoryInfo(h,ctypes.byref(memory),memory.cb):
                    if k.WaitForSingleObject(h,0)==0:break
                    raise OSError(ctypes.get_last_error(),'could not sample owned solver memory')
                assert k.GetProcessTimes(h,*[ctypes.byref(t) for t in times])
                peak_ws=max(peak_ws,memory.peak_ws);peak_commit=max(peak_commit,memory.peak_commit);samples+=1
                f.write(json.dumps(dict(observed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),cpu_seconds=(ticks(times[2])+ticks(times[3]))/1e7,
                    working_set_bytes=memory.ws,observed_peak_working_set_bytes=memory.peak_ws,observed_peak_commit_bytes=memory.peak_commit))+'\n');f.flush()
                k.WaitForSingleObject(h,1000)
        assert k.WaitForSingleObject(h,0)==0 and k.GetProcessTimes(h,*[ctypes.byref(t) for t in times])
        assert samples>0
        report=dict(subdivisions=n,elements=6*n**3,elapsed_seconds=(ticks(times[1])-created)/1e7,cpu_seconds=(ticks(times[2])+ticks(times[3]))/1e7,
            observed_peak_working_set_bytes=peak_ws,observed_peak_commit_bytes=peak_commit,samples=samples,solver_sha256=expected,
            actual_process_image_and_creation_held=True,no_process_signalled=True,concurrent_other_jobs=True,isolated_throughput_benchmark=False,p6_complete=False)
        with (out/'process-terminal.json').open('x') as f:json.dump(report,f,indent=2)
        reports.append(report);print(json.dumps(report),flush=True)
    finally:k.CloseHandle(h)
with (p/'telemetry-terminal.json').open('x') as f:json.dump(dict(passed=True,read_only=True,cases=reports,no_process_signalled=True,p6_complete=False),f,indent=2)
