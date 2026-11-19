"""Bounded child-only runner. Exit 2 stays incomplete and never means a pass."""
import argparse
import ctypes
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time


def windows_job(memory_bytes):
    from ctypes import wintypes as w
    class Basic(ctypes.Structure):
        _fields_=[('process_time',ctypes.c_int64),('job_time',ctypes.c_int64),('flags',w.DWORD),
            ('min_ws',ctypes.c_size_t),('max_ws',ctypes.c_size_t),('active',w.DWORD),
            ('affinity',ctypes.c_size_t),('priority',w.DWORD),('scheduling',w.DWORD)]
    class IO(ctypes.Structure):
        _fields_=[(name,ctypes.c_uint64) for name in ('read_ops','write_ops','other_ops','read_bytes','write_bytes','other_bytes')]
    class Extended(ctypes.Structure):
        _fields_=[('basic',Basic),('io',IO),('process_memory',ctypes.c_size_t),('job_memory',ctypes.c_size_t),
                  ('peak_process_memory',ctypes.c_size_t),('peak_job_memory',ctypes.c_size_t)]
    kernel=ctypes.WinDLL('kernel32',use_last_error=True)
    kernel.CreateJobObjectW.argtypes=[ctypes.c_void_p,w.LPCWSTR];kernel.CreateJobObjectW.restype=w.HANDLE
    kernel.SetInformationJobObject.argtypes=[w.HANDLE,ctypes.c_int,ctypes.c_void_p,w.DWORD]
    kernel.AssignProcessToJobObject.argtypes=[w.HANDLE,w.HANDLE]
    kernel.QueryInformationJobObject.argtypes=[w.HANDLE,ctypes.c_int,ctypes.c_void_p,w.DWORD,ctypes.c_void_p]
    kernel.TerminateJobObject.argtypes=[w.HANDLE,w.UINT];kernel.CloseHandle.argtypes=[w.HANDLE]
    handle=kernel.CreateJobObjectW(None,None)
    if not handle:raise ctypes.WinError(ctypes.get_last_error())
    limits=Extended();limits.basic.flags=0x2000|0x200 # kill on close + job commit memory
    limits.job_memory=memory_bytes
    if not kernel.SetInformationJobObject(handle,9,ctypes.byref(limits),ctypes.sizeof(limits)):
        kernel.CloseHandle(handle);raise ctypes.WinError(ctypes.get_last_error())
    return kernel,handle,Extended


def wait_linux_child(pid, timeout):
    deadline = time.monotonic() + timeout
    # Keep the leader waitable until after killpg, so its PID cannot be reused.
    while os.waitid(os.P_PID, pid, os.WEXITED | os.WNOHANG | os.WNOWAIT) is None:
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            return False
        time.sleep(min(0.05, remaining))
    return True


def stop_linux_group(child):
    pgid = child.pid  # start_new_session makes the child the group leader.
    try:
        os.killpg(pgid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    code = child.wait()
    # The subreaper adopts orphaned members, including zombies on systems
    # where PID 1 does not reap them. Reap only children of this group.
    while True:
        try:
            os.waitpid(-pgid, 0)
        except ChildProcessError:
            break
    # After reaping, this is only a presence check. Never signal a recycled PGID.
    try:
        os.killpg(pgid, 0)
    except ProcessLookupError:
        return code
    raise RuntimeError(f'could not verify disappearance of process group {pgid}')


def run_linux_group(command, timeout, **options):
    pr_set_child_subreaper = 36
    pr_get_child_subreaper = 37
    libc = ctypes.CDLL(None, use_errno=True)
    previous = ctypes.c_int()
    if libc.prctl(pr_get_child_subreaper, ctypes.byref(previous), 0, 0, 0) != 0:
        error = ctypes.get_errno()
        raise OSError(error, os.strerror(error))

    def set_subreaper(value):
        if libc.prctl(pr_set_child_subreaper, value, 0, 0, 0) != 0:
            error = ctypes.get_errno()
            raise OSError(error, os.strerror(error))

    interrupt = {signal.SIGINT}
    previous_mask = signal.pthread_sigmask(signal.SIG_BLOCK, interrupt)
    cap = options['preexec_fn']

    def prepare_child():
        signal.pthread_sigmask(signal.SIG_SETMASK, previous_mask)
        cap()

    options['preexec_fn'] = prepare_child
    try:
        set_subreaper(1)
        try:
            child = subprocess.Popen(command, **options)
            try:
                signal.pthread_sigmask(signal.SIG_SETMASK, previous_mask)
                timed_out = not wait_linux_child(child.pid, timeout)
            finally:
                try:
                    signal.pthread_sigmask(signal.SIG_BLOCK, interrupt)
                finally:
                    code = stop_linux_group(child)
        finally:
            set_subreaper(previous.value)
    finally:
        signal.pthread_sigmask(signal.SIG_SETMASK, previous_mask)
    return child.pid, code, timed_out


def run(command,output,timeout=1800,memory_gib=8,cwd=None):
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    if not 0<timeout<=7500 or not 0<memory_gib<=8:raise ValueError('bounded runner requires <=7500s and <=8GiB')
    budget=int(memory_gib*1024**3);job=None;p=None;timed=False;peak=None;start=time.monotonic()
    try:
        options={}
        if os.name=='nt':
            job=windows_job(budget);options['creationflags']=0x00000004 # CREATE_SUSPENDED
        else:
            import resource
            def cap():resource.setrlimit(resource.RLIMIT_AS,(budget,budget))
            options.update(start_new_session=True,preexec_fn=cap)
        with (output/'process.log').open('x',encoding='utf-8') as log:
            if sys.platform == 'linux':
                pid, code, timed = run_linux_group(
                    command, timeout, cwd=cwd, stdout=log, stderr=subprocess.STDOUT, **options)
            else:
                p=subprocess.Popen(command,cwd=cwd,stdout=log,stderr=subprocess.STDOUT,**options)
                pid=p.pid
                if job:
                    kernel,handle,extended=job
                    if not kernel.AssignProcessToJobObject(handle,int(p._handle)):
                        p.kill();p.wait();raise ctypes.WinError(ctypes.get_last_error())
                    ntdll=ctypes.WinDLL('ntdll');ntdll.NtResumeProcess.argtypes=[ctypes.c_void_p]
                    if ntdll.NtResumeProcess(int(p._handle))!=0:
                        kernel.TerminateJobObject(handle,1);p.wait();raise OSError('NtResumeProcess failed')
                try:code=p.wait(timeout=timeout)
                except subprocess.TimeoutExpired:
                    timed=True
                    if job:kernel.TerminateJobObject(handle,124)
                    else:os.killpg(p.pid,signal.SIGKILL)
                    code=p.wait()
                if job:
                    value=extended()
                    if not kernel.QueryInformationJobObject(handle,9,ctypes.byref(value),ctypes.sizeof(value),None):raise ctypes.WinError(ctypes.get_last_error())
                    peak=value.peak_job_memory
        status='incomplete' if timed or code==2 else ('done' if code==0 else 'failed')
        result=dict(status=status,exit_code=code,timed_out=timed,command=command,cwd=str(Path(cwd or Path.cwd()).resolve()),
            elapsed_seconds=time.monotonic()-start,memory_limit_bytes=budget,peak_job_committed_bytes=peak,
            memory_measurement='Windows job committed memory; Linux peak not measured',pid=pid,
            numerical_acceptance=False)
    except BaseException as exc:
        if p is not None and p.poll() is None:p.kill();p.wait()
        result=dict(status='failed',error=repr(exc),numerical_acceptance=False)
    finally:
        if job:job[0].CloseHandle(job[1])
    (output/'process.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    return result


def chain(items,output):
    """Only explicitly listed independent tasks; preserve every nonzero result."""
    output=Path(output);output.mkdir(parents=True,exist_ok=False);results=[]
    for i,item in enumerate(items):
        results.append(run(output=output/str(i),**item))
    report=dict(status='done' if all(r['status']=='done' for r in results) else 'incomplete',results=results,
                numerical_acceptance=False)
    (output/'chain.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--timeout',type=int,default=1800);p.add_argument('--memory-gib',type=float,default=8)
    p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args()
    command=a.command[1:] if a.command[:1]==['--'] else a.command
    if not command:p.error('child command required after --')
    result=run(command,a.output,a.timeout,a.memory_gib)
    print(json.dumps(result,ensure_ascii=False));raise SystemExit(0 if result['status']=='done' else (2 if result['status']=='incomplete' else 1))
