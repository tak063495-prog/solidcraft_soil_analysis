"""Only the child owned by this wrapper is stopped, before unlock."""
import json,os,signal,subprocess,sys,tempfile,time,unittest
from pathlib import Path
from unittest.mock import patch
from test_reference_resume import fixture
import resume_projected_sparse_study as resume
ROOT=Path(__file__).resolve().parents[1]
WORKER="""import os,pathlib,time
root=pathlib.Path(__file__).resolve().parents[1]
(root/'worker-pid').write_text(str(os.getpid()))
time.sleep(60)
"""
def sleeper(folder):
    fixture(folder);script=folder/'reference-sources/mcc_projected_sparse_run.py';script.write_text(WORKER)
    snapshot=resume.read(folder/'source-snapshot.json');snapshot['tools/'+script.name]=resume.sha(script)
    (folder/'source-snapshot.json').write_text(json.dumps(snapshot))
def alive(pid):
    if os.name=='nt':
        import ctypes
        kernel=ctypes.WinDLL('kernel32',use_last_error=True)
        kernel.OpenProcess.restype=ctypes.c_void_p
        kernel.OpenProcess.argtypes=[ctypes.c_uint32,ctypes.c_int,ctypes.c_uint32]
        kernel.WaitForSingleObject.argtypes=[ctypes.c_void_p,ctypes.c_uint32]
        kernel.CloseHandle.argtypes=[ctypes.c_void_p]
        handle=kernel.OpenProcess(0x100000,False,pid)
        if not handle:
            if ctypes.get_last_error()==87:return False
            raise ctypes.WinError(ctypes.get_last_error())
        try:
            result=kernel.WaitForSingleObject(handle,0)
            if result not in (0,258):raise ctypes.WinError(ctypes.get_last_error())
            return result==258
        finally:kernel.CloseHandle(handle)
    try:
        os.kill(pid,0)
        stat=Path('/proc')/str(pid)/'stat'
        return not stat.exists() or stat.read_text().split()[2]!='Z'
    except ProcessLookupError:return False

# Windows cannot send POSIX parent-only SIGTERM with os.kill: it hard-terminates
# the target. Run the same CLI and handlers, delivering a real CRT signal in the
# parent's main thread just before its owned-child wait. The child gets no signal.
WINDOWS_PARENT=r"""import pathlib,runpy,signal,sys,time
script,folder,number=sys.argv[1:]
folder=pathlib.Path(folder)
number=int(number)
def trace(frame,event,arg):
    if (event=='line' and frame.f_code.co_name=='run'
            and frame.f_code.co_filename==script
            and (frame.f_locals.get('progress') or {}).get('child_pid')):
        sys.settrace(None)
        deadline=time.monotonic()+10
        while not (folder/'deliver-parent-signal').exists():
            if time.monotonic()>deadline:raise RuntimeError('signal handshake timeout')
            time.sleep(.01)
        (folder/'parent-signal-delivered').write_text(str(number))
        signal.raise_signal(number)
    return trace
sys.argv=[script,str(folder)]
sys.settrace(trace)
runpy.run_path(script,run_name='__main__')
"""

class OwnedChild(unittest.TestCase):
    def parent_signal(self, parent_signal, cooperative=False):
        cooperative=cooperative or os.name=='nt'
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);sleeper(folder);child=None
            script=str(ROOT/'tools/resume_projected_sparse_study.py')
            command=([sys.executable,'-c',WINDOWS_PARENT,script,str(folder),str(int(parent_signal))]
                     if cooperative else [sys.executable,script,str(folder)])
            parent=subprocess.Popen(command,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
            try:
                deadline=time.monotonic()+10
                while not (folder/'worker-pid').exists() and time.monotonic()<deadline:time.sleep(.02)
                self.assertTrue((folder/'worker-pid').exists());child=int((folder/'worker-pid').read_text())
                self.assertTrue((folder/'.reference-resume.lock').exists())
                if cooperative:(folder/'deliver-parent-signal').write_text(str(int(parent_signal)))
                else:os.kill(parent.pid,parent_signal)
                parent.wait(timeout=15)
                self.assertNotEqual(parent.returncode,0)
                if cooperative:
                    self.assertEqual((folder/'parent-signal-delivered').read_text(),str(int(parent_signal)))
                progress=resume.read(Path(resume.read(folder/'latest-resume.json')['attempt'])/'progress.json')
                self.assertTrue(progress['interrupted'])
                self.assertTrue(progress['child_exit_confirmed'])
                self.assertEqual(progress['child_pid'],child)
                self.assertFalse(alive(child),'owned child survives parent SIGINT')
                self.assertFalse((folder/'.reference-resume.lock').exists())
                # Replace only this temporary study's worker with the short normal fixture.
                from test_reference_resume import WORKER as normal
                script=folder/'reference-sources/mcc_projected_sparse_run.py';script.write_text(normal)
                data=resume.read(folder/'source-snapshot.json');data['tools/'+script.name]=resume.sha(script);(folder/'source-snapshot.json').write_text(json.dumps(data))
                self.assertEqual(resume.run(folder),0)
            finally:
                if parent.poll() is None:parent.kill();parent.wait()
                if child and alive(child):os.kill(child,signal.SIGTERM if os.name=='nt' else signal.SIGKILL)
    def test_parent_only_sigint_reaps_child_before_unlock_and_allows_safe_retry(self):
        self.parent_signal(signal.SIGINT)
    def test_parent_only_sigterm_reaps_child_before_unlock_and_allows_safe_retry(self):
        self.parent_signal(signal.SIGTERM)

    def test_main_thread_runtime_sigint_reaps_real_child(self):
        self.parent_signal(signal.SIGINT,cooperative=True)
    def test_main_thread_runtime_sigterm_reaps_real_child(self):
        self.parent_signal(signal.SIGTERM,cooperative=True)

    def test_interrupt_during_launch_retains_the_handle_for_cleanup(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);sleeper(folder);original=resume.subprocess.Popen;children=[]
            def launch(*args,**kwargs):
                process=original(*args,**kwargs);children.append(process)
                signal.getsignal(signal.SIGINT)(signal.SIGINT,None)
                return process
            try:
                with patch.object(resume.subprocess,'Popen',side_effect=launch):
                    with self.assertRaises(KeyboardInterrupt):resume.run(folder)
                self.assertIsNotNone(children[0].returncode)
                self.assertFalse((folder/'.reference-resume.lock').exists())
            finally:
                for process in children:
                    if process.poll() is None:process.kill();process.wait()

    def injected(self,where,stop_fails=False):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);fixture(folder);events=[]
            class Process:
                pid=424242
                def __init__(self,*a,**kw):self.running=True;events.append('launch')
                def poll(self):return None if self.running else -15
                def terminate(self):
                    events.append('terminate');self_lock=folder/'.reference-resume.lock';self_test.assertTrue(self_lock.exists())
                    with self_test.assertRaisesRegex(ValueError,'lock exists'):resume.run(folder)
                    if stop_fails:raise OSError('cannot stop')
                    self.running=False
                def kill(self):
                    events.append('kill')
                    if stop_fails:raise OSError('cannot kill')
                    self.running=False
                def wait(self,timeout=None):
                    events.append('wait');self_test.assertTrue((folder/'.reference-resume.lock').exists())
                    if self.running:raise RuntimeError('wait injection')
                    return -15
            self_test=self;original=resume.atomic
            def atomic(path,value):
                if where=='status' and value.get('child_pid') and value.get('status')=='initial_running':raise OSError('status injection')
                return original(path,value)
            with patch.object(resume.subprocess,'Popen',Process),patch.object(resume,'atomic',side_effect=atomic):
                with self.assertRaises((OSError,RuntimeError)):resume.run(folder)
            self.assertIn('terminate',events)
            if stop_fails:
                self.assertTrue((folder/'.reference-resume.lock').exists());record=resume.read(folder/'.reference-resume.lock')
                self.assertEqual(record['child_pid'],424242);self.assertEqual(record['parent_pid'],os.getpid());self.assertTrue(record['interrupted'])
            else:self.assertIn('wait',events);self.assertFalse((folder/'.reference-resume.lock').exists())
    def test_wait_exception_stops_and_waits(self):self.injected('wait')
    def test_status_exception_after_launch_stops_and_waits(self):self.injected('status')
    def test_unconfirmed_child_exit_retains_lock_with_owner_record(self):self.injected('wait',True)
    def test_normal_child_failure_unlocks(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);fixture(folder)
            script=folder/'reference-sources/mcc_projected_sparse_run.py';script.write_text('raise SystemExit(7)')
            data=resume.read(folder/'source-snapshot.json');data['tools/'+script.name]=resume.sha(script);(folder/'source-snapshot.json').write_text(json.dumps(data))
            self.assertEqual(resume.run(folder),7);self.assertFalse((folder/'.reference-resume.lock').exists())

if __name__=='__main__':unittest.main()
