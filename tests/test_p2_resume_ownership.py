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
    try:
        os.kill(pid,0)
        stat=Path('/proc')/str(pid)/'stat'
        return not stat.exists() or stat.read_text().split()[2]!='Z'
    except ProcessLookupError:return False

class OwnedChild(unittest.TestCase):
    def parent_signal(self, parent_signal):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);sleeper(folder);child=None
            parent=subprocess.Popen([sys.executable,str(ROOT/'tools/resume_projected_sparse_study.py'),str(folder)],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
            try:
                deadline=time.monotonic()+10
                while not (folder/'worker-pid').exists() and time.monotonic()<deadline:time.sleep(.02)
                self.assertTrue((folder/'worker-pid').exists());child=int((folder/'worker-pid').read_text())
                self.assertTrue((folder/'.reference-resume.lock').exists())
                os.kill(parent.pid,parent_signal);parent.wait(timeout=15)
                self.assertFalse(alive(child),'owned child survives parent SIGINT')
                self.assertFalse((folder/'.reference-resume.lock').exists())
                # Replace only this temporary study's worker with the short normal fixture.
                from test_reference_resume import WORKER as normal
                script=folder/'reference-sources/mcc_projected_sparse_run.py';script.write_text(normal)
                data=resume.read(folder/'source-snapshot.json');data['tools/'+script.name]=resume.sha(script);(folder/'source-snapshot.json').write_text(json.dumps(data))
                self.assertEqual(resume.run(folder),0)
            finally:
                if parent.poll() is None:parent.kill();parent.wait()
                if child and alive(child):os.kill(child,signal.SIGKILL)
    @unittest.skipUnless(os.name=='posix','parent-only signals are POSIX regressions')
    def test_parent_only_sigint_reaps_child_before_unlock_and_allows_safe_retry(self):
        self.parent_signal(signal.SIGINT)
    @unittest.skipUnless(os.name=='posix','parent-only signals are POSIX regressions')
    def test_parent_only_sigterm_reaps_child_before_unlock_and_allows_safe_retry(self):
        self.parent_signal(signal.SIGTERM)

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
