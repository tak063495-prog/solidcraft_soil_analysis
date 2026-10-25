"""Job ownership, immutable inputs and actual owned-tree stop checks."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tools'))
from native_job import atomic,copy_checked,request_cancel,run,status,stop_child,watch_host
from verify_acceptance_packet import sha


class NativeJob(unittest.TestCase):
    def test_live_host_pipe_allows_normal_interpreter_shutdown(self):
        # Regression for the actual native-host run: a daemon must not retain
        # BufferedReader's lock and crash Python when a job finishes normally.
        code=f"import sys,time;sys.path.insert(0,{str(ROOT/'tools')!r});from native_job import watch_host;watch_host(sys.stdin.buffer);time.sleep(.15)"
        options=dict(creationflags=subprocess.CREATE_NEW_PROCESS_GROUP|subprocess.CREATE_NO_WINDOW) if os.name=='nt' else dict(start_new_session=True)
        for _ in range(3):
            child=subprocess.Popen([sys.executable,'-c',code],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,**options)
            try:
                self.assertEqual(child.wait(timeout=10),0)
                self.assertFalse(child.stdin.closed)
            finally:
                child.stdin.close()
                if child.poll() is None:stop_child(child)
                child.stdout.close();child.stderr.close()

    def test_host_pipe_disconnect_requests_stop_and_accepts_no_commands(self):
        incoming,outgoing=os.pipe()
        with os.fdopen(incoming,'rb',buffering=0) as reader:
            event=watch_host(reader)
            os.write(outgoing,b'ignored input')
            self.assertFalse(event.wait(.1))
            os.close(outgoing)
            self.assertTrue(event.wait(3))

    def test_host_input_version_mismatch_creates_no_job_and_launches_nothing(self):
        from test_solidcraft_view import job_fixture
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp);job_fixture(folder/'fixture');capture=folder/'fixture/conformal'
            shutil.copyfile(capture/'native-capture.json',capture/'capture.json')
            with self.assertRaisesRegex(ValueError,'host input snapshot'):
                run(capture,folder/'refused-job',expected_input='0'*64)
            self.assertFalse((folder/'refused-job').exists())

    def test_cancel_is_identity_bound_idempotent_and_terminal_safe(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp)
            atomic(folder/'job.json',dict(schema='soil-fem-native-job/1',job_id='owned-job'))
            current=dict(schema='soil-fem-native-job-status/1',job_id='owned-job',state='running')
            atomic(folder/'status.json',current)
            self.assertTrue(request_cancel(folder)['requested']);self.assertTrue(request_cancel(folder)['requested'])
            self.assertEqual(json.loads((folder/'cancel-request.json').read_text())['job_id'],'owned-job')
            self.assertIn('last worker report',status(folder)['observation'])
            atomic(folder/'status.json',{**current,'state':'succeeded'})
            before=sha(folder/'cancel-request.json');self.assertFalse(request_cancel(folder)['requested'])
            self.assertEqual(sha(folder/'cancel-request.json'),before)

    def test_foreign_cancel_identity_is_refused_without_replacement(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp)
            atomic(folder/'job.json',dict(schema='soil-fem-native-job/1',job_id='owned'))
            atomic(folder/'status.json',dict(schema='soil-fem-native-job-status/1',job_id='owned',state='running'))
            atomic(folder/'cancel-request.json',dict(job_id='other'))
            before=sha(folder/'cancel-request.json')
            with self.assertRaisesRegex(ValueError,'identity'):request_cancel(folder)
            self.assertEqual(sha(folder/'cancel-request.json'),before)

    def test_snapshot_hash_and_existing_output_are_enforced(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp);source=folder/'input';source.write_bytes(b'original')
            digest=sha(source);copy_checked(source,folder/'snapshot',digest)
            source.write_bytes(b'changed')
            self.assertEqual((folder/'snapshot').read_bytes(),b'original')
            with self.assertRaises(FileExistsError):copy_checked(source,folder/'snapshot',sha(source))
            with self.assertRaisesRegex(ValueError,'changed during'):copy_checked(source,folder/'refused',digest)
            self.assertEqual(source.read_bytes(),b'changed')

    def test_actual_owned_child_and_descendant_stop_without_a_name_wide_kill(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp)
            grandchild="from pathlib import Path;import sys,time;p=Path(sys.argv[1]);f=p.open('a');\nwhile True:f.write('tick\\n');f.flush();time.sleep(.05)"
            parent="import json,subprocess,sys,time;from pathlib import Path;p=subprocess.Popen([sys.executable,'-u','-c',sys.argv[1],sys.argv[2]],stdin=subprocess.DEVNULL,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL);Path(sys.argv[3]).write_text(str(p.pid));time.sleep(600)"
            options=dict(creationflags=subprocess.CREATE_NEW_PROCESS_GROUP|subprocess.CREATE_NO_WINDOW) if os.name=='nt' else dict(start_new_session=True)
            child=subprocess.Popen([sys.executable,'-u','-c',parent,grandchild,str(folder/'heartbeat'),str(folder/'ready')],
                stdin=subprocess.DEVNULL,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,**options)
            try:
                deadline=time.monotonic()+10
                while not (folder/'heartbeat').exists() and time.monotonic()<deadline:time.sleep(.05)
                self.assertTrue((folder/'ready').is_file());self.assertTrue((folder/'heartbeat').is_file())
                before=(folder/'heartbeat').stat().st_size;time.sleep(.2)
                self.assertGreater((folder/'heartbeat').stat().st_size,before)
                stop_child(child);self.assertIsNotNone(child.poll())
                after=(folder/'heartbeat').stat().st_size;time.sleep(.4)
                self.assertEqual((folder/'heartbeat').stat().st_size,after)
                stop_child(child)  # Already terminal: do not signal another PID.
            finally:
                if child.poll() is None:stop_child(child)


if __name__=='__main__':unittest.main()
