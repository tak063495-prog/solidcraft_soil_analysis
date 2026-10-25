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
from native_job import atomic,copy_checked,request_cancel,run,status,stop_child,watch_host,frozen_profile
from verify_acceptance_packet import sha


class NativeJob(unittest.TestCase):
    def test_live_solver_survives_transient_status_read_denials(self):
        from unittest.mock import patch
        import native_stage_profile as profile
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp)
            expected={'schema':'soil-fem-stages-job-status/1','job_id':'owned-test','revision':7}
            (folder/'status.json').write_text(json.dumps(expected),encoding='utf-8')
            real=profile.read;attempts=[];reports=[]
            def reader(actual,limit):
                attempts.append((actual,limit))
                if len(attempts)<4:raise PermissionError('transient status sharing')
                return real(actual,limit)
            with patch.object(profile,'read',side_effect=reader):
                code,stopped=profile.command([sys.executable,'-c',"import time; print('ready',flush=True); time.sleep(.5)"],
                    ROOT,folder/'worker.log',lambda **kw:reports.append(kw),lambda:False,
                    lambda child:self.fail('transient read must not stop owned child'),phase='solving',solver_job=folder)
            self.assertEqual((code,stopped),(0,False))
            self.assertTrue(any(report.get('solver_report')==expected for report in reports))
            self.assertGreaterEqual(len(attempts),4)
            self.assertEqual(set(attempts),{(folder/'status.json',256*1024)})
            self.assertEqual(real(folder/'status.json'),expected)

    def test_solver_status_denial_is_bounded_and_corrupt_or_missing_progress_is_not_retried(self):
        from unittest.mock import patch
        import native_stage_profile as profile
        path=Path('owned/status.json')
        with patch.object(profile,'read',side_effect=PermissionError('permanent denial')) as reader,patch.object(profile.time,'sleep') as sleep:
            with self.assertRaises(PermissionError):profile.read_solver_status(path)
            self.assertEqual(reader.call_count,5)
            self.assertEqual(sleep.call_count,4)
            self.assertEqual({call.args for call in reader.call_args_list},{(path,256*1024)})
        for error in (ValueError('invalid JSON'),FileNotFoundError('missing progress')):
            with self.subTest(error=type(error).__name__),patch.object(profile,'read',side_effect=error) as reader,patch.object(profile.time,'sleep') as sleep:
                with self.assertRaises(type(error)):profile.read_solver_status(path)
                self.assertEqual(reader.call_count,1)
                sleep.assert_not_called()

    def test_persistent_status_failure_stops_owned_child_and_closes_output(self):
        from unittest.mock import patch
        import native_stage_profile as profile
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp);(folder/'status.json').write_text('{}',encoding='utf-8');stopped=[]
            def stop(child):
                stopped.append(child);stop_child(child)
            with patch.object(profile,'read',side_effect=PermissionError('persistent status denial')) as reader:
                with self.assertRaises(PermissionError):
                    profile.command([sys.executable,'-c',"import time; print('ready',flush=True); time.sleep(.5)"],
                        ROOT,folder/'worker.log',lambda **kw:None,lambda:False,stop,phase='solving',solver_job=folder)
            self.assertEqual(reader.call_count,5)
            self.assertEqual(len(stopped),1)
            self.assertIsNotNone(stopped[0].poll())
            self.assertTrue(stopped[0].stdout.closed)

    def test_atomic_status_retries_only_its_own_complete_temporary_on_windows_reader_contention(self):
        from unittest.mock import patch
        import native_workflow_snapshot as snapshot
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp);path=folder/'status.json';atomic(path,{'old':True})
            real_replace=os.replace;attempts=[]
            def replace(source,destination):
                attempts.append((source,destination))
                self.assertEqual(json.loads(source.read_text(encoding='utf-8')),{'new':True})
                if len(attempts)<4:raise PermissionError('reader held old status')
                real_replace(source,destination)
            with patch.object(snapshot.os,'replace',side_effect=replace):atomic(path,{'new':True})
            self.assertEqual(len(attempts),4);self.assertEqual(len(set(attempts)),1)
            self.assertEqual(json.loads(path.read_text(encoding='utf-8')),{'new':True})
            self.assertFalse(list(folder.glob('*.part')))
            with patch.object(snapshot.os,'replace',side_effect=PermissionError('persistent refusal')) as denied:
                with self.assertRaises(PermissionError):atomic(path,{'unpublished':True})
                self.assertEqual(denied.call_count,5)
            self.assertEqual(json.loads(path.read_text(encoding='utf-8')),{'new':True})
            self.assertEqual(len(list(folder.glob('*.part'))),1)

    def test_frozen_workflow_survives_source_edits_and_loads_its_own_later_profile(self):
        import native_workflow_snapshot as snapshot
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'source';(root/'tools').mkdir(parents=True)
            source=root/'tools/native_stage_profile.py';source.write_text('VALUE = 17\n',encoding='utf-8')
            bindings={'tools/native_stage_profile.py':sha(source)}
            runtime=snapshot.capture(root,Path(tmp)/'runtime',bindings)
            source.write_text('raise RuntimeError("edited workflow must never run")\n',encoding='utf-8')
            snapshot.verify(runtime,bindings)
            self.assertEqual(frozen_profile(runtime).VALUE,17)
            (runtime/'tools/native_stage_profile.py').write_text('VALUE = 99\n',encoding='utf-8')
            with self.assertRaisesRegex(ValueError,'captured workflow'):snapshot.verify(runtime,bindings)

    def test_workflow_inventory_cannot_copy_outside_either_root(self):
        import native_workflow_snapshot as snapshot
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'source';root.mkdir()
            secret=Path(tmp)/'outside';secret.write_text('outside',encoding='utf-8')
            for i,name in enumerate(('../outside','/outside','C:/outside','tools\\outside')):
                with self.subTest(name=name),self.assertRaisesRegex(ValueError,'safe relative'):
                    snapshot.capture(root,Path(tmp)/('runtime-'+str(i)),{name:sha(secret)})

    def saved_fixture(self, folder):
        import native_stage_profile as profile
        from test_solidcraft_view import job_fixture
        workflow=folder/'workflow';job_fixture(workflow)
        shutil.copytree(workflow/'conformal',folder/'capture')
        shutil.copyfile(folder/'capture/native-capture.json',folder/'capture/capture.json')
        atomic(workflow/'analysis/prepared-stage-input.json',dict(prepared=True))
        atomic(workflow/'conformal/prepared-model.json',dict(prepared=True))
        artifact=workflow/'analysis/solver-job/checkpoint.json';artifact.parent.mkdir()
        artifact.write_bytes(b'{"checkpoint":"data only, not executable"}')
        capture=json.loads((folder/'capture/capture.json').read_text(encoding='utf-8'))
        owner=dict(schema='soil-fem-native-job/2',job_id='saved-owned',input_sha256=capture['input_sha256'],
                   solver_sha256='a'*64,tool_sources_sha256={'current':'trusted'})
        atomic(folder/'job.json',owner)
        report=dict(schema='soil-fem-native-job-status/2',job_id=owner['job_id'])
        profile.saved(folder,owner,'solver_checkpoint',artifact,lambda **changes:report.update(changes))
        atomic(folder/'status.json',report)
        return owner,report

    def test_saved_resume_requires_current_host_solver_workflow_and_cad_identity(self):
        import native_stage_profile as profile
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp);owner,report=self.saved_fixture(folder)
            actual,receipt=profile.previous(folder,'a'*64,{'current':'trusted'},owner['input_sha256'])
            self.assertEqual(actual,owner);self.assertEqual(receipt['kind'],'solver_checkpoint')
            for solver,sources,cad in [('b'*64,{'current':'trusted'},owner['input_sha256']),
                                       ('a'*64,{'current':'changed'},owner['input_sha256']),
                                       ('a'*64,{'current':'trusted'},'f'*64)]:
                with self.subTest(solver=solver,sources=sources,cad=cad),self.assertRaises(ValueError):
                    profile.previous(folder,solver,sources,cad)

    def test_checkpoint_receipt_artifact_and_analysis_input_tampering_are_refused(self):
        import native_stage_profile as profile
        for name in ('checkpoint-receipt.json','workflow/analysis/solver-job/checkpoint.json',*profile.INPUT_FILES):
            with self.subTest(name=name),tempfile.TemporaryDirectory() as tmp:
                folder=Path(tmp);self.saved_fixture(folder)
                with (folder/name).open('ab') as stream:stream.write(b' changed')
                with self.assertRaises(ValueError):profile.previous(folder,'a'*64,{'current':'trusted'})

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
