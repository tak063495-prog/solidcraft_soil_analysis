import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

folder=Path(__file__).resolve().parent
root=folder.parents[1]
sys.path.insert(0,str(root/'tools'))
from native_job import stop_child

def module(name):
    spec=importlib.util.spec_from_file_location(name,folder/(name+'-profile.py'))
    value=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value

original=module('original')
candidate=module('candidate')

class StatusSharing(unittest.TestCase):
    def test_old_live_command_aborts_on_first_transient_status_open_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)
            (path/'status.json').write_text(json.dumps({'schema':'soil-fem-stages-job-status/1'}),encoding='utf-8')
            stopped=[]
            def stop(child):
                stopped.append(child)
                stop_child(child)
            with patch.object(original,'read',side_effect=PermissionError('transient status sharing')):
                with self.assertRaises(PermissionError):
                    original.command([sys.executable,'-c',"import time; print('ready',flush=True); time.sleep(.5)"],
                        root,path/'old.log',lambda **kw:None,lambda:False,stop,phase='solving',solver_job=path)
            self.assertEqual(len(stopped),1)
            stopped[0].stdout.close()

    def test_live_command_survives_transient_read_denials_without_stopping_child(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)
            expected={'schema':'soil-fem-stages-job-status/1','job_id':'owned-test','revision':7}
            (path/'status.json').write_text(json.dumps(expected),encoding='utf-8')
            real=candidate.read
            attempts=[]
            reports=[]
            def reader(actual,limit):
                attempts.append((actual,limit))
                if len(attempts)<4:
                    raise PermissionError('transient status sharing')
                return real(actual,limit)
            with patch.object(candidate,'read',side_effect=reader):
                code,stopped=candidate.command([sys.executable,'-c',"import time; print('ready',flush=True); time.sleep(.5)"],
                    root,path/'new.log',lambda **kw:reports.append(kw),lambda:False,
                    lambda child:self.fail('transient read must not stop owned child'),phase='solving',solver_job=path)
            self.assertEqual((code,stopped),(0,False))
            self.assertTrue(any(report.get('solver_report')==expected for report in reports))
            self.assertGreaterEqual(len(attempts),4)
            self.assertEqual(set(attempts),{(path/'status.json',256*1024)})
            self.assertEqual(real(path/'status.json'),expected)

    def test_permanent_denial_is_bounded_and_invalid_or_missing_progress_is_not_retried(self):
        path=Path('owned/status.json')
        with patch.object(candidate,'read',side_effect=PermissionError('permanent denial')) as reader,patch.object(candidate.time,'sleep') as sleep:
            with self.assertRaises(PermissionError):
                candidate.read_solver_status(path)
            self.assertEqual(reader.call_count,5)
            self.assertEqual(sleep.call_count,4)
            self.assertEqual({call.args for call in reader.call_args_list},{(path,256*1024)})
        for error in (ValueError('non-finite or invalid JSON'),FileNotFoundError('missing progress')):
            with self.subTest(error=type(error).__name__),patch.object(candidate,'read',side_effect=error) as reader,patch.object(candidate.time,'sleep') as sleep:
                with self.assertRaises(type(error)):
                    candidate.read_solver_status(path)
                self.assertEqual(reader.call_count,1)
                sleep.assert_not_called()

    def test_permanent_read_failure_stops_own_child_and_closes_its_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)
            (path/'status.json').write_text('{}',encoding='utf-8')
            stopped=[]
            def stop(child):
                stopped.append(child)
                stop_child(child)
            with patch.object(candidate,'read',side_effect=PermissionError('persistent status denial')) as reader:
                with self.assertRaises(PermissionError):
                    candidate.command([sys.executable,'-c',"import time; print('ready',flush=True); time.sleep(.5)"],
                        root,path/'persistent.log',lambda **kw:None,lambda:False,stop,phase='solving',solver_job=path)
            self.assertEqual(reader.call_count,5)
            self.assertEqual(len(stopped),1)
            self.assertIsNotNone(stopped[0].poll())
            self.assertTrue(stopped[0].stdout.closed)

if __name__=='__main__':
    unittest.main(verbosity=2)
