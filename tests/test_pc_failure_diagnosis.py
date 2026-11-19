"""Light PC regressions. No equilibrium recalculation or acceptance relaxation."""
import copy
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import zipfile
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
HERE=ROOT/'interim/2027-01-08-pc-failure-diagnosis'
ISSUE_36=ROOT/'interim/2027-01-10-issue36-pcfail-runner-process-group'
sys.path.insert(0,str(HERE))
import pcfail_b02 as b02
import pcfail_evidence as evidence
from pcfail_hydraulic import Precision, load
from pcfail_run import chain,run
from pcfail_beam_build import generate


class PCFailureDiagnosis(unittest.TestCase):
    def test_reused_pid_reproduces_old_assert_and_preserves_exact_restart_checks(self):
        native=dict(pid=1156,invocation_id='first',fresh_popen=True,exit_code=0)
        resumed=dict(native,invocation_id='second')
        old=(ROOT/'interim/2026-12-24-b02-offset-composite-cause-diagnosis/b02cause_driver.py').read_text(encoding='utf-8')
        self.assertIn("assert process['pid'] != native['pid']",old)
        with self.assertRaises(AssertionError):
            exec("assert process['pid'] != native['pid']",{'process':resumed,'native':native})
        state={'physics':[1.,-0.], 'elapsed_seconds':1.}
        later=[dict(accepted_index=11,stress=[4.])]
        b02.check_resume(native,resumed,state,dict(state,elapsed_seconds=2.),later,later)
        for changed,changes in [(dict(state,physics=[1.,1.]),later),(state,[dict(accepted_index=11,stress=[5.])])]:
            with self.assertRaises(AssertionError):b02.check_resume(native,resumed,state,changed,changes,later)
        with self.assertRaises(AssertionError):b02.check_resume(native,native,state,state,later,later)

    def test_fresh_popen_records_distinct_invocations_even_if_pid_reused(self):
        class Process:
            pid=1156
            def wait(self,**kwargs):return 0
        with tempfile.TemporaryDirectory() as tmp,patch.object(b02.subprocess,'Popen',return_value=Process()):
            a=b02.command(['unused'],Path(tmp)/'a.log');b=b02.command(['unused'],Path(tmp)/'b.log')
        self.assertEqual(a['pid'],b['pid']);self.assertNotEqual(a['invocation_id'],b['invocation_id'])

    def fixture(self,folder):
        stream=io.BytesIO()
        with zipfile.ZipFile(stream,'w') as z:z.writestr('selected/input.json',b'{"value":1}\n')
        raw=stream.getvalue();part=folder/'evidence.zip.part000'
        index=dict(parts=[dict(name=part.name,size=len(raw),sha256=hashlib.sha256(raw).hexdigest())],
            archive_bytes=len(raw),archive_sha256=hashlib.sha256(raw).hexdigest(),
            entries={'selected/input.json':dict(size=12,sha256=hashlib.sha256(b'{"value":1}\n').hexdigest())})
        return part,raw,index

    def test_missing_evidence_is_incomplete_nonpass_and_no_extract(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp);part,raw,index=self.fixture(folder)
            result=evidence.verify(folder,index,folder/'extract')
            self.assertEqual(result['status'],'incomplete');self.assertFalse(result['passed'])
            self.assertEqual(result['parts'][0]['status'],'missing');self.assertFalse((folder/'extract').exists())
            part.write_bytes(raw[:-1]);self.assertEqual(evidence.verify(folder,index)['status'],'failed')

    def test_valid_bytes_are_verified_without_changing_acceptance_and_mutation_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp);part,raw,index=self.fixture(folder);part.write_bytes(raw)
            result=evidence.verify(folder,index,folder/'out')
            self.assertEqual(result['status'],'verified');self.assertFalse(result['passed'])
            self.assertEqual((folder/'out/selected/input.json').read_bytes(),b'{"value":1}\n')
            bad=copy.deepcopy(index);bad['entries']['selected/input.json']['sha256']='0'*64
            with self.assertRaisesRegex(ValueError,'member SHA'):evidence.verify(folder,bad)

    def test_selective_fixed_blob_restore_stays_sparse_with_posix_nul_paths(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);folder=root/'interim/fixture';folder.mkdir(parents=True)
            part,raw,index=self.fixture(folder);part.write_bytes(raw)
            (root/'src.py').write_bytes(b'print(1)\n')
            def git(*args,**kw):return subprocess.check_output(['git','-C',str(root),*args],stderr=subprocess.PIPE,**kw)
            git('init');git('config','user.name','test');git('config','user.email','test@example.invalid')
            git('config','core.autocrlf','false');git('add','.');git('commit','-m','fixture')
            commit=git('rev-parse','HEAD').decode().strip()
            git('sparse-checkout','init','--no-cone')
            git('sparse-checkout','set','--no-cone','--stdin',input=b'/*\r\n!/interim/**/*.part*\r\n')
            self.assertFalse(part.exists())
            with patch.object(evidence,'ROOT',root):
                result=evidence.fetch(folder,index,commit)
            self.assertEqual(result['status'],'available');self.assertEqual(part.read_bytes(),raw)
            paths={p.decode('utf-8').replace('\\','/') for p in git('ls-files','-z').split(b'\0') if p}
            self.assertIn('interim/fixture/evidence.zip.part000',paths)
            self.assertEqual(git('diff','--name-only','-z'),b'')

    def test_cli_missing_dependency_records_before_solver(self):
        import pcfail_cli
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp);_,_,index=self.fixture(folder)
            with patch.object(pcfail_cli,'manifest',return_value=(folder,index)):
                self.assertEqual(pcfail_cli.run(folder/'result',folder/'nonexistent-solver'),2)
            result=json.loads((folder/'result/evidence-status.json').read_text())
            self.assertEqual(result['status'],'incomplete');self.assertFalse((folder/'result/controls').exists())

    def test_chain_records_nonzero_and_continues_independent_tasks(self):
        with tempfile.TemporaryDirectory() as tmp:
            items=[dict(command=[sys.executable,'-c','raise SystemExit(2)']),
                   dict(command=[sys.executable,'-c','print("independent task")'])]
            result=chain(items,Path(tmp)/'chain')
            self.assertEqual([r['status'] for r in result['results']],['incomplete','done'])
            self.assertEqual(result['status'],'incomplete');self.assertFalse(result['numerical_acceptance'])

    def test_timeout_retains_child_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            result=run([sys.executable,'-c','import time; time.sleep(5)'],Path(tmp)/'run',timeout=.1)
            self.assertTrue(result['timed_out']);self.assertEqual(result['status'],'incomplete')

    def test_precision_override_is_isolated_and_reference_policy_unchanged(self):
        original=np.longdouble;module=load(np.float64)
        self.assertIs(np.longdouble,original);self.assertIs(module.np.longdouble,np.float64)
        self.assertIs(module.hyd.np.longdouble,np.float64)
        self.assertEqual(Precision(np.float64).finfo(np.float64).nmant,52)
        source=(ROOT/'interim/2026-11-22-wet-joint-transient-water/hydraulic_reference.py').read_text()
        for marker in ('range(60)','range(20)','if merit < .01:','lin > 1e-10'):
            self.assertIn(marker,source)

    def test_beam_observation_changes_only_failure_branch_in_isolated_copy(self):
        original=(ROOT/'src/beams.rs').read_bytes()
        with tempfile.TemporaryDirectory() as tmp:
            output=Path(tmp)/'observer';binding=generate(output)
            self.assertTrue(binding['production_unchanged'])
            s=(output/'src/beams.rs').read_text()
            self.assertIn('PCFAIL_BEAM_TRACE',s)
            self.assertEqual(s.count('beam interior true residual exceeded tolerance'),1)
            self.assertIn('relative > residual_tolerance',s)
        self.assertEqual((ROOT/'src/beams.rs').read_bytes(),original)

    def test_sealed_numerical_failures_are_not_accepted_or_recomputed(self):
        from pcfail_beam_report import norm,residual
        self.assertEqual(norm([3.,4.]),5.)
        np.testing.assert_array_equal(residual(np.array([[1e16,1.,-1e16]]),[1.,1.,1.],[0.]),[1.])
        value=json.loads((HERE/'measurements.json').read_text(encoding='utf-8'))
        self.assertFalse(value['passed']);self.assertFalse(value['all_acceptance_complete'])
        roof=value['roof']
        self.assertEqual([r['complete'] for r in roof['paths']],[True,True,False])
        self.assertGreater(roof['linear_observation']['native_relative'],1e-10)
        self.assertGreater(roof['linear_observation']['fixed_order_original_relative'],1e-10)
        self.assertEqual([r['execution_complete'] for r in value['hydraulic']['precision_cross']],
                         [True,False,False,True])
        index=json.loads((HERE/'archive-index.json').read_text(encoding='utf-8'))
        digest=hashlib.sha256()
        for part in index['parts']:
            path=HERE/part['name']
            self.assertEqual(path.stat().st_size,part['size'])
            with path.open('rb') as stream:
                self.assertEqual(hashlib.file_digest(stream,'sha256').hexdigest(),part['sha256'])
            with path.open('rb') as stream:
                for raw in iter(lambda:stream.read(1024*1024),b''):digest.update(raw)
        self.assertEqual(digest.hexdigest(),index['archive_sha256'])

    def test_ci_explicit_archive_selection_uses_posix_nul_and_crlf_rules(self):
        from pcfail_sparse import rules
        names=[p.relative_to(ROOT).as_posix() for p in sorted(HERE.glob('evidence.zip.part*'))]
        self.assertTrue(names)
        inputs=('\0'.join(names)+'\0').encode('utf-8')
        for group in ('python-core-0','python-core-1','rust-integration-1','python-c01-cause'):
            selected=rules(ROOT/'.github/workflows/fem-ci.yml',group)
            expected=inputs if group.startswith('python-core-') else b''
            for newline in ('\n','\r\n'):
                with tempfile.TemporaryDirectory() as tmp:
                    path=Path(tmp)/'rules';path.write_bytes((newline.join(selected)+newline).encode('utf-8'))
                    result=subprocess.run(['git','sparse-checkout','check-rules','-z','--no-cone','--rules-file',str(path)],
                        cwd=ROOT,input=inputs,capture_output=True)
                self.assertEqual(result.returncode,0,result.stderr)
                self.assertEqual(result.stdout,expected)

    def test_document_newlines_and_unresolved_status(self):
        marker='今回の限定確認を全体完了へ変更しない。\n'
        text=(HERE/'README.ja.md').read_text(encoding='utf-8')
        self.assertTrue(text.endswith(marker))
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'README.ja.md';path.write_bytes(text.replace('\n','\r\n').encode())
            self.assertTrue(path.read_text(encoding='utf-8').endswith(marker))

    def test_issue36_old_and_repaired_runner_digests(self):
        record = json.loads((ISSUE_36/'acceptance.json').read_text(encoding='utf-8'))
        runner = HERE/'pcfail_run.py'
        name = runner.relative_to(ROOT).as_posix()
        fix = record['runner']
        self.assertEqual(fix['path'], name)
        revision = f"{fix['before_commit']}:{name}"
        args = ['git', '-C', str(ROOT)]
        blob = subprocess.check_output([*args, 'rev-parse', revision]).decode().strip()
        self.assertEqual(blob, fix['before_blob'])
        old = subprocess.check_output([*args, 'show', revision])
        before = hashlib.sha256(old).hexdigest()
        self.assertEqual(before, fix['before_sha256'])
        for manifest in ('final-code-binding.json', 'validation-code-binding.json'):
            with self.subTest(manifest=manifest):
                archived = json.loads((HERE/manifest).read_text(encoding='utf-8'))
                self.assertEqual(archived[name], before)
        self.assertEqual(hashlib.sha256(runner.read_bytes()).hexdigest(), fix['after_sha256'])
        self.assertNotEqual(before, fix['after_sha256'])

    def test_issue36_documents_end_with_unresolved_status_on_lf_and_crlf(self):
        ending = 'A-01〜A-04、B-01〜B-03、C-01〜C-04、NC-01、P6 は未達、passed=false。\n'
        with tempfile.TemporaryDirectory(prefix='issue 36 docs ') as tmp:
            for name in ('README.ja.md', 'REPRODUCE.ja.md'):
                text = (ISSUE_36/name).read_text(encoding='utf-8')
                self.assertTrue(text.endswith(ending), name)
                for newline in ('\n', '\r\n'):
                    with self.subTest(name=name, newline=repr(newline)):
                        copy_path = Path(tmp)/name
                        copy_path.write_bytes(text.replace('\n', newline).encode('utf-8'))
                        self.assertTrue(copy_path.read_text(encoding='utf-8').endswith(ending))

if __name__=='__main__':unittest.main()
