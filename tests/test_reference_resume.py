"""Frozen phase reuse must not change inputs, sources, results or backend."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
import resume_projected_sparse_study as resume

WORKER='''import hashlib,json,pathlib,sys
p=sys.argv;phase=p[1];inputs=list(map(pathlib.Path,p[2:6]));output=pathlib.Path(p[6]);root=output.parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
with (root/'called.txt').open('a') as f:f.write(phase+'\\n')
r=dict(schema='soil-fem-projected-volume-mcc-diagnostic/1' if phase=='initial' else 'soil-fem-projected-mcc-cycle-diagnostic/1',
 reference_steps=int(p[p.index('--steps')+1]),load_fraction=1.,
 reference_backend='mcc_projected_sparse_reference.SparseProjectedSystem',
 input_and_source_sha256={str(q):sha(q) for q in inputs+[pathlib.Path(__file__)]})
if phase=='cycle':r['initial_reference_sha256']=sha(root/'initial.json')
with output.open('x') as f:json.dump(r,f)
'''

def fixture(folder):
    source=folder/'reference-sources';source.mkdir();(folder/'inputs').mkdir()
    for n in resume.INPUTS:(folder/'inputs'/n).write_text('{}',encoding='utf-8')
    (source/'mcc_projected_sparse_run.py').write_text(WORKER,encoding='utf-8')
    (source/'requirements-reference.lock.txt').write_text('scipy==1.18.1\n',encoding='utf-8')
    (source/'requirements.lock.txt').write_text('numpy==2.5.4\n',encoding='utf-8')
    for name in ('mcc_projected_sparse_reference.py','mcc_projected_volume_diagnostic.py',
                 'mcc_projected_cycle_diagnostic.py','mcc_rate_fem.py','mcc_rate_midpoint.py',
                 'fem_reference.py','verify_high_order.py','high_order_io.py'):
        (source/name).write_text('# frozen fixture module\n',encoding='utf-8')
    snapshot={('tools/'+p.name if p.suffix=='.py' else p.name):resume.sha(p) for p in source.iterdir()}
    (folder/'source-snapshot.json').write_text(json.dumps(snapshot),encoding='utf-8')
    (folder/'progress.json').write_text(json.dumps(dict(status='failed',steps=2)),encoding='utf-8')

class ReferenceResume(unittest.TestCase):
    def test_new_initial_and_cycle_are_verified_and_completed_phases_are_reused(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);fixture(p)
            self.assertEqual(resume.run(p),0)
            before={n:resume.sha(p/n) for n in ('initial.json','cycle.json','progress.json')}
            self.assertEqual(resume.run(p),0)
            self.assertEqual(before,{n:resume.sha(p/n) for n in before})
            self.assertEqual((p/'called.txt').read_text(),'initial\ncycle\n')
            self.assertFalse((p/'.reference-resume.lock').exists())

    def test_saved_initial_is_reused_without_replacing_its_history(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);fixture(p);resume.run(p)
            initial=(p/'initial.json').read_bytes();(p/'cycle.json').unlink()
            self.assertEqual(resume.run(p),0)
            self.assertEqual((p/'initial.json').read_bytes(),initial)
            self.assertEqual((p/'called.txt').read_text(),'initial\ncycle\ncycle\n')

    def test_changed_source_input_and_backend_are_rejected_before_launch(self):
        for change in ('source','input','backend'):
            with self.subTest(change=change),tempfile.TemporaryDirectory() as d:
                p=Path(d);fixture(p);resume.run(p);(p/'cycle.json').unlink()
                called=(p/'called.txt').read_bytes()
                if change=='source':
                    with (p/'reference-sources/mcc_projected_sparse_run.py').open('a') as f:f.write('\n# changed\n')
                elif change=='input':(p/'inputs/ground.json').write_text('{"changed":true}',encoding='utf-8')
                else:
                    data=resume.read(p/'initial.json');data['reference_backend']='unknown.Reference'
                    (p/'initial.json').write_text(json.dumps(data),encoding='utf-8')
                with self.assertRaises(ValueError):resume.run(p)
                self.assertEqual((p/'called.txt').read_bytes(),called)

    def test_partial_result_and_existing_lock_are_preserved(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);fixture(p);(p/'initial.json').write_text('{partial',encoding='utf-8')
            with self.assertRaises(ValueError):resume.run(p)
            self.assertEqual((p/'initial.json').read_text(),'{partial')
            (p/'.reference-resume.lock').write_text('{"pid":123}',encoding='utf-8')
            with self.assertRaises(ValueError):resume.run(p)
            self.assertTrue((p/'.reference-resume.lock').exists())

    def test_missing_frozen_module_is_rejected_before_running(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);fixture(p);(p/'reference-sources/fem_reference.py').unlink()
            with self.assertRaisesRegex(ValueError,'missing'):resume.run(p)
            self.assertFalse((p/'called.txt').exists())

    def test_explicit_historical_program_requires_its_completed_cycle_binding(self):
        with tempfile.TemporaryDirectory() as d:
            prior=Path(d)/'prior';current=Path(d)/'current';prior.mkdir();current.mkdir()
            fixture(prior);fixture(current);resume.run(prior)
            program=prior/'reference-sources/mcc_projected_sparse_run.py'
            version=dict(same_programs_as_completed_200=True,steps=2,
                previous_completed_cycle_sha256=resume.sha(prior/'cycle.json'),
                files_sha256={str(program):resume.sha(program)})
            (current/'reference-version-binding.json').write_text(json.dumps(version),encoding='utf-8')
            snapshot=resume.read(current/'source-snapshot.json');snapshot['tools/mcc_projected_sparse_run.py']='0'*64
            (current/'source-snapshot.json').write_text(json.dumps(snapshot),encoding='utf-8')
            self.assertEqual(resume.run(current),0)
            with (prior/'cycle.json').open('a') as f:f.write('\n ')
            with self.assertRaisesRegex(ValueError,'cycle changed'):resume.run(current)

if __name__=='__main__':unittest.main()
