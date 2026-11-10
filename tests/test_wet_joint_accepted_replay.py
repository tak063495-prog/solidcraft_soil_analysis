"""Reject plausible, rehashed interior histories that endpoint checks miss."""
import copy
import io
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import audit_wet_joint_accepted_history as replay


class AcceptedJointReplay(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.temp.cleanup)
        cls.folder = Path(cls.temp.name)/'original'
        cls.values = replay.history.inputs()
        sink = replay.history.HistoryWriter(cls.folder, chunk_size=3)
        cls.result = replay.history.cut_volume_reference(*cls.values, 4, 8, 7, history_sink=sink)
        sink.close()

    def mutate(self, folder, kind):
        manifest = json.loads((folder/'manifest.json').read_text())
        chunk = manifest['chunks'][0]
        path = folder/chunk['file']
        with np.load(io.BytesIO(path.read_bytes()), allow_pickle=False) as source:
            fields = {key: replay.history.decode(source[key], manifest['encoding']) for key in source.files}
        if kind == 'slip':
            fields['plastic_slip'][1, 0, 0] += 1e-5
            fields['accumulated'][1, 0] += 1.
        elif kind == 'accumulated':
            fields['accumulated'][1, 0] += 1.
        else:
            fields['plastic_dissipation'][1, 0] += 1.
            fields['endpoint_work'][1, 0] += 1.
        with path.open('wb') as destination:
            np.savez_compressed(destination, **{key: replay.history.encode(value) for key, value in fields.items()})
        chunk['sha256'] = replay.history.sha(path.read_bytes())
        chunk['bytes'] = path.stat().st_size
        (folder/'manifest.json').write_text(json.dumps(manifest))

    def test_all_accepted_states_replay_without_claiming_convergence(self):
        report = replay.audit(self.folder, self.result, self.values)
        self.assertTrue(report['history_consistent'])
        self.assertEqual(report['accepted_states'], 28)
        self.assertEqual(report['contact_points'], 32)
        self.assertFalse(report['passed'])
        self.assertFalse(report['full_history_convergence_accepted'])

    def test_rehashed_plastic_vector_change_is_rejected(self):
        self.check_mutation('slip')

    def test_rehashed_accumulation_change_is_rejected(self):
        self.check_mutation('accumulated')

    def test_rehashed_work_split_change_is_rejected(self):
        self.check_mutation('work')

    def check_mutation(self, kind):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)/'mutated'
            shutil.copytree(self.folder, folder)
            self.mutate(folder, kind)
            # The original fields, chunks and endpoints still pass. The new
            # check must detect lost history propagation, not a hash mismatch.
            replay.history.check_history(folder, self.result, self.values)
            with self.assertRaisesRegex(AssertionError, 'accepted history propagation'):
                replay.audit(folder, self.result, self.values)


if __name__ == '__main__':
    unittest.main()
