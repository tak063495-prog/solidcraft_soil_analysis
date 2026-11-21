"""Reject rehashed interior fields that contact-only and endpoint checks accept."""
import copy
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import audit_wet_joint_bulk_history as bulk
from wet_joint_cycle_reference import serializable


class BulkAcceptedHistory(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(); cls.addClassCleanup(cls.temp.cleanup)
        cls.folder = Path(cls.temp.name)/'history'
        cls.values = bulk.history.inputs()
        sink = bulk.history.HistoryWriter(cls.folder, chunk_size=4)
        cls.result = serializable(bulk.history.cut_volume_reference(*cls.values, 4, 8, 7, history_sink=sink))
        sink.close()

    def test_all_accepted_bulk_water_fields_reassemble_without_acceptance(self):
        report = bulk.audit(self.folder, self.result, self.values)
        self.assertTrue(report['bulk_water_consistent'])
        self.assertEqual(len(report['reassembled_rows']), 28)
        self.assertLess(report['bulk_water_maximum_absolute_differences']['force'], 1e-10)
        self.assertFalse(report['passed'])
        self.assertFalse(report['native_comparison'])
        self.assertFalse(report['independent_bulk_plastic_solution'])

    def mutation(self, key):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)/'history'; shutil.copytree(self.folder, folder)
            manifest = json.loads((folder/'manifest.json').read_text())
            result = copy.deepcopy(self.result)
            chunk = manifest['chunks'][0]; path = folder/chunk['file']
            with np.load(path, allow_pickle=False) as source:
                fields = {k: bulk.history.decode(source[k], manifest['encoding']) for k in source.files}
            if key == 'factor':
                chunk['rows'][1]['factor'] += .1
                result['path'][1]['factor'] += .1
            else:
                fields[key].reshape(len(fields[key]), -1)[1, 0] += 1.
                with path.open('wb') as stream:
                    np.savez_compressed(stream, **{k: bulk.history.encode(v) for k, v in fields.items()})
                chunk['sha256'] = bulk.history.sha(path.read_bytes()); chunk['bytes'] = path.stat().st_size
            (folder/'manifest.json').write_text(json.dumps(manifest))
            # Hash, all old field invariants, endpoints and contact propagation
            # still pass: rejection must come from the added bulk/water audit.
            bulk.contact_audit(folder, result, self.values)
            with self.assertRaises(AssertionError):
                bulk.audit(folder, result, self.values)

    def test_rehashed_interior_stress_is_rejected(self):
        self.mutation('stress')

    def test_rehashed_interior_strain_is_rejected(self):
        self.mutation('strain')

    def test_rehashed_interior_pore_pressure_is_rejected(self):
        self.mutation('pore')

    def test_rehashed_interior_water_load_is_rejected(self):
        self.mutation('water_load')

    def test_rehashed_interior_reaction_is_rejected(self):
        self.mutation('reaction')

    def test_consistently_relabelled_interior_factor_is_rejected(self):
        self.mutation('factor')

    def test_other_volume_metadata_is_not_certified_by_this_audit(self):
        result = copy.deepcopy(self.result)
        result['volume_rule'] = 'input_gauss_jacobi'
        with self.assertRaisesRegex(AssertionError, 'fixed affine-cut'):
            bulk.audit(self.folder, result, self.values)


if __name__ == '__main__':
    unittest.main()
