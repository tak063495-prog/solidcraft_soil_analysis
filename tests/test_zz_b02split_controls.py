"""Fixed-seal replay and input-only checks for B-02, without a full solver run."""
import hashlib
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
import zipfile

FOLDER = Path(__file__).resolve().parents[1] / 'interim/2027-01-14-b02-fixed-offset-thickness-declaration-order'


class B02SplitControls(unittest.TestCase):
    def test_one_factor_inputs(self):
        previous = sys.path[:]
        try:
            sys.path.insert(0, str(FOLDER))
            from b02split_run import CONTROLS, inputs
            from b02follow_run import input_case
            for level in (0, 1):
                base_mesh, base_ground, base_stages, _ = inputs(level, .05, 'default')
                for thickness, order in CONTROLS:
                    mesh, ground, stages, mapping = inputs(level, thickness, order)
                    self.assertEqual(ground, base_ground)
                    self.assertEqual(stages['stages'], base_stages['stages'])
                    system = stages['beam_system']
                    self.assertEqual({p['attachment_offset_m'] for p in system['plates']}, {.025})
                    self.assertEqual({s['thickness_m'] for s in system['plate_sections']}, {thickness})
                    for key in base_mesh:
                        expected = base_mesh[key]
                        if order == 'elements' and key in ('tetrahedra','tetrahedra10','material_ids','region_ids'):
                            expected = expected[::-1]
                        if order == 'faces' and key == 'face_sets':
                            expected = {k: v[::-1] for k,v in expected.items()}
                        self.assertEqual(mesh[key], expected)
                    plates = base_stages['beam_system']['plates']
                    self.assertEqual(system['plates'], plates[::-1] if order == 'plates' else plates)
                    self.assertEqual(mapping['offset_over_thickness'], .025/thickness)
                    if order == 'default':
                        old_mesh, old_ground, old_stages, _ = input_case(level, thickness/.1, 'default')
                        self.assertEqual((mesh, ground), (old_mesh, old_ground))
                        for plate in old_stages['beam_system']['plates']:
                            self.assertEqual(plate['attachment_offset_m'], thickness/2)
                            plate['attachment_offset_m'] = .025
                        self.assertEqual(stages, old_stages)
        finally:
            sys.path[:] = previous

    def test_saved_thin96_original_gate(self):
        summary = json.loads((FOLDER/'measurements.json').read_bytes())
        index = json.loads((FOLDER/'archive-index.json').read_bytes())
        self.assertIs(summary['passed'], False)
        self.assertEqual(summary['accepted_states'], 480)
        h = hashlib.sha256()
        with tempfile.TemporaryFile() as stream:
            for part in index['parts']:
                data = (FOLDER/part['name']).read_bytes()
                self.assertEqual(len(data), part['size'])
                self.assertEqual(hashlib.sha256(data).hexdigest(), part['sha256'])
                stream.write(data)
                h.update(data)
            self.assertEqual(h.hexdigest(), index['archive_sha256'])
            self.assertEqual(h.hexdigest(), json.loads((FOLDER/'acceptance.json').read_bytes())['archive_sha256'])
            with zipfile.ZipFile(stream) as z:
                for order in ('default','elements','faces','plates'):
                    key = 'tet96-t0.05-' + order
                    path = f'controls/{key}/native-reassembly.json'
                    data = z.read(path)
                    self.assertEqual(hashlib.sha256(data).hexdigest(), index['entries'][path]['sha256'])
                    rows = json.loads(data)
                    peak = max(rows, key=lambda r: r['original_relative'])
                    recorded = summary['cases'][key]['residuals']['native']
                    self.assertEqual(peak['original_relative'], recorded['max_relative'])
                    self.assertGreater(recorded['max_relative'], 1e-10)
                    self.assertEqual(peak['accepted_index'], recorded['peak_accepted_index'])
                    magnitude = math.sqrt(math.fsum(x*x for x in peak['free_residual_components_kn']))
                    self.assertEqual(magnitude, recorded['peak_l2_kn'])
                    self.assertEqual(magnitude/peak['denominator_kn'], peak['relative_fixed_order'])
