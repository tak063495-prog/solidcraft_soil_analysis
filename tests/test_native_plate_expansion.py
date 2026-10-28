"""Compact CAD declarations stay bound through meshing, restart and display."""
import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
import fem_reference as ref
from plate_surfaces import attach
from prepare_model_record import record
from prepare_native_plates import derive, prepare, verify


class NativePlateExpansion(unittest.TestCase):
    def inputs(self):
        root = Path(__file__).resolve().parents[1]
        mesh = json.loads((root/'examples/hardening-column-tet10-mesh.json').read_text())
        mesh['face_sets'] = {'wall': [list(f) for f, owners in ref.topology(mesh).items()
                                     if len(owners) == 1 and all(mesh['nodes'][n][2] == 1. for n in f)]}
        spec = dict(schema='soil-fem-plate-surfaces/1',
                    sections=[dict(id='wall-section', young_kpa=3e7, poisson=.2, thickness_m=.2,
                                   unit_weight_kn_m3=24., drilling_factor=1e-4)],
                    surfaces=[dict(id='wall', face_set='wall', owner_body='UpperClay',
                                   section_id='wall-section', initially_active=False)])
        case = dict(beam_system=None, stages=[dict(kind='install_structures', members=['wall'])],
                    plate_surfaces=spec)
        return mesh, case

    def source(self, folder):
        source = folder/'source'
        source.mkdir()
        mesh, case = self.inputs()
        for name, value in [('mesh.json', mesh), ('stages.json', case), ('ground.json', {}),
                            ('native-capture.json', dict(input_sha256='a'*64))]:
            (source/name).write_text(json.dumps(value), encoding='utf-8')
        record(source)
        return source, mesh, case

    def test_compact_case_derives_exact_existing_plate_binding_without_source_mutation(self):
        mesh, case = self.inputs()
        before = copy.deepcopy((mesh, case))
        plain = copy.deepcopy(case)
        spec = plain.pop('plate_surfaces')
        expected = attach(mesh, plain, spec)
        result = derive(mesh, case)
        self.assertEqual(result[0], mesh)
        self.assertEqual(result[2], expected[2])
        old_nodes = expected[1]['beam_system'].pop('nodes')
        new_nodes = copy.deepcopy(result[1]['beam_system']['nodes'])
        new_case = copy.deepcopy(result[1])
        del new_case['beam_system']['nodes']
        self.assertEqual(new_case, expected[1])
        for old, new in zip(old_nodes, new_nodes, strict=True):
            soil = expected[0]['node_sets'][old.pop('soil_node_set')][0]
            self.assertTrue(new.pop('match_soil_node'))
            self.assertEqual(new, old)
            self.assertEqual([i for i, point in enumerate(mesh['nodes'])
                              if sum((a-b)**2 for a, b in zip(point, new['position_m'])) <= 1e-16], [soil])
        self.assertEqual((mesh, case), before)
        self.assertNotIn('plate_surfaces', result[1])
        self.assertEqual(result[1]['stages'][0]['members'], result[2]['member_groups']['wall'])

    def test_generated_connections_do_not_consume_the_native_user_group_budget(self):
        mesh, case = self.inputs()
        mesh['node_sets']['plate:mesh-node:user-authored'] = [0]
        while len(mesh['node_sets']) < 256:
            mesh['node_sets'][f'user-group-{len(mesh["node_sets"])}'] = [0]
        plain = copy.deepcopy(case)
        spec = plain.pop('plate_surfaces')
        expanded, _, _ = attach(mesh, plain, spec)
        self.assertGreater(len(expanded['node_sets']), 256)
        bound, structure, _ = derive(mesh, case)
        self.assertEqual(bound, mesh)
        self.assertEqual(len(bound['node_sets']), 256)
        self.assertTrue(all(n.get('match_soil_node') for n in structure['beam_system']['nodes']))

    def test_coordinate_ambiguity_is_rejected_before_a_native_attachment_can_change(self):
        mesh, case = self.inputs()
        _, _, binding = derive(mesh, case)
        index = binding['generated_soil_attachment_nodes'][0]
        mesh['nodes'].append(copy.deepcopy(mesh['nodes'][index]))
        with self.assertRaisesRegex(ValueError, 'coordinate attachment is ambiguous'):
            derive(mesh, case)

    def test_prepared_inventory_and_restart_copy_keep_reproducible_source_linkage(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            source, mesh, case = self.source(folder)
            before = {p.name: p.read_bytes() for p in source.iterdir()}
            output = folder/'bound'
            receipt = prepare(source, output)
            expanded_mesh = json.loads((output/'mesh.json').read_text())
            expanded_case = json.loads((output/'stages.json').read_text())
            self.assertEqual(verify(source, output, expanded_mesh, expanded_case, case), receipt)
            self.assertEqual({p.name: p.read_bytes() for p in source.iterdir()}, before)
            inventory = json.loads((output/'prepared-model.json').read_text())
            self.assertFalse(inventory['equilibrium_verified'])
            for name, digest in inventory['files_sha256'].items():
                self.assertEqual(hashlib.sha256((output/name).read_bytes()).hexdigest(), digest)
            # Resume copies this inventory rather than relying on the original expansion directory.
            copied = folder/'resume-ground'
            copied.mkdir()
            for name in inventory['files_sha256']:
                (copied/name).write_bytes((output/name).read_bytes())
            self.assertEqual(verify(source, copied, expanded_mesh, expanded_case, case), receipt)

    def test_regeneration_rejects_mesh_stage_binding_and_receipt_tampering(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            source, _, case = self.source(folder)
            output = folder/'bound'
            prepare(source, output)
            mesh = json.loads((output/'mesh.json').read_text())
            expanded = json.loads((output/'stages.json').read_text())
            wrong_mesh = copy.deepcopy(mesh)
            wrong_mesh['nodes'][0][0] += .01
            with self.assertRaisesRegex(ValueError, 'regenerated'):
                verify(source, output, wrong_mesh, expanded, case)
            wrong_case = copy.deepcopy(expanded)
            wrong_case['beam_system']['plate_sections'][0]['young_kpa'] *= 2
            with self.assertRaisesRegex(ValueError, 'regenerated'):
                verify(source, output, mesh, wrong_case, case)
            for name in ('unexpanded-stages.json', 'plate-surface-binding.json', 'plate-expansion.json'):
                with self.subTest(name=name):
                    path = output/name
                    original = path.read_bytes()
                    path.write_text('{}', encoding='utf-8')
                    with self.assertRaises(ValueError):
                        verify(source, output, mesh, expanded, case)
                    path.write_bytes(original)
            wrong_source = copy.deepcopy(case)
            wrong_source['plate_surfaces']['sections'][0]['thickness_m'] = .4
            with self.assertRaisesRegex(ValueError, 'captured source'):
                verify(source, output, mesh, expanded, wrong_source)

    def test_changed_source_inventory_and_unsafe_names_fail_before_creating_output(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            source, _, _ = self.source(folder)
            original = (source/'mesh.json').read_bytes()
            (source/'mesh.json').write_text('{}', encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'input changed'):
                prepare(source, folder/'bad')
            self.assertFalse((folder/'bad').exists())
            (source/'mesh.json').write_bytes(original)
            manifest = json.loads((source/'prepared-model.json').read_text())
            manifest['files_sha256']['../escape.json'] = '0'*64
            (source/'prepared-model.json').write_text(json.dumps(manifest), encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'safe plain'):
                prepare(source, folder/'unsafe')
            self.assertFalse((folder/'unsafe').exists())

    def test_invalid_declarations_and_unverified_couplings_are_not_silently_dropped(self):
        mesh, case = self.inputs()
        for invalid in (None, [], 'wall-spec.json'):
            wrong = copy.deepcopy(case)
            wrong['plate_surfaces'] = invalid
            with self.assertRaisesRegex(ValueError, 'specification object'):
                derive(mesh, wrong)
        for coupling in ('joint_system', 'anchor_system'):
            wrong = copy.deepcopy(case)
            wrong[coupling] = {}
            with self.assertRaisesRegex(ValueError, 'separate verification'):
                derive(mesh, wrong)

    def test_legacy_case_is_unchanged_and_cannot_claim_an_unrequested_expansion(self):
        mesh, case = self.inputs()
        del case['plate_surfaces']
        expected = copy.deepcopy((mesh, case))
        output_mesh, output_case, binding = derive(mesh, case)
        self.assertEqual((output_mesh, output_case), expected)
        self.assertIsNone(binding)
        with tempfile.TemporaryDirectory() as temp:
            prepared = Path(temp)
            self.assertIsNone(verify(prepared, prepared, mesh, case, case))
            (prepared/'plate-expansion.json').write_text('{}', encoding='utf-8')
            with self.assertRaises(ValueError):
                verify(prepared, prepared, mesh, case, case)


if __name__ == '__main__':
    unittest.main()
