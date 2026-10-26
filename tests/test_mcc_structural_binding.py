"""Source preservation cannot replace changed physical reference inputs."""
from pathlib import Path
import hashlib
import copy
import json
import sys
import tempfile
import unittest
import zipfile
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from verify_mcc_structural_study import verify_reference_bindings, physical_case


class ReferenceBinding(unittest.TestCase):
    def physical_fixture(self):
        root=Path(__file__).resolve().parents[1]
        return json.loads((root/'examples/stages-mixed-mcc.json').read_text(encoding='utf-8'))

    def test_numerical_integration_and_load_steps_preserve_physical_reference_case(self):
        case=self.physical_fixture();changed=copy.deepcopy(case)
        changed['settings'].update(mcc_integration='midpoint',initial_increment=.03125,
            material_maximum_substep=2.5e-5,relative_force_tolerance=1e-10)
        changed['stages'][0]['initial_increment']=.03125
        self.assertEqual(physical_case(case),physical_case(changed))

    def test_kinematics_unknown_settings_and_pressure_policy_are_physical_changes(self):
        case=self.physical_fixture();canonical=physical_case(case)
        for key,value in [('kinematics','oedometer'),('unknown_physics',True)]:
            changed=copy.deepcopy(case);changed['settings'][key]=value
            self.assertNotEqual(canonical,physical_case(changed))
        changed=copy.deepcopy(case)
        changed['pressure_on_original_surface_only']=not case.get('pressure_on_original_surface_only',False)
        self.assertNotEqual(canonical,physical_case(changed))

    def fixture(self,folder):
        root=Path(__file__).resolve().parents[1]
        source=root/'tools/nonexistent-preserved-reference-test.py'
        self.assertFalse(source.exists())
        original=b'original independent reference source\n'
        mesh=folder/'mesh.json';mesh.write_text('{}')
        archive=folder/'sources.zip'
        with zipfile.ZipFile(archive,'x') as out:
            out.writestr('tools/'+source.name,original);out.writestr('mesh.json','{}')
        reference=dict(input_and_source_sha256={str(source):hashlib.sha256(original).hexdigest(),
            str(mesh):hashlib.sha256(mesh.read_bytes()).hexdigest()})
        return reference,archive,mesh

    def test_exact_original_source_archive_resolves_code_that_has_since_changed(self):
        with tempfile.TemporaryDirectory() as tmp:
            reference,archive,_=self.fixture(Path(tmp))
            verify_reference_bindings(reference,archive)
            with self.assertRaises(ValueError):verify_reference_bindings(reference)

    def test_archive_never_substitutes_a_changed_mesh(self):
        with tempfile.TemporaryDirectory() as tmp:
            reference,archive,mesh=self.fixture(Path(tmp));mesh.write_text('{"changed":true}')
            with self.assertRaises(ValueError):verify_reference_bindings(reference,archive)

    def test_different_archived_source_bytes_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            reference,archive,_=self.fixture(Path(tmp))
            bad=Path(tmp)/'different.zip'
            with zipfile.ZipFile(bad,'x') as out:
                out.writestr('tools/nonexistent-preserved-reference-test.py','different equations')
            with self.assertRaises(ValueError):verify_reference_bindings(reference,bad)


if __name__=='__main__':unittest.main()
