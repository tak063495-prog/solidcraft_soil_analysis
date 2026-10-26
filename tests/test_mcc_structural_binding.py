"""Source preservation cannot replace changed physical reference inputs."""
from pathlib import Path
import hashlib
import json
import sys
import tempfile
import unittest
import zipfile
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from verify_mcc_structural_study import verify_reference_bindings


class ReferenceBinding(unittest.TestCase):
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
