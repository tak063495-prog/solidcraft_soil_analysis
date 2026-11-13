"""Strict append-only compatibility: old document bytes still fail if changed."""
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('nc_physical_p601_appendix', ROOT/'interim/2026-12-06-p6-01-status-reconciliation/audit_status.py')
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


class PhysicalDocumentAppendix(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        source = ROOT/audit.PHYSICAL_APPENDICES_PATH
        self.manifest = self.root/audit.PHYSICAL_APPENDICES_PATH
        self.manifest.parent.mkdir(parents=True); self.manifest.write_bytes(source.read_bytes())
        self.tails = json.loads(source.read_bytes())['appendices']
        for path, tail in self.tails.items():
            p = self.root/path; p.parent.mkdir(parents=True, exist_ok=True)
            p.write_bytes(b'original prefix\n'+tail.encode('utf-8'))

    def test_exact_pinned_tails_are_accepted(self):
        self.assertEqual(audit.physical_document_appendices(self.root), self.tails)

    def test_resealed_or_whitespace_changed_manifest_is_rejected(self):
        self.manifest.write_bytes(self.manifest.read_bytes()+b' ')
        with self.assertRaisesRegex(ValueError, 'SHA-256 differs'):
            audit.physical_document_appendices(self.root)

    def test_added_document_text_cannot_hide_after_the_pinned_tail(self):
        path = self.root/'docs/LOW_PRELOAD_NC.ja.md'
        path.write_bytes(path.read_bytes()+b'forged')
        with self.assertRaisesRegex(ValueError, 'suffix differs'):
            audit.physical_document_appendices(self.root)

    def test_missing_manifest_is_rejected(self):
        self.manifest.unlink()
        with self.assertRaises(FileNotFoundError):
            audit.physical_document_appendices(self.root)

    def test_old_prefix_mutation_still_fails_original_hash_check(self):
        path = 'interim/2026-11-16-issue-triage/todo.ja.md'; original = b'original prefix\n'
        ledger = dict(base_commit=audit.BASE, sources={'old': dict(path=path, commit=audit.BASE,
            kind='markdown', sha256=audit.sha(original), status_keys=[],
            git_blob_sha1=hashlib.sha1(b'blob '+str(len(original)).encode()+b'\0'+original).hexdigest())},
            rows=[], evidence={}, claims=[], ci_records=[])
        def codes():
            return {f['code'] for f in audit.audit(self.root, ledger, verify_git=False)['findings']}
        self.assertNotIn('SOURCE_SHA256', codes())
        p = self.root/path; p.write_bytes(b'changed prefix\n'+self.tails[path].encode('utf-8'))
        self.assertIn('SOURCE_SHA256', codes())
        self.assertIn('SOURCE_BLOB', codes())


if __name__ == '__main__':
    unittest.main()
