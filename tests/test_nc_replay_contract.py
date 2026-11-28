"""Small archives exercise historical source isolation and fail-closed bindings.

The real mean/spatial archives are still exercised by their existing four tests
each. These synthetic archives certify no numerical result.
"""
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tarfile
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]


def replay_module(kind):
    path = ROOT/f'interim/2026-11-10-nc-{kind}-release/remeasure.py'
    spec = importlib.util.spec_from_file_location('contract_'+kind, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fixture(root, kind, corruption=None):
    here = root/f'interim/2026-11-10-nc-{kind}-release'
    here.mkdir(parents=True)
    script = f'interim/2026-11-10-nc-{kind}-release/measure.py'
    code = """import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'tools'))
from contract_dependency import VALUE
def measure(evidence, allow_platform_float_differences=False):
    return dict(marker=VALUE, metric=1.25, numerical_acceptance=False)
"""
    entries = {
        'compressed-jsonl.json': b'{}',
        'sources/'+script: code.encode('utf-8'),
        'sources/tools/contract_dependency.py': "VALUE='封印済み λ=1 未達'\n".encode('utf-8'),
        'sources/tools/binding_eol.py': b'# sealed LF source\n',
    }
    if corruption == 'missing_source':
        entries.pop('sources/'+script)
    manifest = dict(entries={name: dict(size=len(data), sha256=hashlib.sha256(data).hexdigest())
                             for name, data in entries.items()})
    if corruption == 'entry_hash':
        manifest['entries']['sources/tools/contract_dependency.py']['sha256'] = '0'*64
    if corruption == 'inventory':
        entries['unlisted.txt'] = b'unlisted'
    if corruption == 'missing_entry':
        entries.pop('sources/tools/binding_eol.py')
    stream = io.BytesIO()
    with tarfile.open(fileobj=stream, mode='w:gz') as archive:
        for name, data in entries.items():
            entry = tarfile.TarInfo(name)
            entry.size = len(data)
            archive.addfile(entry, io.BytesIO(data))
    data = stream.getvalue()
    part = here/'evidence.tar.gz.part000'
    part.write_bytes(data + (b'changed' if corruption == 'part' else b''))
    manifest.update(archive_size=len(data), archive_sha256=hashlib.sha256(data).hexdigest(),
                    parts=[dict(name=part.name, size=len(data), sha256=hashlib.sha256(data).hexdigest())])
    (here/'archive-sha256.json').write_text(json.dumps(manifest), encoding='utf-8')
    expected = dict(marker='封印済み λ=1 未達', metric=1.25, numerical_acceptance=False)
    if corruption == 'discrete':
        expected['numerical_acceptance'] = True
    if corruption == 'float':
        expected['metric'] = 2.5
    (here/'measurements.json').write_text(json.dumps(expected, indent=2)+'\n', encoding='utf-8')
    (here/'source-fingerprints.json').write_text('{}', encoding='utf-8')
    # Both live imports would raise if the historical subprocess used them.
    (here/'measure.py').write_text("raise RuntimeError('live measurement executed')\n", encoding='utf-8')
    (root/'tools').mkdir()
    (root/'tools/contract_dependency.py').write_text("raise RuntimeError('live dependency executed')\n", encoding='utf-8')
    (root/'tools/binding_eol.py').write_bytes(b'# sealed LF source\r\n')
    return here


class HistoricalReplayContract(unittest.TestCase):
    def test_changed_checkout_uses_verified_snapshot_and_reports_eol_separately(self):
        for kind in ('mean-stress', 'spatial'):
            with self.subTest(kind=kind), tempfile.TemporaryDirectory() as folder:
                root = Path(folder)
                here = fixture(root, kind)
                module = replay_module(kind)
                with patch.object(module, 'ROOT', root), patch.object(module, 'HERE', here):
                    result = module.run(root/'replay')
                self.assertTrue(result['archive_verified'])
                self.assertEqual(result['verified_entries'], 4)
                self.assertTrue(result['measurements_byte_identical'])
                self.assertTrue(result['discrete_evidence_exact'])
                self.assertFalse(result['numerical_acceptance'])
                self.assertFalse(result['current_checkout_numerical_certification'])
                self.assertEqual(len(result['checkout_source_differences']), 2)
                self.assertEqual([x['path'] for x in result['checkout_source_EOL_differences']],
                                 ['sources/tools/binding_eol.py'])
                self.assertEqual(json.loads((root/'replay/measurements.json').read_text())['marker'],
                                 '封印済み λ=1 未達')

    def test_corrupt_parts_sources_inventory_and_measurements_are_rejected(self):
        for kind in ('mean-stress', 'spatial'):
            for corruption, message in (
                    ('part', 'archive part changed'),
                    ('entry_hash', 'archive entry differs'),
                    ('inventory', 'unsafe/unlisted archive entry'),
                    ('missing_entry', 'archive inventory differs'),
                    ('missing_source', 'sealed measurement source missing'),
                    ('discrete', 'remeasured discrete evidence differs'),
                    ('float', 'remeasured JSON differs')):
                with self.subTest(kind=kind, corruption=corruption), tempfile.TemporaryDirectory() as folder:
                    root = Path(folder)
                    here = fixture(root, kind, corruption)
                    module = replay_module(kind)
                    with patch.object(module, 'ROOT', root), patch.object(module, 'HERE', here):
                        with self.assertRaisesRegex(ValueError, message):
                            module.run(root/'replay')
                    self.assertFalse((root/'replay/verification.json').exists())


if __name__ == '__main__':
    unittest.main()
