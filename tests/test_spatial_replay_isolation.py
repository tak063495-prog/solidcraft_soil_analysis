"""Historical spatial evidence must use pinned sources, never new JOINT code."""
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
HERE=ROOT/'interim/2026-11-10-nc-spatial-release'
SPEC=importlib.util.spec_from_file_location('spatial_replay_isolation',HERE/'remeasure.py')
replay=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(replay)


class SpatialReplayIsolation(unittest.TestCase):
    def test_utf8_isolated_child_imports_only_the_selected_snapshot(self):
        with tempfile.TemporaryDirectory() as folder:
            snapshot=Path(folder);(snapshot/'joint_reference.py').write_text("marker='sealed'\n")
            code="import sys;sys.path.insert(0,sys.argv[1]);import joint_reference;print(joint_reference.marker)"
            result=subprocess.run(replay.isolated_command('-c',code,snapshot),capture_output=True,text=True,check=True)
            self.assertEqual(result.stdout,'sealed\n')
            self.assertEqual(replay.isolated_command('script.py')[1:4],['-I','-X','utf8'])

    def test_missing_dependencies_are_copied_only_at_the_historical_sha(self):
        expected=json.loads((HERE/'source-fingerprints.json').read_text())
        with tempfile.TemporaryDirectory() as folder:
            snapshot=Path(folder);copied=replay.complete_snapshot(snapshot)
            self.assertIn('tools/fem_reference.py',copied)
            self.assertIn('historical-spatial-dependencies',copied['tools/fem_reference.py']['source'])
            for name,record in copied.items():
                self.assertEqual(hashlib.sha256((snapshot/name).read_bytes()).hexdigest(),expected[name])
                self.assertEqual(record['sha256'],expected[name])
            # The changed live JOINT must be supplied by the sealed archive;
            # dependency completion cannot silently fill it with current code.
            self.assertFalse((snapshot/'tools/joint_reference.py').exists())

    def test_tampered_fallback_is_never_copied(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);(root/'tools').mkdir();(root/'tools/fem_reference.py').write_bytes(b'changed')
            retained=root/'interim/2026-11-11-joint-crack-tip/historical-spatial-dependencies';retained.mkdir(parents=True)
            (retained/'fem_reference.py').write_bytes(b'also changed')
            with patch.object(replay,'ROOT',root):
                copied=replay.complete_snapshot(root/'snapshot')
            self.assertNotIn('tools/fem_reference.py',copied)
            self.assertFalse((root/'snapshot/tools/fem_reference.py').exists())


if __name__=='__main__':unittest.main()
