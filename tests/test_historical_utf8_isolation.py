"""Isolated historical replay must read UTF-8 even on an ANSI Windows locale."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('historical_utf8_replay', ROOT/'interim/2026-11-10-nc-mean-stress-release/remeasure.py')
replay = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(replay)


class HistoricalUTF8Isolation(unittest.TestCase):
    def test_explicit_utf8_retains_isolation_and_ignores_conflicting_environment(self):
        with tempfile.TemporaryDirectory() as temp:
            env = dict(os.environ, PYTHONUTF8='0', PYTHONIOENCODING='ascii', PYTHONPATH=temp)
            source = 'import json,sys;print(json.dumps([sys.flags.isolated,sys.flags.utf8_mode,sys.argv[1] in sys.path]))'
            child = subprocess.run(replay.isolated_command('-c', source, temp), env=env,
                                   capture_output=True, text=True, encoding='utf-8', check=True)
            self.assertEqual(json.loads(child.stdout), [1, 1, False])

    def test_default_path_read_roundtrips_the_archived_diagnostic_characters(self):
        text = 'NC “非収束” λ=1 未達'
        with self.assertRaises(UnicodeDecodeError):
            text.encode('utf-8').decode('cp1252')
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/'reference.json'
            path.write_text(json.dumps(dict(failure_reason=text), ensure_ascii=False), encoding='utf-8')
            source = 'import json,sys;from pathlib import Path;print(json.dumps(json.loads(Path(sys.argv[1]).read_text())))'
            child = subprocess.run(replay.isolated_command('-c', source, path),
                                   env=dict(os.environ, PYTHONUTF8='0'), capture_output=True,
                                   text=True, encoding='utf-8', check=True)
            self.assertEqual(json.loads(child.stdout), dict(failure_reason=text))


if __name__ == '__main__':
    unittest.main()
