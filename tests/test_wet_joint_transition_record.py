from pathlib import Path
import tempfile
import unittest


RECORD = (
    Path(__file__).resolve().parents[1]
    / "interim/2026-11-21-issue37-38-wet-joint-transition"
)
STATUS = "FEM の全項目は未達、passed=false。\n"


class WetJointTransitionRecordTests(unittest.TestCase):
    def test_document_status_with_lf_and_crlf(self):
        with tempfile.TemporaryDirectory() as directory:
            for name in ("README.ja.md", "REPRODUCE.ja.md"):
                with self.subTest(document=name):
                    source = RECORD / name
                    raw = source.read_bytes()
                    self.assertNotIn(b"\r", raw)
                    self.assertTrue(raw.endswith(STATUS.encode("utf-8")))
                    windows_copy = Path(directory) / name
                    windows_copy.write_bytes(raw.replace(b"\n", b"\r\n"))
                    self.assertTrue(windows_copy.read_text(encoding="utf-8").endswith(STATUS))
