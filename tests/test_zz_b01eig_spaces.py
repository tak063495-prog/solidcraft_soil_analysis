"""Sealed 150-DOF results and one small binary64 eigenspace regression."""
import hashlib
import io
import json
from pathlib import Path
import sys
import unittest
import zipfile

import numpy as np

FOLDER = Path(__file__).resolve().parents[1] / "interim/2027-01-13-b01-highprec-eigenspace"


class B01Eigenspaces(unittest.TestCase):
    def test_sealed_eigenpairs(self):
        summary = json.loads((FOLDER / "measurements.json").read_bytes())
        index = json.loads((FOLDER / "archive-index.json").read_bytes())
        interruption, = index["supplemental_archives"]
        data = (FOLDER / interruption["name"]).read_bytes()
        self.assertEqual(len(data), interruption["size"])
        self.assertEqual(hashlib.sha256(data).hexdigest(), interruption["sha256"])
        archive = io.BytesIO()
        for part in index["parts"]:
            data = (FOLDER / part["name"]).read_bytes()
            self.assertEqual(len(data), part["size"])
            self.assertEqual(hashlib.sha256(data).hexdigest(), part["sha256"])
            archive.write(data)
        self.assertEqual(hashlib.sha256(archive.getbuffer()).hexdigest(), summary["archive_sha256"])
        self.assertEqual(summary["archive_sha256"], index["archive_sha256"])
        with zipfile.ZipFile(archive) as z:
            self.assertEqual(sorted(z.namelist()), sorted(index["entries"]))
            for name, member in sorted(index["entries"].items()):
                data = z.read(name)
                self.assertEqual(len(data), member["size"], name)
                self.assertEqual(hashlib.sha256(data).hexdigest(), member["sha256"], name)
            self.assertEqual(sum(p.startswith("raw/") for p in z.namelist()),32)
            self.assertEqual(sum(p.startswith("observations/") for p in z.namelist()),32)

    def test_saved_matrix_binary64(self):
        previous = sys.path[:]
        try:
            sys.path.insert(0, str(FOLDER))
            from b01eig_solve import load_inputs
            cases = {c["id"]: c for c in load_inputs()}
        finally:
            sys.path[:] = previous
        k = np.array(cases["regular--identity"]["stiffness"])
        # Decoupled tangential/drilling pair in the original flat patch.
        # The eigenspace projector is exact for this diagonal two-DOF slice.
        sub = k[np.ix_([3,5],[3,5])]
        lam,v = np.linalg.eigh(sub)
        np.testing.assert_array_equal(lam,[k[5,5],k[3,3]])
        np.testing.assert_array_equal(np.outer(v[:,0],v[:,0]),[[0,0],[0,1]])
        k = np.array(cases["interior-0.06--rx37"]["stiffness"])
        lam = np.linalg.eigvalsh((k+k.T)*0.5)
        summary = json.loads((FOLDER / "measurements.json").read_bytes())
        saved = next(s for s in summary["cases"] if s["id"]=="interior-0.06--rx37")
        self.assertEqual(np.count_nonzero(np.abs(lam)<=lam[-1]*1e-10),saved["mp_zero_count"])
        # The isolated positive mode remains between its two
        # adjacent saved modes. No componentwise eigenvector comparison.
        self.assertLess(float(saved["mp_lambda6"]),lam[6])
        self.assertLess(lam[6],float(saved["mp_lambda8"]))
