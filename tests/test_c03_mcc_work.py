"""MCC work decomposition, sealed finite comparison, and exact-time alignment."""
import csv
import hashlib
import json
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FOLDER = ROOT / "interim/2026-11-23-c03-mcc-plastic-work"
sys.path.insert(0, str(FOLDER))
from compare_work import audit_force, summarize
from independent_work import components


class MccWork(unittest.TestCase):
    def test_pressure_and_deviatoric_exact_cases(self):
        stress = [[100,100,100,0,0,0], [180,60,60,0,0,0],
                  [120,120,60,60,0,0], [60,120,120,0,60,0]]
        dep = np.array([[.5,.5,.5,0,0,0], [4,-2,-2,0,0,0],
                        [1,1,-2,6,0,0], [-2,1,1,0,6,0]])/512
        terms = components(stress, dep)
        np.testing.assert_array_equal(terms, [[300/1024,300/1024,0]]+[[.9375,0,.9375]]*3)
        with self.assertRaisesRegex(ValueError, "accepted history"):
            summarize([], 1e-10)

    def test_records_and_split_scopes(self):
        folder = FOLDER / "results"
        v=json.loads((folder/"validation.json").read_bytes())
        for key,name in (("increments_sha256","increments.csv"),("summary_sha256","summary.json")):
            self.assertEqual(hashlib.sha256((folder/name).read_bytes()).hexdigest(),v[key])
        with (folder/"increments.csv").open(encoding="utf-8",newline="") as f:
            rows=list(csv.DictReader(f))
        self.assertEqual(len(rows),6*224)
        self.assertEqual(len({(r['case'],r['solver'],r['stage'],r['progress']) for r in rows}),len(rows))
        self.assertEqual(len({r['id'] for r in v['scopes']}),5)
        self.assertEqual([r['status'] for r in v['scopes']],['達成']*3+['未達']*2)
        self.assertIsNone(v['work_tolerance'])
        self.assertEqual(audit_force(folder),json.loads((folder/"force-audit.json").read_bytes()))
