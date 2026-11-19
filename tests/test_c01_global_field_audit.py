"""Exercise the real global-step audit wrapper with a single synthetic Tet4."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

from test_c01_global_step_control import st as study

ROOT = Path(__file__).resolve().parents[1]
ISSUE_35 = ROOT / "interim/2027-01-11-issue35-c01-audit-fem-import"


class GlobalStepFieldAudit(unittest.TestCase):
    def audit_fixture(self, name, stress_error=0.0):
        mesh = {
            "nodes": [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
            "tetrahedra": [[0, 1, 2, 3]],
            "node_sets": {"fixed": [0]},
            "regions": [{"id": 0, "body": "soil"}],
            "region_ids": [0],
            "material_ids": [0],
        }
        material = dict(young_kpa=1000.0, poisson=0.25, cohesion_kpa=10.0,
                        friction_deg=30.0, dilation_deg=0.0, flow="nonassociated")
        initial = dict(effective_stress_compression_positive_kpa=[0.0] * 6, void_ratio=0.5)
        ground = dict(
            materials=[dict(id=0, mohr_coulomb=material)],
            initialization=dict(geostatic=dict(seed=dict(states=dict(soil=initial)))),
        )
        case = dict(
            stages=[dict(target_factor=0.0) for _ in study.STAGES[1:]],
            boundary=dict(constraints=[dict(node_set="fixed", value_m=0.0, direction=[1, 0, 0])]),
        )

        volume = study.fem.volume(mesh, 0)
        u = np.zeros((4, 3))
        joint_vector = np.zeros((1, 3))
        scalar = np.zeros(1)
        strain = np.zeros((1, 6))
        stress = strain.copy()
        stress[0, 0] = stress_error
        fields = dict(
            u=u, stress=stress, strain=strain, void=np.array([0.5]), pore=scalar,
            plastic_strain=strain, jump=joint_vector, traction=joint_vector, plastic_slip=joint_vector,
            accumulated=scalar, elastic_energy=scalar, plastic_dissipation=scalar,
            numerical_dissipation=scalar, endpoint_work=scalar, physical_pressure=scalar,
            water_load=u, reaction=u, gap=scalar,
        )
        load = u.ravel()
        # Supply equilibrium operators without running a contact or volume solve.
        op = SimpleNamespace(
            system=SimpleNamespace(
                c=np.zeros((1, 3, u.size)), normal=np.array([[0.0, 0.0, 1.0]]),
                kn=1.0, mu=0.5, weight=np.ones(1),
            ),
            water=SimpleNamespace(gamma=10.0), level=0.0,
            k=np.zeros((u.size, u.size)), seed=load, base=load, holding=load, boundary=load,
            ops=[(volume, np.arange(u.size))], interface=dict(pressure=scalar, load=u),
            reduce=lambda force: force,
        )
        with tempfile.TemporaryDirectory(prefix="issue 35 audit ") as tmp:
            folder = Path(tmp) / name
            writer = study.Writer(folder)
            for stage in study.STAGES:
                writer(dict(stage=stage, progress=1.0, factor=0.0), fields)
            writer.close()
            with patch.dict(sys.modules, {"study": study.prior}), patch.object(
                study.prior, "operators", return_value=op
            ):
                return study.audit(folder, mesh, ground, case, steps=1)

    def assert_minimal_history(self, name, tolerance):
        self.assertEqual(self.audit_fixture(name), dict(
            accepted_states=len(study.STAGES), maximum_reassembled_force_residual=0.0,
            existing_force_audit_tolerance=tolerance, maximum_pressure_change_kpa=0.0,
            maximum_pore_change_kpa=0.0, maximum_elastic_stress_compatibility_error_kpa=0.0,
            all_saved_states_audited=True, numerical_acceptance=False,
        ))

    def test_minimal_native_history_returns_residuals(self):
        self.assert_minimal_history("native", 1e-8)

    def test_minimal_independent_history_returns_residuals(self):
        self.assert_minimal_history("independent", 1e-10)

    def test_inconsistent_stress_is_rejected(self):
        with self.assertRaisesRegex(AssertionError, r"Max absolute difference among violations: 1"):
            self.audit_fixture("native", stress_error=1.0)


class GlobalStepAuditEvidence(unittest.TestCase):
    def test_old_and_repaired_audit_digests(self):
        record = json.loads((ISSUE_35 / "acceptance.json").read_text(encoding="utf-8"))
        path = study.HERE / "audit.py"
        name = path.relative_to(ROOT).as_posix()
        binding = record["audit"]
        self.assertEqual(binding["path"], name)

        git = ["git", "-C", str(ROOT)]
        revision = f"{binding['before_commit']}:{name}"
        old_blob = subprocess.check_output([*git, "rev-parse", revision], text=True).strip()
        before = subprocess.check_output([*git, "show", revision])
        before_sha = hashlib.sha256(before).hexdigest()
        self.assertEqual(old_blob, binding["before_blob"])
        self.assertEqual(before_sha, binding["before_sha256"])

        sealed = ROOT / "interim/2027-01-08-pc-failure-diagnosis"
        for manifest in ("final-code-binding.json", "validation-code-binding.json"):
            with self.subTest(manifest=manifest):
                historical = json.loads((sealed / manifest).read_text(encoding="utf-8"))
                self.assertEqual(historical[name], before_sha)

        after = path.read_bytes()
        after_blob = subprocess.check_output([*git, "hash-object", "--no-filters", str(path)],
                                            text=True).strip()
        self.assertEqual(after_blob, binding["after_blob"])
        self.assertEqual(hashlib.sha256(after).hexdigest(), binding["after_sha256"])
        self.assertNotEqual(before_sha, binding["after_sha256"])

    def test_documents_end_with_unresolved_status_on_lf_and_crlf(self):
        ending = "A-01〜A-04、B-01〜B-03、C-01〜C-04、NC-01、P6 は未達、passed=false。\n"
        with tempfile.TemporaryDirectory(prefix="issue 35 docs ") as tmp:
            for name in ("README.ja.md", "REPRODUCE.ja.md"):
                text = (ISSUE_35 / name).read_text(encoding="utf-8")
                self.assertTrue(text.endswith(ending), name)
                for newline in ("\n", "\r\n"):
                    with self.subTest(name=name, newline=repr(newline)):
                        copy = Path(tmp) / name
                        copy.write_bytes(text.replace("\n", newline).encode("utf-8"))
                        self.assertTrue(copy.read_text(encoding="utf-8").endswith(ending))
