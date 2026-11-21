"""Explicit offset joint contracts; each comparison is tied to its physical point."""

import copy
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "interim/2026-11-21-impl-t6-offset-kink-joint"))
from t6_patch import patch
from t6_reference import five, solve
from t6_attached import attached

CLI = (
    ROOT / "target/release" / ("soil-fem-3d.exe" if os.name == "nt" else "soil-fem-3d")
)


class OffsetKink(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name)
        self.serial = 0

    def execute(self, case, *args, fail=None):
        self.serial += 1
        source = self.path / "case.json"
        source.write_text(json.dumps(case), encoding="utf-8")
        out = self.path / f"run-{self.serial}.json"
        cmd = [CLI, "offset-kink-joint", source, out, *args]
        p = subprocess.run(cmd, capture_output=True, text=True)
        if fail:
            self.assertNotEqual(p.returncode, 0)
            self.assertIn(fail, p.stderr)
            self.assertFalse(out.exists())
            return None
        self.assertEqual(p.returncode, 0, p.stderr)
        return out

    def frame(self, case):
        return json.loads(self.execute(case).read_bytes())["frames"][-1]

    def test_signed_arms_same_physical_points_all_rigid_modes(self):
        c = patch(1)
        shell = c["shell"]
        shell["steps"][1]["free_dofs"] = [
            6 * n + k for n in c["joint"]["wings"][1]["edge_nodes"] for k in range(6)
        ]
        for k in range(6):
            v = np.eye(6)[k] * 0.001
            q = np.array(
                [np.r_[v[:3] + np.cross(v[3:], x), v[3:]] for x in shell["nodes"]]
            )
            shell["steps"][1]["prescribed_configuration"] = q.ravel().tolist()
            c["joint"]["loads"] = {}
            f = self.frame(c)
            self.assertLess(np.linalg.norm(f["reactions"]), 1e-9)
            for s in f["active"]["joint"]["stations"]:
                for name in ["joint", "attachment", "middle"]:
                    np.testing.assert_allclose(
                        s["u_" + name + "_m"],
                        v[:3] + np.cross(v[3:], s[name + "_m"]),
                        atol=1e-15,
                    )

    def test_released_and_spring_axial_moments(self):
        for law in ["rigid", "hinge", "spring"]:
            c = patch(1, law)
            f = self.frame(c)
            s = f["active"]["joint"]
            left, right = np.array_split(s["stations"], 2)
            angles = np.array(
                [
                    np.dot(
                        a["line_axis"],
                        np.subtract(b["rotation_rad"], a["rotation_rad"]),
                    )
                    for a, b in zip(left, right)
                ]
            )
            moment = np.array(
                [np.dot(b["line_axis"], b["moment_at_joint_kn_m"]) for b in right]
            )
            mass = np.array([[4, 2, -1], [2, 16, 2], [-1, 2, 4]]) / 30
            if law == "rigid":
                np.testing.assert_allclose(angles, 0.0, atol=1e-14)
            elif law == "hinge":
                np.testing.assert_allclose(moment, 0.0, atol=1e-10)
            else:
                np.testing.assert_allclose(moment, -25 * mass @ angles, atol=1e-10)
                self.assertAlmostEqual(
                    s["spring_energy_kn_m"], 25 * angles @ mass @ angles / 2, 13
                )

            c["joint"]["line_m"].reverse()
            for wing in c["joint"]["wings"]:
                wing["edge_nodes"].reverse()
            c["joint"]["loads"][2][0]["interval_m"] = [0.0, 0.75]
            reverse = self.frame(c)
            np.testing.assert_allclose(
                reverse["configuration"], f["configuration"], atol=1e-12
            )
            np.testing.assert_allclose(
                reverse["active"]["joint"]["relative_angle_increment_rad"],
                -angles[::-1],
                atol=1e-12,
            )

    def test_projected_five_dof_and_artificial_work(self):
        values = []
        for alpha in [1e-5, 1e-4, 1e-3]:
            c = patch(1, "rigid", alpha)
            frame = self.frame(c)
            ref = five(c)
            artificial = sum(
                e["drilling_energy_kn_m"] for e in frame["active"]["elements"]
            )
            physical = sum(
                e["physical_energy_kn_m"] for e in frame["active"]["elements"]
            )
            values.append(artificial)
            self.assertGreater(physical, 0.0)
            self.assertAlmostEqual(
                physical + artificial, frame["cumulative_applied_work_kn_m"], 12
            )
            self.assertGreater(ref["physical"], 0.0)
            self.assertEqual(ref["nullity"], 0)
        self.assertLess(values[0], values[-1])

    def test_one_sided_sections_and_line_densities(self):
        c = patch(1)
        frame = self.frame(c)
        data = frame["active"]["joint"]
        self.assertEqual({s["group"] for s in data["cuts"]}, {"wing-0", "wing-1"})
        for row in data["cuts"]:
            sec = row["section"]
            e = next(
                f["offset_m"] for f in c["shell"]["faces"] if f["id"] == row["face"]
            )
            np.testing.assert_allclose(
                sec["m_attachment_kn"],
                np.array(sec["m_kn"]) - e * np.array(sec["n_kn_m"]),
                atol=1e-13,
            )
            axes = np.array(sec["mid_axes"])
            np.testing.assert_allclose(axes @ axes.T, np.eye(3), atol=1e-14)
        for wing in np.array_split(data["stations"], 2):
            for density, total in [
                ("force_density_kn_m", "force_kn"),
                ("moment_density_kn_m_per_m", "moment_at_joint_kn_m"),
            ]:
                np.testing.assert_allclose(
                    np.array([1, 4, 1]) / 6 @ np.array([s[density] for s in wing]),
                    np.sum([s[total] for s in wing], axis=0),
                    atol=1e-13,
                )

    def test_surface_body_and_partial_line_wrench_once(self):
        c = patch(1)
        shell = c["shell"]
        shell["section"]["unit_weight_kn_m3"] = 2.0
        f = shell["faces"][0]
        z = 0.04
        shell["face_loads"] = {
            2: [
                dict(
                    face=f["id"], z_m=z, pressure_kpa=3.0, traction_kpa=[1.0, 0.0, 0.0]
                )
            ]
        }
        out = self.frame(c)
        xyz = np.array(shell["nodes"])
        load = np.array(out["applied"]).reshape(-1, 6)
        force = load[:, :3].sum(axis=0)
        moment = (np.cross(xyz, load[:, :3]) + load[:, 3:]).sum(axis=0)
        F = np.zeros(3)
        M = np.zeros(3)
        for face in shell["faces"]:
            a, b, d = xyz[face["nodes"][:3]]
            cross = np.cross(b - a, d - a)
            area = np.linalg.norm(cross) / 2
            normal = cross / (2 * area)
            mid = (a + b + d) / 3 - face["offset_m"] * normal
            weight = np.array([0.0, 0.0, -2.0 * 0.08 * area])
            F += weight
            M += np.cross(mid, weight)
            if face == f:
                traction = area * (np.array([1.0, 0.0, 0.0]) - 3 * normal)
                F += traction
                M += np.cross(mid + z * normal, traction)
        l = c["joint"]["loads"][2][0]
        line = 0.75 * np.array(l["force_kn_m"])
        F += line
        M += np.cross([0.625, 0, 0], line) + 0.75 * np.array(l["couple_kn_m_per_m"])
        np.testing.assert_allclose(force, F, atol=1e-13)
        np.testing.assert_allclose(moment, M, atol=1e-13)
        c["joint"]["loads"][2].append(copy.deepcopy(l))
        self.execute(c, fail="duplicate line load ID")

        # Use the six independent load bases with arbitrary nodal virtual
        # motion. Evaluate the displacement at the load itself for the work.
        c = patch(1)
        line = c["joint"]["loads"][2][0]
        line["interval_m"] = [0.0, 1.0]
        nodes = c["joint"]["wings"][0]["edge_nodes"]
        xyz = np.asarray(c["shell"]["nodes"])
        virtual = np.random.default_rng(614).normal(size=(len(xyz), 6))
        arm = xyz[nodes] - c["joint"]["line_m"]
        motion = virtual[nodes].copy()
        motion[:, :3] -= np.cross(motion[:, 3:], arm)
        mean = np.array([1.0, 4.0, 1.0]) @ motion / 6
        for base in np.eye(6):
            line["force_kn_m"] = base[:3].tolist()
            line["couple_kn_m_per_m"] = base[3:].tolist()
            applied = self.frame(c)["applied"]
            self.assertAlmostEqual(virtual.ravel() @ applied, mean @ base, 13)

    def test_independent_multiplier_and_coordinate_covariance(self):
        axes = Rotation.from_rotvec([0.27, -0.18, 0.43]).as_matrix()
        for angle in [0.0, 90.0, 135.0]:
            original = None
            for basis in [np.eye(3), axes]:
                c = patch(1, "hinge", angle=angle, axes=basis)
                f = self.frame(c)
                r = solve(c)
                np.testing.assert_allclose(
                    f["configuration"], r["q"], atol=1e-10, rtol=1e-8
                )
                np.testing.assert_allclose(f["reactions"], r["support"], atol=1e-9)
                for s in f["active"]["joint"]["stations"]:
                    j = 6 * s["node"]
                    force = r["transfer"][j : j + 3]
                    arm = np.subtract(s["attachment_m"], s["joint_m"])
                    np.testing.assert_allclose(s["force_kn"], force, atol=1e-9)
                    np.testing.assert_allclose(
                        s["moment_at_joint_kn_m"],
                        r["transfer"][j + 3 : j + 6] + np.cross(arm, force),
                        atol=1e-9,
                    )
                q = np.array(f["configuration"]).reshape(-1, 3)
                if original is None:
                    original = q
                else:
                    np.testing.assert_allclose(q, original @ axes.T, atol=1e-10)

    def test_generation_freeze_soil_tie_and_process_resume(self):
        c, nodes = attached()
        full = self.execute(c)
        frames = json.loads(full.read_bytes())["frames"]
        self.assertEqual(frames[3]["active"], frames[7]["archived"][0])
        self.assertEqual(frames[7]["active"]["generation"], 2)
        for f in frames:
            np.testing.assert_array_equal(
                np.array(f["configuration"][: 6 * len(nodes)]).reshape(-1, 6)[:, :3],
                np.array(f["soil_translations_m"])[nodes],
            )
        for stop in [2, 3, 4, 6]:
            saved = self.execute(c, "--until", str(stop))
            self.assertEqual(
                full.read_bytes(), self.execute(c, "--resume", saved).read_bytes()
            )
        bad = json.loads(saved.read_bytes())
        bad["frames"][-1]["archived"][0]["joint"]["stations"][0]["joint_m"][2] += 0.01
        saved.write_text(json.dumps(bad), encoding="utf-8")
        self.execute(c, "--resume", saved, fail="replay mismatch")

    def test_input_rank_trim_and_edge_ownership(self):
        c = patch(1)
        c["shell"]["steps"][1]["free_dofs"] = []
        self.execute(c, fail="constraint rank")
        c = patch(1)
        c["joint"]["separators"][0]["point_m"] = [0.0, 1.0, 0.0]
        self.execute(c, fail="uncertified trimmed")
        c = patch(1)
        c["joint"]["wings"][0]["edge_nodes"].reverse()
        self.execute(c, fail="nonconforming")

        # Two exterior faces of Tet10 cell 10 share soil edge 8--41--9.
        # Give each shell independent IDs. When J is displaced from their
        # common attachment line, the shared translations would lock the hinge.
        case, _ = attached()
        shell = case["shell"]
        soil = shell["soil"]
        mesh = soil["mesh"]
        owners = {"top": [8, 9, 11, 41, 43, 40], "side": [4, 9, 8, 42, 41, 39]}
        shell["nodes"] = []
        shell["faces"] = []
        joint = case["joint"]
        joint["wings"] = []
        for group, ids in owners.items():
            start = len(shell["nodes"])
            shell["nodes"].extend(mesh["nodes"][n] for n in ids)
            shell["faces"].append(
                dict(
                    id=group,
                    group=group,
                    nodes=list(range(start, start + 6)),
                    offset_m=-0.05,
                )
            )
            joint["wings"].append(
                dict(group=group, edge_nodes=[start + ids.index(n) for n in [8, 41, 9]])
            )
        soil["owners"] = {g: dict(cell=10, nodes=ids) for g, ids in owners.items()}
        soil["exposed_side"] = {g: 1 for g in owners}
        joint.update(
            rotation=dict(kind="hinge"),
            line_m=[[x, -0.02, 1.0] for x in [0.0, 0.5, 1.0]],
            loads={},
            separators=[
                dict(
                    negative="side",
                    positive="top",
                    normal=[0.0, 1.0, 0.0],
                    point_m=[0.0, 0.0, 1.0],
                )
            ],
        )
        shell["steps"] = []
        self.execute(case, fail="shared tie adds rank 1/3")

    def test_three_wing_sections_and_declaration_order(self):
        c = patch(1)
        shell = c["shell"]
        joint = c["joint"]
        count = len(shell["nodes"]) // 2
        shift = 2 * count
        shell["nodes"] += (
            np.asarray(shell["nodes"][:count]) + [0.0, 0.0, 1.5]
        ).tolist()
        for f in list(shell["faces"]):
            if f["group"] == "wing-0":
                shell["faces"].append(
                    dict(
                        f,
                        id="third-" + f["id"],
                        group="third",
                        nodes=[n + shift for n in f["nodes"]],
                        offset_m=0.0,
                    )
                )
        third = copy.deepcopy(joint["wings"][0])
        third.update(
            group="third",
            edge_nodes=[n + shift for n in third["edge_nodes"]],
            section=dict(shell["section"], thickness_m=0.12),
        )
        joint["wings"].append(third)
        joint["separators"].extend(
            [
                dict(
                    negative=g,
                    positive="third",
                    point_m=[0.0, 0.0, 1.3],
                    normal=[0.0, 0.0, 1.0],
                )
                for g in ["wing-0", "wing-1"]
            ]
        )
        step = shell["steps"][1]
        n = 6 * count
        step["free_dofs"] += [j + 6 * shift for j in step["free_dofs"] if j < n]
        for key in ["prescribed_configuration", "nodal_load"]:
            step[key].extend([0.0] * n)
        a = self.frame(c)
        joint["wings"].reverse()
        shell["faces"].reverse()
        b = self.frame(c)
        np.testing.assert_allclose(a["configuration"], b["configuration"], atol=1e-11)
        # Sign change of both the director and the offset keeps the physical middle point.
        for f in shell["faces"]:
            if f["group"] == "wing-0":
                f["nodes"] = [f["nodes"][i] for i in (0, 2, 1, 5, 4, 3)]
                f["offset_m"] = -f["offset_m"]
        rev = self.frame(c)
        np.testing.assert_allclose(a["configuration"], rev["configuration"], atol=1e-11)
