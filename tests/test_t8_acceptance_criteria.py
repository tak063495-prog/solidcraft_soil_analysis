from collections import Counter
import copy
import io
import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
import zipfile

import numpy as np

STUDY = Path(__file__).resolve().parents[1] / "interim/2026-11-22-impl-t8-acceptance-criteria"
sys.path.insert(0, str(STUDY))
import assess as assessment
import water
from branches import paired, support_reason
from refinement import richardson
from two_mesh import metrics
from observations import verify


class AcceptanceCriteria(unittest.TestCase):
    def estimate(self, values, widths=(1., .5, .25), **changes):
        conditions = dict(similar=True, complete=True, reference_converged=True,
                          other_axes_fixed=True, smooth=True)
        return richardson(values, widths, **(conditions | changes))

    def test_manufactured_second_order_scalar_two_successive_triples(self):
        for widths in ((1., .5, .25), (.5, .25, .125)):
            with self.subTest(widths=widths):
                result = self.estimate([2 + h*h for h in widths], widths)
                self.assertEqual(result["holds"], [])
                self.assertEqual(result["observed_order"], 2)
                self.assertEqual(result["extrapolated"], 2)
                self.assertAlmostEqual(result["gci_fine"], 1.25 * widths[-1]**2 / (2+widths[-1]**2))

    def test_hold_conditions_keep_extrapolation_empty(self):
        cases = [((1, 1.5, 1.25), {}, "reverse"), ((1, 2, 4), {}, "non-positive"),
                 ((1, 1, 1), {}, "roundoff"), ((1, 2, None), {}, "missing"),
                 ((1, 2), {}, "three levels"), ((1, 2, 2.5), {"similar":False}, "non-similar"),
                 ((1, 2, 2.5), {"widths":(1., .4, .25)}, "unequal"),
                 ((1, 2, 2.5), {"complete":False}, "incomplete"),
                 ((1, 2, 2.5), {"reference_converged":False}, "unresolved"),
                 ((1, 2, 2.5), {"other_axes_fixed":False}, "mixed"),
                 ((1, 2, 2.5), {"smooth":False}, "non-smooth")]
        for q, options, reason in cases:
            with self.subTest(reason=reason):
                result = self.estimate(q, **options)
                self.assertIn(reason, " ".join(result["holds"]))
                self.assertIsNone(result["extrapolated"])
                self.assertIsNone(result["gci_fine"])

    def test_zero_fine_observable_retains_dimensional_correction(self):
        result = self.estimate([3., 1., 0.])
        self.assertAlmostEqual(result["dimensional_correction"], 1)
        self.assertIsNone(result["gci_fine"])

    def test_branch_correspondence_precedes_field_difference(self):
        fields = dict(stage="open", generations=[1,1], event=["stick","slip"])
        point = dict(actual=.02, bracket=[.01,.03], matched_side=True, fields=fields)
        self.assertIsNone(paired(point, point, .02))
        for altered in (dict(matched_side=False), dict(bracket=[.03,.01]), dict(actual=.025),
                        dict(fields=fields|{"stage":"birth"}), dict(fields=fields|{"generations":[1,2]}),
                        dict(fields=fields|{"event":["stick","stick"]})):
            with self.subTest(altered=altered):
                self.assertIsNotNone(paired(point, point | altered, .02))

        mesh = dict(nodes=[[0, 0, 0]], node_sets={"top":[0]})
        case = dict(stages=[dict(name="shear",kind="boundary",target_factor=0.),
                            dict(name="open",kind="boundary",target_factor=1.)],
                    boundary={"constraints":[dict(node_set="top",direction=[1,0,0],value_m=1.)]})
        correct = dict(fields=fields | {"u":[.02,0.,0.]})
        wrong = dict(fields=fields | {"u":[.01,0.,0.]})
        self.assertIsNone(support_reason((correct, correct), case, mesh, .02))
        self.assertIsNotNone(support_reason((correct, wrong), case, mesh, .02))

    def test_integral_cancellation_preserves_point_error_and_location(self):
        x = np.array([[[1.,0.],[-1.,0.]]])
        r = metrics(x, np.zeros_like(x), np.array([.5,.5]), [("open",1)],
                    np.array([[0.,0.,0.],[1.,0.,0.]]))
        self.assertEqual(r["surface_integral"]["absolute_l2"], 0)
        self.assertEqual(r["area_weighted"]["absolute_l2"], 1)
        self.assertEqual(r["point"]["maximum_absolute"], 1)
        self.assertEqual(r["maximum_at"]["accepted"], ("open",1))


    def test_all_original_rows_and_verified_limited_scopes(self):
        self.assertGreater(verify(STUDY), 0)
        result = json.loads((STUDY/"acceptance.json").read_bytes())
        self.assertEqual(set(result["items"]), {f"A-0{i}" for i in range(1,5)} |
                         {f"B-0{i}" for i in range(1,4)} | {f"C-0{i}" for i in range(1,5)} | {"NC-01"})
        self.assertTrue(all(r["original"] == "未達" for r in result["items"].values()))
        self.assertFalse(result["p6_complete"])

    def test_finite_references_and_strict_category_inventory(self):
        with zipfile.ZipFile(STUDY/"observations.zip") as z:
            c1 = json.loads(z.read("c01-results.json"))
            c2 = json.loads(z.read("c02-extended.json"))
            c3 = json.loads(z.read("c03/results.json"))
            c4 = json.loads(z.read("c04-results.json"))
        self.assertLess(c1["directions"]["independent"]["plastic_dissipation"]["observed_order"], 0)
        self.assertIsNone(c1["gci_fine"])
        self.assertEqual(c2["common_points"], 160)
        self.assertEqual(len(c2["time_mesh_exchange"]), 12)
        self.assertIsNone(c2["observed_order"])
        sequence = [c3["sensitivity"][f"m0c{i}-native--m0c{i+1}-native"]["all_fields"]["plastic_slip"]["absolute_l2"]
                    for i in range(3)]
        self.assertGreater(sequence[1],sequence[0])
        self.assertFalse(c3["reference_assessed_as_exact"])
        self.assertEqual(c4["preserved_mismatch_counts"], [36,63,133])
        fine = next(r for r in c4["measured"] if (r["mechanical_steps"],r["water_steps"]) == (32,32))
        self.assertEqual(fine["strict_mismatches"], 133)
        self.assertEqual(len(fine["boundary"]), 133)
        self.assertEqual(sum(fine["confusion"].values()), 2496)
        self.assertEqual(sum(r["zero_saved_plastic_increment_both"] for r in fine["boundary"]), 133)
        self.assertFalse(c4["boundary_reclassification"])

    def test_admissible_projection_and_moving_water_plane(self):
        with zipfile.ZipFile(STUDY/"observations.zip") as z:
            projection = json.loads(z.read("a02/results.json"))
            placements = json.loads(z.read("t8-water/placement.json"))
            water = json.loads(z.read("t8-water/affine.json"))
        self.assertTrue(projection["projection/good"]["accepted"])
        for name in ("tension", "outside-circle", "balanced-unreproducible-tip"):
            self.assertFalse(projection["projection/"+name]["accepted"])
        self.assertEqual([r["level_m"] for r in placements], [.625,.75,.875])
        for row in placements:
            self.assertIsNone(row["whole_domain_budget"])
            for key, value in row["expected"].items():
                self.assertLess(abs(value-row["measured"][key]), 1e-12)
        self.assertEqual(len(water), 6)
        for row in water:
            self.assertTrue(all(v["passed"] for v in row["water"]["profiles"].values()))
            for point in row["traces"]:
                self.assertLess(abs(point["a_flow_error_m3_s"]), 1e-14)
                self.assertLess(abs(point["b_flow_error_m3_s"]), 1e-14)
                self.assertLess(point["max_nodal_force_error_kn"], 1e-10)
                self.assertEqual(point["expected_a_outflow_m3_s"] == 0, row["slope"] == 0)

    def test_affine_solver_output_requires_each_face_point_and_time(self):
        mesh = water.controls.read(water.controls.EX / "joint-water-column-mesh.json")
        frames = []
        for seconds in (.125, .25):
            traces = []
            forces = []
            for face, (a_nodes, b_nodes) in enumerate(zip(
                    mesh["face_sets"]["Seam:a"], mesh["face_sets"]["Seam:b"], strict=True)):
                owners = [next(i for i, tet in enumerate(mesh["tetrahedra"])
                               if set(nodes) <= set(tet)) for nodes in (a_nodes, b_nodes)]
                z = mesh["nodes"][a_nodes[0]][2]
                pressure = 10 * (2 - z)
                cell_z = sum(mesh["nodes"][n][2] for n in mesh["tetrahedra"][owners[0]]) / 4
                outward = 1 if z > cell_z else -1
                for point in range(3):
                    traces.append(dict(joint_id="Seam", face_index=face, point_index=point,
                                       time_s=seconds, elevation_m=z, applied_pressure_kpa=pressure,
                                       **{f"side_{side}": dict(face=face, owner=owner,
                                                              total_head_m=2., outward_flow_m3_s=0.)
                                          for side, owner in zip("ab", owners, strict=True)}))
                forces.append([[0., 0., -outward * pressure * .5 / 3]] * 3)
            frames.append(dict(phase="transient", state=dict(
                hydraulic=dict(time_s=seconds), fluid=dict(joint_traces=traces),
                mechanical=dict(joint_water=dict(face_nodal_force_a_kn=forces)))))

        def record(name, mesh, run, folder, budget):
            budget[name] = {}

        with patch.object(water.controls, "record", side_effect=record), \
             patch.object(water, "solve") as solve:
            solve.return_value = dict(frames=frames)
            result = water.affine_column(STUDY, None, 0., .125)
            self.assertEqual(len(result["traces"]), 12)

            for mutation in ("missing", "duplicate face", "duplicate point", "extra", "late"):
                with self.subTest(mutation=mutation):
                    invalid = copy.deepcopy(frames)
                    for frame in invalid:
                        traces = frame["state"]["fluid"]["joint_traces"]
                        if mutation == "missing":
                            del traces[3:]
                        elif mutation == "duplicate face":
                            traces[3:] = copy.deepcopy(traces[:3])
                        elif mutation == "duplicate point":
                            traces[1] = copy.deepcopy(traces[0])
                        elif mutation == "extra":
                            traces.append(copy.deepcopy(traces[0]))
                        else:
                            frame["state"]["hydraulic"]["time_s"] += .25
                            for trace in traces:
                                trace["time_s"] += .25
                    solve.return_value = dict(frames=invalid)
                    with self.assertRaises(ValueError):
                        water.affine_column(STUDY, None, 0., .125)

    def test_affine_assessment_requires_both_profiles_and_the_complete_matrix(self):
        work = Path("sealed-study")
        archived_read = assessment.read
        with zipfile.ZipFile(STUDY / "observations.zip") as archive:
            original = json.loads(archive.read("t8-water/affine.json"))
            for row in original:
                counts = Counter()
                for trace in row["traces"]:
                    key = trace["time_s"], trace["face"]
                    trace["point_index"] = counts[key]
                    counts[key] += 1

            def read(path):
                if path == work / "t8-water/affine.json":
                    return rows
                if path.is_relative_to(work):
                    return json.loads(archive.read(path.relative_to(work).as_posix()))
                return archived_read(path)

            rows = original
            branches = archived_read(STUDY / "comparisons.json")["branches"]
            with patch.object(assessment, "read", side_effect=read), \
                 patch.object(assessment, "evaluate", return_value=branches):
                result, _ = assessment.assess(work)
                self.assertEqual(result["items"]["A-04"]["limited"][0]["status"], "限定受入")
                for mutation in ("no profiles", "one profile", "no traces", "missing face",
                                 "duplicate face", "duplicate point", "extra trace", "duplicate path"):
                    with self.subTest(mutation=mutation):
                        rows = copy.deepcopy(original)
                        row = rows[0]
                        if mutation == "no profiles":
                            row["water"]["profiles"].clear()
                        elif mutation == "one profile":
                            row["water"]["profiles"].pop(next(iter(row["water"]["profiles"])))
                        elif mutation == "no traces":
                            row["traces"].clear()
                        elif mutation == "missing face":
                            row["traces"] = [trace for trace in row["traces"] if trace["face"] == 0]
                        elif mutation == "duplicate face":
                            row["traces"] = [trace | {"face": 0} for trace in row["traces"]]
                        elif mutation == "duplicate point":
                            row["traces"][1] = copy.deepcopy(row["traces"][0])
                        elif mutation == "extra trace":
                            row["traces"].append(row["traces"][0])
                        else:
                            rows[1] = copy.deepcopy(row)
                        result, _ = assessment.assess(work)
                        self.assertEqual(result["items"]["A-04"]["limited"][0]["status"], "未達")
                        self.assertEqual(result["items"]["A-04"]["original"], "未達")

    def test_shell_references_and_reassembly_failures_remain_visible(self):
        with zipfile.ZipFile(STUDY/"observations.zip") as z:
            standard = json.loads(z.read("b01/standard/measurements.json"))
            joints = json.loads(z.read("b03/controls.json"))
        fine = [r for r in standard if r["n"] == 8]
        self.assertEqual(len(fine), 6)
        self.assertEqual(sum("displacement" not in r for r in fine), 2)
        self.assertTrue(all(not r["passed"] for r in standard))
        for r in joints:
            self.assertEqual(r["modes"]["free"]["nullity"], {"rigid":6,"hinge":7,"spring":6}[r["law"]])
        table = json.loads((STUDY/"comparisons.json").read_bytes())
        self.assertEqual(sum(r["native"]["max_relative"]>1e-10 for r in table["thin_original"].values()), 4)
        self.assertIsNone(table["nc01"]["physical_collapse_load"])
        self.assertEqual(table["nc01"]["physical_classification"], "D")

    def test_status_suffix_under_lf_and_crlf(self):
        ending = "A-01～C-04、NC-01、P6 の元項目は未達。\n"
        for name in ("README.ja.md","REPRODUCE.ja.md"):
            data = (STUDY/name).read_bytes().replace(b"\r\n", b"\n")
            for newline in (b"\n", b"\r\n"):
                with io.TextIOWrapper(io.BytesIO(data.replace(b"\n",newline)), encoding="utf-8") as f:
                    self.assertTrue(f.read().endswith(ending), name)


if __name__ == "__main__":
    unittest.main()
