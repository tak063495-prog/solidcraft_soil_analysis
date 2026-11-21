"""Explicitly run here; outside existing tests/ and tools/ collection roots.

Numerical assertions qualify these small mathematical implementations only.
They are not inherited acceptance tolerances for P6-03 or the production solver.
"""
import copy
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np
from elastic_reference import elastic_matrix, lame, ring, foundation, solve, signs
from vibration_reference import exact, newmark, rod
from compare_outputs import compare, load, validate, FIELDS
from run_reference import dataset, record
from collect_legacy import collect


def sample():
    return dataset("example",{"E_Pa":30e6,"load":"explicit"},[
        record("p0",[0,0,0],"stress","zz","Pa",-12.,t) for t in [0.,.1]])


class ElasticTests(unittest.TestCase):
    def test_lame_tractions_and_radial_equilibrium(self):
        _,sra,_,_ = lame(1.)
        _,srb,_,_ = lame(2.)
        self.assertAlmostEqual(sra,-1e5)
        self.assertAlmostEqual(srb,0.)
        r,h = 1.4,1e-5
        _,sr,st,_ = lame(r)
        derivative = (lame(r+h)[1]-lame(r-h)[1])/(2*h)
        self.assertAlmostEqual(derivative+(sr-st)/r,0.,places=3)

    def test_affine_patch_interior_2d_3d(self):
        for dim in [2,3]:
            shape = (3,)*dim
            nodes = np.array(list(np.ndindex(shape)),float)
            lookup = {tuple(x):i for i,x in enumerate(nodes)}
            cells = []
            for low in np.ndindex((2,)*dim):
                cells.append([lookup[tuple(np.array(low)+(s+1)/2)] for s in signs(dim)])
            gradient = np.array([[.001,.0002,.0003],[-.0001,-.002,.0004],[.0005,-.0006,.0007]])[:dim,:dim]
            truth = nodes@gradient.T
            prescribed = {dim*i+a:truth[i,a] for i,x in enumerate(nodes) if np.any((x==0)|(x==2)) for a in range(dim)}
            result = solve(nodes,np.array(cells),np.zeros_like(nodes),prescribed)
            np.testing.assert_allclose(result["u"],truth,atol=1e-16,rtol=1e-12)
            strain = [gradient[0,0],gradient[1,1],gradient[2,2] if dim==3 else 0.,
                      gradient[0,1]+gradient[1,0],gradient[1,2]+gradient[2,1] if dim==3 else 0.,
                      gradient[0,2]+gradient[2,0] if dim==3 else 0.]
            stress = elastic_matrix(30e6,.25)@strain
            for g in result["gp"]:
                np.testing.assert_allclose(g["stress"],stress,atol=1e-9,rtol=1e-12)

    def test_plane_strain_and_generalized_3d(self):
        for eps in [0.,1e-4]:
            a,b = ring(ezz=eps),ring(dim=3,ezz=eps)
            self.assertAlmostEqual(a["displacement_relative_L2"],b["displacement_relative_L2"],places=12)
            self.assertAlmostEqual(a["stress_relative_L2"],b["stress_relative_L2"],places=12)
            self.assertLess(b["sigma_zz_relation_max_Pa"],1e-8)
            self.assertAlmostEqual(b["sigma_zz_exact_Pa"],16666.666666666668+30e6*eps)

    def test_ring_observed_error_orders(self):
        a,b = ring(4,8),ring(8,16)
        pu = np.log2(a["displacement_relative_L2"]/b["displacement_relative_L2"])
        ps = np.log2(a["stress_relative_L2"]/b["stress_relative_L2"])
        self.assertTrue(1.8 < pu < 2.2)
        self.assertTrue(.8 < ps < 1.2)

    def test_ring_balance_and_imposed_strain_work(self):
        for dim in [2,3]:
            r = ring(dim=dim,ezz=1e-4)
            self.assertLess(r["free_residual_relative"],1e-12)
            self.assertLess(np.linalg.norm(r["force_balance_N"]),1e-7)
            self.assertAlmostEqual(r["strain_energy_J"],r["boundary_work_half_J"],places=9)

    def test_finite_end_effect_distinguished_from_limit(self):
        a,b,plane = foundation(L=2.),foundation(L=8.),foundation(dim=2)
        limit = plane["profile"][0]["sigma_zz_Pa"]
        self.assertLess(abs(b["profile"][0]["sigma_zz_Pa"]-limit),abs(a["profile"][0]["sigma_zz_Pa"]-limit))
        self.assertGreater(abs(b["profile"][-1]["sigma_zz_Pa"]-limit),abs(b["profile"][0]["sigma_zz_Pa"]-limit))
        self.assertLess(b["free_residual_relative"],1e-11)
        self.assertAlmostEqual(b["strain_energy_J"],b["boundary_work_half_J"],places=8)
        self.assertEqual([p["x_m"] for p in b["profile"]],[0.]*9)

    def test_invalid_material_and_geometry_rejected(self):
        with self.assertRaises(ValueError):
            elastic_matrix(1.,.5)
        with self.assertRaises(ValueError):
            solve(np.array([[0,0],[0,1],[1,1],[1,0.]]),np.array([[0,1,2,3]]),np.zeros((4,2)),{0:0.})


class VibrationTests(unittest.TestCase):
    def test_exact_derivatives_and_initial_state(self):
        w,z = 13.,.03
        u,v,a = exact(np.array([0.,.123]),w,z)
        self.assertAlmostEqual(u[0],.01)
        self.assertAlmostEqual(v[0],.02)
        h,t = 1e-6,.123
        up,vp,_ = exact(t+h,w,z)
        um,vm,_ = exact(t-h,w,z)
        self.assertAlmostEqual((up-um)/(2*h),v[1],places=9)
        self.assertAlmostEqual((vp-vm)/(2*h),a[1],places=8)

    def test_newmark_second_order_and_damped_energy(self):
        w = np.sqrt(200.)
        errors = []
        for per in [40,80,160]:
            r = newmark(np.array([[2.]]),np.array([[2*.03*w*2]]),np.array([[400.]]),
                        np.array([.01]),np.array([.02]),2*np.pi/w/per,5*per)
            truth = exact(r["time_s"],w)[0]
            errors.append(np.linalg.norm(r["u"][:,0]-truth))
            self.assertLess(r["energy_balance_relative"],1e-11)
            self.assertTrue(np.all(np.diff(r["dissipation_J"])>=0))
        # RMS normalization needed because number of samples doubles.
        p = np.log2((errors[-2]/np.sqrt(401))/(errors[-1]/np.sqrt(801)))
        self.assertTrue(1.9 < p < 2.1)

    def test_undamped_energy_and_rod_spatial_frequency(self):
        errors = []
        for n in [4,8,16]:
            M,K,mode,w = rod(n)
            errors.append(w/(250*np.pi)-1)
            r = newmark(M,np.zeros_like(M),K,.01*mode,.02*mode,2*np.pi/w/40,80)
            self.assertLess(r["energy_balance_relative"],1e-10)
        self.assertTrue(1.9 < np.log2(errors[-2]/errors[-1]) < 2.1)


class ExchangeTests(unittest.TestCase):
    def test_cli_returns_two_and_refuses_existing_output(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp)
            (p/"input.json").write_text(json.dumps(sample()))
            cmd = [sys.executable,str(Path(__file__).with_name("compare_outputs.py")),str(p/"input.json"),str(p/"input.json"),str(p/"out.json")]
            first = subprocess.run(cmd,capture_output=True,text=True)
            self.assertEqual(first.returncode,2,first.stderr)
            self.assertIsNone(json.loads((p/"out.json").read_text())["accepted"])
            before = (p/"out.json").read_bytes()
            self.assertEqual(subprocess.run(cmd,capture_output=True).returncode,1)
            self.assertEqual((p/"out.json").read_bytes(),before)

    def test_roundtrip_json_csv(self):
        with tempfile.TemporaryDirectory() as temp:
            d = sample()
            p = Path(temp)
            (p/"a.json").write_text(json.dumps(d))
            (p/"m.json").write_text(json.dumps({k:v for k,v in d.items() if k!="records"}))
            with (p/"a.csv").open("w",newline="") as f:
                writer = csv.DictWriter(f,fieldnames=FIELDS)
                writer.writeheader()
                writer.writerows(d["records"])
            r = compare(load(p/"a.json"),load(p/"a.csv",p/"m.json"))
            self.assertIsNone(r["accepted"])
            self.assertEqual(r["groups"]["stress/zz"]["max_abs_difference"],0.)

    def test_undefined_even_for_identical_values(self):
        self.assertIsNone(compare(sample(),sample())["accepted"])

    def test_difference_preserved_without_sign_correction(self):
        a,b = sample(),sample()
        b["records"][0]["value"] *= -1
        self.assertEqual(compare(a,b)["groups"]["stress/zz"]["max_abs_difference"],24.)

    def test_unit_refused(self):
        d = sample(); d["records"][0]["unit"]="kPa"
        with self.assertRaises(ValueError): validate(d)

    def test_sign_refused(self):
        d = sample(); d["conventions"]=dict(d["conventions"],stress_sign="compression_positive")
        with self.assertRaises(ValueError): validate(d)

    def test_basis_refused(self):
        d = sample(); d["conventions"]=dict(d["conventions"],stress_basis="effective")
        with self.assertRaises(ValueError): validate(d)

    def test_component_missing_refused(self):
        a,b = sample(),sample(); b["records"][0]["component"]="xx"
        with self.assertRaises(ValueError): compare(a,b)

    def test_missing_id_refused(self):
        a,b = sample(),sample(); b["records"].pop()
        with self.assertRaises(ValueError): compare(a,b)

    def test_missing_time_no_interpolation(self):
        a,b = sample(),sample(); b["records"][1]["time_s"] += .0001
        with self.assertRaises(ValueError): compare(a,b)

    def test_coordinate_refused(self):
        a,b = sample(),sample(); b["records"][0]["z_m"]=1.
        with self.assertRaises(ValueError): compare(a,b)

    def test_duplicate_record_refused(self):
        d = sample(); d["records"].append(copy.deepcopy(d["records"][0]))
        with self.assertRaises(ValueError): validate(d)

    def test_nan_and_bool_refused(self):
        for value in [float("nan"),float("inf"),True]:
            d = sample(); d["records"][0]["value"]=value
            with self.assertRaises(ValueError): validate(d)

    def test_incomplete_refused(self):
        d = sample(); d["complete"]=False
        with self.assertRaises(ValueError): validate(d)

    def test_condition_change_refused(self):
        a,b = sample(),sample(); b["case"]["physical_conditions"]["E_Pa"]=20e6
        with self.assertRaises(ValueError): compare(a,b)

    def test_json_duplicate_keys_and_csv_no_manifest(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp)/"a.json"; p.write_text('{"records":[],"records":[]}')
            with self.assertRaises(ValueError): load(p)
            with self.assertRaises(ValueError): load(p.with_suffix(".csv"))


class LegacyTests(unittest.TestCase):
    def test_original_three_indeterminate_rows_preserved(self):
        d = collect()
        self.assertEqual(len(d["rows"]),3)
        self.assertFalse(d["p6_03_complete"])
        self.assertTrue(all(r["accepted"] is None for r in d["rows"]))
        self.assertGreater(d["rows"][0]["difference_ratio"],.7)
        self.assertTrue(all(p is None for p in d["rows"][2]["original_orders"]))

    def test_static_monitor_does_not_define_dynamic_time(self):
        d = collect()["rows"][1]
        self.assertLess(np.prod(d["signed_changes_m"]),0.)
        self.assertEqual(d["conditions"]["physical_time_s"],"未定義")
        self.assertTrue(d["period_s"].startswith("未定義"))

    def test_legacy_tamper_rejected_without_solver_fallback(self):
        here = Path(__file__).resolve().parent
        manifest = json.loads((here/"legacy-sources.json").read_text())
        first = manifest["sources"][0]["path"]
        with tempfile.TemporaryDirectory() as temp:
            changed = Path(temp)/first
            changed.parent.mkdir(parents=True)
            changed.write_text("changed")
            with self.assertRaisesRegex(ValueError,"missing/changed"):
                collect(Path(temp))


if __name__ == "__main__":
    unittest.main(verbosity=2)
