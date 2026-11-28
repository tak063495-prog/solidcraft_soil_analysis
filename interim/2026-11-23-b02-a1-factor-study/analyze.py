"""Reproduce the one A1 investigation from the repository root, on the same 18 x 9 mesh."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent / "2026-11-23-b01-b02-refined-cases"))
from refined import nodal_values, projected
from t10_compare import relative, shape

MODES = ("baseline", "shear_direct", "membrane_direct", "gauss_edge", "full_clamp",
         "duffy_shear", "single_section", "si_units")
REFERENCE = np.array([-500., 5000.])


def build(directory):
    output = subprocess.check_output(["cargo", "build", "--locked", "--release", "--example",
                                      "t10_layered_cantilever", "--message-format=json"], cwd=ROOT)
    libraries = {}
    for line in output.splitlines():
        item = json.loads(line)
        if item["reason"] == "compiler-artifact":
            for name in item["filenames"]:
                if Path(name).suffix == ".rlib":
                    libraries[item["target"]["name"]] = Path(name)
    deps = libraries["soil_fem_3d"].parent
    exe = directory / ("probe.exe" if sys.platform == "win32" else "probe")
    args = ["rustc", "--edition=2024", "-Dwarnings", "-Copt-level=3", "-L", f"dependency={deps}"]
    for name in ("soil_fem_3d", "serde_json", "nalgebra"):
        args.extend(["--extern", f"{name}={libraries[name]}"])
    subprocess.run([*args, str(HERE / "probe.rs"), "-o", str(exe)], check=True)
    for mode in MODES:
        with (directory / f"{mode}.json").open("xb") as output:
            subprocess.run([str(exe), mode], stdout=output, check=True)


def corner_edge(vector, offset):
    # On the lower edge x in [0,h], tied tangential shear is linear through the
    # two Gauss samples. With nu=0 about z=0, Mxx uses the two axial gradients.
    h = 10/18
    u, w, beta = np.asarray(vector).reshape(6, 6)[[0, 3, 1]][:, [0, 2, 4]].T
    def slope(v, t):
        return np.dot(v, [4*t-3, 4-8*t, 4*t-1]) / h
    def gamma(t):
        return slope(w, t) + np.dot(beta, [(1-t)*(1-2*t), 4*t*(1-t), t*(2*t-1)])
    a, b = (1 + np.array([-1, 1])/np.sqrt(3))/2
    q = 1000 * 5/6 * 1e8 * .4 * (b*gamma(a) - a*gamma(b))/(b-a)
    m = 1000 * 2e8 * ((.4**3/12 + offset**2*.4)*slope(beta, 0) - offset*.4*slope(u, 0))
    return [q, m]


def read_layer(row, mode, unit):
    # A single equivalent layer carries the whole section; halving applies only
    # to this symmetric pair under transverse bending.
    factor = (0.5 if mode == "single_section" else 1.) / unit
    report = {"offset_m": row["offset"]}
    samples = row["solved"]
    for field in ("solved", "swap_readout", "exact_translations", "exact_rotations", "exact_both"):
        if field in row:
            report[field] = (factor * nodal_values(row[field])[0]).tolist()
    if mode == "baseline":
        report["direct_edge_endpoint"] = corner_edge(row["solution_at_nodes"], row["offset"])
        # Same DOFs/sites/axes, changing only the plane for the moment.
        mid = [dict(p, m_attachment_kn=p["m_kn"]) for p in samples]
        report["layer_midplane"] = nodal_values(mid)[0].tolist()
        nearest = min(range(len(samples)), key=lambda i: np.linalg.norm(samples[i]["rs"]))
        x = (10/18)*samples[nearest]["rs"][0]
        report["nearest"] = {"x_m": x, "actual": projected([samples[nearest]])[0].tolist(),
                             "analytical_at_x": [-500., 500.*(10.-x)]}
        exact = nodal_values(row["exact_both"])
        along = (1 + np.array([-1, 1])/np.sqrt(3))/2
        report["exact_nodes_tying_QX"] = {
            "lower_edge": [float(shape(t, 0)[0] @ exact[:, 0]) for t in along],
            "centroid": float(shape(1/3, 1/3)[0] @ exact[:, 0]),
            "centroid_minus_continuum_closed_form": 1000*(5/6)*1e8*.4*(10/18)**2 /
                                                  (36*(2e8*.8**3/12))}
        report["exact_nodes_gauss_QX"] = projected(row["exact_both"])[:, 0].tolist()
        q, m = report["solved"]
        report["qx"] = relative(q, -500., 0.001)
        report["mxx"] = relative(m, 5000., 0.001)
    return report


def summarize(directory):
    modes = {}
    for mode in MODES:
        raw = json.loads((directory / f"{mode}.json").read_bytes())
        if raw["mode"] != mode:
            raise ValueError("mode mismatch")
        layers = [read_layer(layer, mode, raw["force_units_per_kn"]) for layer in raw["layers"]]
        modes[mode] = {"layers": layers, "A1_per_layer": np.mean([v["solved"] for v in layers],axis=0).tolist(),
                       "tip_dz_m": raw["tip_dz_m"], "reaction_kN_kNm": raw["root_reactions_kn_knm"]}
        if mode == "baseline":
            modes[mode]["free_component_max_residual"] = raw["free_component_max_residual"]
    base = np.array(modes["baseline"]["A1_per_layer"])
    for row in modes.values():
        value = np.array(row["A1_per_layer"])
        row["change_from_baseline"] = (value-base).tolist()
        row["signed_percent"] = (100*(value-REFERENCE)/abs(REFERENCE)).tolist()
    old = json.loads((HERE.parent / "2026-11-23-b01-b02-refined-cases/results/b02-a1.json").read_bytes())
    previous = next(r for r in old["grids"] if r["grid"] == "18x9")["layers"]
    discrepancies = []
    for measured, stored in zip(modes["baseline"]["layers"], previous, strict=True):
        if measured["offset_m"] != stored["offset_m"]:
            raise ValueError("layer mismatch")
        discrepancies.append((np.array(measured["solved"]) -
                              [stored["qx_n_m"]["measured"], stored["mxx_nm_m"]["measured"]]).tolist())
    return {"case": "B-02 A1 final single-factor investigation", "grid": [18,9],
            "reference_QX_MXX": REFERENCE.tolist(), "relative_width": 0.001,
            "baseline_delta_against_published_layers": discrepancies, "modes": modes,
            "a1_passed": all(g[k]["passed"] for g in modes["baseline"]["layers"] for k in ("qx","mxx"))}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("raw_directory", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--read-existing", action="store_true")
    args = parser.parse_args()
    if not args.read_existing:
        args.raw_directory.mkdir(parents=True, exist_ok=False)
        build(args.raw_directory.resolve())
    with args.output.open("x", encoding="utf-8", newline="\n") as f:
        json.dump(summarize(args.raw_directory), f, indent=2, ensure_ascii=False, allow_nan=False)
        f.write("\n")
