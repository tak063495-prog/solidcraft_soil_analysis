"""Run structural series and retain like-point comparisons against a KKT solve."""

import argparse
import json
from pathlib import Path
import subprocess
import time

import numpy as np
from t6_patch import patch
from t6_reference import solve, five, modes


def run(cli, case, folder):
    folder.mkdir(parents=True, exist_ok=False)
    source = folder / "case.json"
    source.write_text(json.dumps(case, allow_nan=False), encoding="utf-8")
    dest = folder / "run.json"
    subprocess.run([cli, "offset-kink-joint", source, dest], check=True)
    return json.loads(dest.read_bytes())


def symmetric(v):
    return np.array([[v[0], v[2]], [v[2], v[1]]])


def norms(d, w, ref):
    d = np.asarray(d)
    ref = np.asarray(ref)
    d = d.reshape(len(w), -1)
    ref = ref.reshape(d.shape)
    w = np.asarray(w)
    l2 = float(np.sqrt(w @ (d * d).sum(axis=1)))
    den = float(np.sqrt(w @ (ref * ref).sum(axis=1)))
    return dict(
        max_abs=float(abs(d).max()),
        signed_min=float(d.min()),
        signed_max=float(d.max()),
        weighted_l2=l2,
        reference_l2=den,
        relative_l2=l2 / den if den else None,
    )


def compare(case, run_data):
    r = solve(case)
    frame = run_data["frames"][-1]
    gen = frame["active"]
    q = np.array(frame["configuration"])
    delta = q - r["q"]
    stations = gen["joint"]["stations"]
    t = np.array(stations[0]["line_axis"])
    xyz = np.asarray(case["shell"]["nodes"])
    jf = []
    ref_forces = []
    ju = []
    ref_u = []
    lw = []
    for row in stations:
        i = row["node"]
        arm = xyz[i] - row["joint_m"]
        f = r["transfer"][6 * i : 6 * i + 3]
        m = r["transfer"][6 * i + 3 : 6 * i + 6] + np.cross(arm, f)
        u = r["q"][6 * i : 6 * i + 3] - np.cross(r["q"][6 * i + 3 : 6 * i + 6], arm)
        jf.append(np.r_[row["force_kn"], row["moment_at_joint_kn_m"]] - np.r_[f, m])
        ref_forces.append(np.r_[f, m])
        ju.append(row["u_joint_m"] - u)
        ref_u.append(u)
    for wing in case["joint"]["wings"]:
        arc = [s["arc_m"] for s in stations if s["group"] == wing["group"]]
        w = np.zeros(len(arc))
        for i in range(0, len(w) - 2, 2):
            w[i : i + 3] += (arc[i + 2] - arc[i]) * np.array([1, 4, 1]) / 6
        lw.extend(w)
    diffs = {
        k: [] for k in ["u", "rotation", "N", "M", "Q", "stress", "transverse_stress"]
    }
    refs = {k: [] for k in diffs}
    area = []
    vol = []
    split = np.zeros(3)

    def save(k, val, ref):
        diffs[k].append(np.asarray(val) - ref)
        refs[k].append(ref)

    from t6_reference import shape

    for element, obs in zip(r["elements"], gen["elements"], strict=True):
        u = r["q"][element.dofs]
        du = delta[element.dofs].reshape(6, 6)
        for j, (x, y, wt, B, D, S) in enumerate(element.samples):
            sec = obs["resultants"][j]
            L = np.array(sec["mid_axes"])[:2] @ element.basis.T
            mid = (B - element.e * D) @ u
            curv = D @ u
            shear = S @ u
            sig = element.elastic @ mid
            grad = element.elastic @ curv
            tau = element.g * 5 / 6 * shear

            def components(v):
                tensor = L @ symmetric(v) @ L.T
                return tensor[[0, 1, 0], [0, 1, 1]]

            save("N", sec["n_kn_m"], components(element.t * sig))
            save("M", sec["m_kn"], components(element.t**3 / 12 * grad))
            save("Q", sec["q_kn_m"], L @ tau * element.t)
            h = shape(x, y)[0]
            value = h @ u.reshape(6, 6)
            save("u", value[:3] + (h @ du)[:3], value[:3])
            save("rotation", value[3:] + (h @ du)[3:], value[3:])
            split += (
                wt
                / 2
                * np.array(
                    [
                        element.t * mid @ sig,
                        element.t**3 / 12 * curv @ grad,
                        element.t * shear @ tau,
                    ]
                )
            )
            area.append(wt)
            for field in obs["fields"][j::6]:
                save(
                    "stress", field["stress_kpa"], components(sig + field["z_m"] * grad)
                )
                save("transverse_stress", field["shear_kpa"], L @ tau)
                vol.append(field["weight"])
    cuts = {}
    for row in gen["joint"]["cuts"]:
        idx = next(
            i for i, f in enumerate(case["shell"]["faces"]) if f["id"] == row["face"]
        )
        e = r["elements"][idx]
        u = r["q"][e.dofs]
        sec = row["section"]
        B, D, S, _ = e.strains(*sec["rs"])
        axes = np.array(sec["mid_axes"])[:2] @ e.basis.T

        def tensor_result(v):
            mat = axes @ symmetric(v) @ axes.T
            return mat[[0, 1, 0], [0, 1, 1]]

        ref_n = tensor_result(e.t * e.elastic @ ((B - e.e * D) @ u))
        ref_m = tensor_result(e.t**3 / 12 * e.elastic @ (D @ u))
        ref_q = e.t * e.g * 5 / 6 * axes @ (S @ u)
        cuts[row["face"]] = dict(
            arc_m=row["arc_m"],
            group=row["group"],
            outward=row["outward_sign"],
            N=sec["n_kn_m"],
            M=sec["m_kn"],
            Q=sec["q_kn_m"],
            N_difference=(sec["n_kn_m"] - ref_n).tolist(),
            M_difference=(sec["m_kn"] - ref_m).tolist(),
            Q_difference=(sec["q_kn_m"] - ref_q).tolist(),
        )
    return dict(
        cuts=cuts,
        dofs=len(q),
        elements=len(gen["elements"]),
        constraint_residual=r["constraint"],
        nodal_max_abs=float(abs(delta).max()),
        applied_max_abs=float(max(abs(frame["applied"] - r["load"]))),
        support_max_abs=float(max(abs(frame["reactions"] - r["support"]))),
        joint_displacement=norms(ju, lw, ref_u),
        joint_nodal_wrench=norms(jf, np.ones(len(jf)), ref_forces),
        fields={
            k: norms(diffs[k], vol if "stress" in k else area, refs[k]) for k in diffs
        },
        native_energies=[
            sum(e[k] for e in gen["elements"])
            for k in ("physical_energy_kn_m", "drilling_energy_kn_m")
        ]
        + [gen["joint"]["spring_energy_kn_m"]],
        kkt_energies=[r["physical"], r["artificial"], r["spring"]],
        kkt_membrane_bending_shear=split.tolist(),
        compliance=float(r["load"] @ q),
        midpoint_global_displacement=stations[(len(case["joint"]["line_m"]) - 1) // 2][
            "u_joint_m"
        ],
        nodal_action_reaction_max=float(
            np.abs(
                np.array(
                    [np.r_[s["force_kn"], s["moment_at_joint_kn_m"]] for s in stations]
                )
                .reshape(len(case["joint"]["wings"]), -1, 6)
                .sum(axis=0)
            ).max()
        ),
    )


def series(cli, out):
    out.mkdir(parents=True, exist_ok=False)
    rows = []
    for distorted in [False, True]:
        for n in [2, 4, 8]:
            case = patch(n, distorted=distorted)
            start = time.monotonic()
            data = run(cli, case, out / f"n{n}-distort{int(distorted)}")
            row = compare(case, data)
            row.update(n=n, distorted=distorted, seconds=time.monotonic() - start)
            rows.append(row)
            (out / "series.json").write_text(
                json.dumps(rows, indent=2, allow_nan=False) + "\n", encoding="utf-8"
            )
            print(
                n,
                distorted,
                row["nodal_max_abs"],
                row["compliance"],
                row["seconds"],
                flush=True,
            )
    controls = []
    for law in ["rigid", "hinge", "spring"]:
        for angle in [0.0, 0.1, 90.0]:
            for alpha in [1e-5, 1e-4, 1e-3]:
                case = patch(1, law, alpha, angle)
                data = run(cli, case, out / f"{law}-{angle}-{alpha}")
                row = compare(case, data)
                controls.append(
                    dict(
                        law=law,
                        angle=angle,
                        alpha=alpha,
                        five=five(case),
                        modes=modes(case),
                        six=row,
                    )
                )
    (out / "controls.json").write_text(
        json.dumps(controls, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("cli", type=Path)
    p.add_argument("output", type=Path)
    args = p.parse_args()
    series(args.cli.resolve(), args.output)
