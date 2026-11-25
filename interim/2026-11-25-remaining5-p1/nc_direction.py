"""Assemble actual forward-increment MCC tangents on the replayed NC history."""

import gzip
import hashlib
import json
from pathlib import Path
import sys
import zipfile
import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from verify_high_order import prepare


def acoustic_matrix(c, n):
    x, y, z = n
    # Engineering-strain rank-one map. Q = B(n)^T C B(n).
    b = np.array([[x, 0, 0], [0, y, 0], [0, 0, z], [y, x, 0], [0, z, y], [z, 0, x]])
    return b.T @ c @ b


def acoustic(c, n):
    q = acoustic_matrix(c, n)
    return np.linalg.eigvalsh((q + q.swapaxes(-1, -2)) / 2)[:, 0]


def main(work, out):
    out.mkdir()
    with zipfile.ZipFile(
        HERE.parent / "2026-11-23-nc01-coordinate-polish/evidence.zip"
    ) as z:
        archived = json.loads(z.read("run.json"))
    raw_run = (work / "nc-run/run.json").read_text()
    raw_state = raw_run.split('"last_committed":', 1)[1]
    _, stop = json.JSONDecoder().raw_decode(raw_state)
    expected_parent = hashlib.sha256(raw_state[:stop].encode()).hexdigest()
    run = json.loads(raw_run)
    assert run["last_committed"] == archived["last_committed"]

    def strip_time(a):
        if isinstance(a, dict):
            return {k: strip_time(v) for k, v in a.items() if k != "elapsed_seconds"}
        if isinstance(a, list):
            return list(map(strip_time, a))
        return a

    assert [strip_time(a) for a in run["attempts"]] == [
        strip_time(a) for a in archived["attempts"]
    ]
    assert run["failure_reason"] == archived["failure_reason"]
    lam = run["last_committed"]["stage_progress"]
    assert lam == 0.5647809505462646
    mesh = json.loads((work / "inputs/mesh.json").read_bytes())
    quadrature, offsets = prepare(mesh, run["integration_points"])
    live = np.array(run["last_committed"]["active_elements"])
    hold = np.array(run["excavation_audits"][0]["initial_holding_load_kn"])
    h = run["case"]["settings"]["tangent_perturbation"]
    pairs = [
        a for a in run["attempts"] if a["stage"] == "excavate_upper" and a["accepted"]
    ][-5:] + [run["attempts"][-1]]
    prior = {}
    table = []
    terminal = json.loads(
        gzip.decompress(
            (
                HERE.parent
                / "2026-11-25-remaining5-diagnosis/nc-terminal-system.json.gz"
            ).read_bytes()
        )
    )
    for idx, attempt in enumerate(pairs):
        a, b = attempt["progress_from"], attempt["progress_to"]
        key = lambda x: np.array(x, dtype=np.float64).view(np.uint64).item()
        files = list(work.glob(f"nc-{key(a):016x}-{key(b):016x}-*.json"))
        path = max(files, key=lambda f: int(f.stem.split("-")[-1]))
        data = json.loads(path.read_bytes())
        assert (data["from"], data["to"]) == (a, b)
        order = sorted((r[0][0], g) for g, r in enumerate(data["rows"]) if r)
        assert all(len(r) == 1 and r[0][1] == 1 for r in data["rows"] if r)
        free = np.array([g for _, g in order])
        assert len(free) == 348
        du = (np.array(data["u"]) - data["old_u"]).ravel()[free]
        delta_load = (np.array(data["target_load"]) - data["old_load"])[free]
        hf = hold[free]
        np.testing.assert_allclose(delta_load, -(b - a) * hf, rtol=1e-8, atol=3e-12)
        points = data["points"]
        positive, negative, steps, ids = (
            np.array([r[k] for r in points])
            for k in ("positive", "negative", "steps", "id")
        )
        row = dict(
            parent=a,
            target=b,
            accepted=attempt["accepted"],
            iteration=data["iteration"],
            original_force=data["relative"],
            parent_sha256=data["parent_hash"],
            du_m=np.linalg.norm(du),
            load_increment_work_kn_m=du @ delta_load,
            material_points=len(ids),
            positive_plastic=int(sum(positive > 0)),
            negative_plastic=int(sum(negative > 0)),
            forward_only=int(sum((positive > 0) & (negative == 0))),
            backward_only=int(sum((negative > 0) & (positive == 0))),
            tangents={},
        )
        row["adaptive_width_range"] = [
            min(r["width"] for r in points),
            max(r["width"] for r in points),
        ]
        row["one_sided_columns"] = sum(r["single"] for r in points)
        row["mixed_columns"] = sum(r["unresolved"] for r in points)
        archive = dict(
            free=free,
            du=du,
            holding=hf,
            point_id=ids,
            positive=positive,
            negative=negative,
            partitions=steps,
            strain=np.array([r["strain"] for r in points]),
        )
        for name, width in (("c", h), ("narrow", "adaptive_branch_stencil")):
            tensors = (
                np.array([p[name] for p in points]).reshape(-1, 6, 6).transpose(0, 2, 1)
            )
            k = np.zeros((hold.size, hold.size))
            for ei in np.flatnonzero(live):
                q = quadrature[ei]
                source = np.arange(offsets[ei], offsets[ei + 1])
                take = np.searchsorted(ids, source)
                np.testing.assert_array_equal(ids[take], source)
                dofs = (3 * q["nodes"][:, None] + np.arange(3)).ravel()
                k[np.ix_(dofs, dofs)] += np.einsum(
                    "qai,qab,qbj,q->ij",
                    q["b"],
                    tensors[take],
                    q["b"],
                    q["weight"],
                    optimize=True,
                )
            k = k[np.ix_(free, free)]
            left, s, vh = np.linalg.svd(k)
            w = left[:, -1]
            ori = w @ prior[name] if name in prior else w[np.argmax(abs(w))]
            if ori < 0:
                w = -w
            prior[name] = w
            scan = []
            active = tensors[positive > 0]
            for count in (64, 256):
                j = np.arange(count)
                z = (j + 0.5) / count
                angle = j * np.pi * (3 - np.sqrt(5.0))
                normals = np.column_stack(
                    (
                        np.sqrt(1 - z * z) * np.cos(angle),
                        np.sqrt(1 - z * z) * np.sin(angle),
                        z,
                    )
                )
                best = (np.inf, None, None)
                for direction in np.vstack((np.eye(3), normals)):
                    values = acoustic(active, direction)
                    p = int(np.argmin(values))
                    if values[p] < best[0]:
                        best = (
                            values[p],
                            int(ids[positive > 0][p]),
                            direction.tolist(),
                        )
                scan.append(
                    dict(
                        directions=count + 3,
                        min_symmetric_kpa=best[0],
                        point=best[1],
                        normal=best[2],
                    )
                )
            point = scan[-1]["point"]
            local = tensors[np.searchsorted(ids, point)]
            directions = np.vstack((np.eye(3), normals))
            determinants = np.array(
                [np.linalg.det(acoustic_matrix(local, n)) for n in directions]
            )
            hi, lo = np.argmax(determinants), np.argmin(determinants)
            root = dict(
                point=point,
                bracket_determinants=[determinants[lo], determinants[hi]],
                bracket_normals=directions[[lo, hi]].tolist(),
            )
            if determinants[lo] < 0 < determinants[hi]:
                n0, n1 = directions[[lo, hi]]
                if n0 @ n1 < 0:
                    n1 = -n1

                def normal(t):
                    n = (1 - t) * n0 + t * n1
                    return n / np.linalg.norm(n)

                zero = brentq(
                    lambda t: np.linalg.det(acoustic_matrix(local, normal(t))),
                    0,
                    1,
                    xtol=1e-14,
                )
                q0 = acoustic_matrix(local, normal(zero))
                _, sq, vq = np.linalg.svd(q0)
                root.update(
                    normal=normal(zero).tolist(),
                    singular_values=sq.tolist(),
                    polarization=vq[-1].tolist(),
                )
            info = dict(
                acoustic_zero=root,
                stencil=width,
                smin=s[-1],
                second=s[-2],
                smax=s[0],
                w_dot_holding_kn=w @ hf,
                normalized_projection=w @ hf / np.linalg.norm(hf),
                min_symmetric=np.linalg.eigvalsh((k + k.T) / 2)[0],
                incremental_quadratic_work_kn_m=du @ k @ du,
                acoustic_scan=scan,
            )
            if not attempt["accepted"] and name == "c":
                original = np.array(terminal["a_column_major"]).reshape(
                    348, 348, order="F"
                )
                info["terminal_matrix_relative_difference"] = np.linalg.norm(
                    k - original
                ) / np.linalg.norm(original)
                assert info["terminal_matrix_relative_difference"] < 1e-12
            row["tangents"][name] = info
            archive.update(
                {name + "_matrix": k, name + "_left": w, name + "_local": tensors}
            )
        if a == lam:
            assert data["parent_hash"] == expected_parent
            assert data["old_u"] == archived["last_committed"]["displacement_m"]
        np.savez_compressed(out / f"increment-{idx}.npz", **archive)
        table.append(row)
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                original_failure=run["failure_reason"],
                identical_attempts=len(run["attempts"]),
                final_state_equal=True,
                rows=table,
            ),
            indent=2,
            allow_nan=False,
        )
        + "\n"
    )


if __name__ == "__main__":
    main(*map(Path, sys.argv[1:]))
