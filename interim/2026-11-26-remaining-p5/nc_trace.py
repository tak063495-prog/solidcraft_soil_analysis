"""Replay two archived branches up to their material diagnostics, without extending either."""
import gzip
import hashlib
import inspect
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
P4 = HERE.parent / "2026-11-25-remaining-p4"
P3 = HERE.parent / "2026-11-25-remaining-p3"
sys.path.insert(0, str(P4))
import nc_extended
import nc_branch
from nc_restart import read


def save(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def point_trace(binary, row, mode, out):
    data = dict(row, mode=mode)
    proc = subprocess.run([binary], input=json.dumps(data), text=True, capture_output=True,
                          env=dict(os.environ, P5_POINT_TRACE="1"), check=True)
    steps = list(map(json.loads, proc.stderr.splitlines()))
    save(out, dict(result=json.loads(proc.stdout), trace=steps))
    return json.loads(proc.stdout)


def retry_series(work, out):
    """Check the two next rejections and the first successful cutback at this parent."""
    binary = work / "probe/target/release/examples"
    capture = out / "retry-inputs.jsonl"
    os.environ.update(NC_WORKER=str(binary / "p5_homogeneous"), NC_FAILURES=str(capture))
    m = nc_extended.Extended(work / "inputs", P4 / "nc_capture.py")
    parent = read(out / "parent.json.gz")
    state, lam, t = parent["active_material_states"], parent["load"], np.array(parent["secant"])
    rows = []
    for row in map(json.loads, (P4 / "nc/extended-004/trials.jsonl").open()):
        if row["attempt"] not in (107, 108, 109):
            continue
        assert row["parent_lambda"] == lam
        result, it, reason = nc_branch.correct(m, state, lam, t, row["ds_m"], 1e-12)
        assert it == row["iterations"] and reason == row.get("reason")
        assert row["accepted"] == (result is not None)
        if result is None:
            failed = json.loads(capture.read_text().splitlines()[-1])
            assert failed["state"] == state[351] and failed["point"] == 351
            check = point_trace(binary/"p5_point", failed, "homogeneous",
                                out/f"line-search-{row['attempt']}.json")
            assert "midpoint local line search failed" in check["reason"]
        else:
            dz, _, diag = result
            assert lam+dz[-1] == row["load"] and all(diag[k] == row[k] for k in diag)
        rows.append(row)
        if result is not None:
            break
    assert [r["accepted"] for r in rows] == [False,False,True]
    save(out/"matched-retries.json", rows)


def run(work, out):
    out.mkdir(parents=True)
    binaries = work / "probe/target/release/examples"
    capture = out / "failed-inputs.jsonl"
    os.environ.update(NC_WORKER=str(binaries / "p5_homogeneous"), NC_FAILURES=str(capture))
    m = nc_extended.Extended(work / "inputs", P4 / "nc_capture.py")
    parent_file = P3 / "nc/0.02/last.json.gz"
    parent = read(parent_file)
    state, lam = parent["active_material_states"], parent["load"]
    u = np.ravel(parent["displacement_m"])
    t = np.load(P4 / "nc/fine-secant.npy")
    matched = []
    history = P4 / "nc/extended-004/trials.jsonl"
    for line in history.open():
        row = json.loads(line)
        assert row["attempt"] < 1600 and lam == row["parent_lambda"]
        ds = row["ds_m"]
        result, it, reason = nc_branch.correct(m, state, lam, t, ds, 1e-12)
        assert row["accepted"] == (result is not None) and it == row["iterations"]
        assert reason == row.get("reason")
        matched.append(dict(attempt=row["attempt"], parent_lambda=lam, ds_m=ds,
                            accepted=row["accepted"], reason=reason))
        if result is None:
            if "active point 351" in reason and "midpoint local line search failed" in reason:
                with gzip.open(out / "parent.json.gz", "wt") as stream:
                    json.dump(dict(displacement_m=u.tolist(), active_material_states=state,
                                   load=lam, secant=t.tolist(), ds=ds), stream, allow_nan=False)
                break
            continue
        dz, state, diag = result
        u[m.free] += dz[:-1]
        lam += dz[-1]
        assert lam == row["load"] and all(diag[k] == row[k] for k in diag)
        t = dz / np.linalg.norm(dz)
    else:
        raise ValueError("missing point 351 failure")
    failed = json.loads(capture.read_text().splitlines()[-1])
    assert failed["point"] == 351 and failed["state"] == state[351]
    result = point_trace(binaries / "p5_point", failed, "homogeneous", out / "line-search.json")
    assert "midpoint local line search failed" in result["reason"]
    save(out / "matched-trials.json", dict(source=str(history.relative_to(HERE.parents[1])),
        start_sha256=hashlib.sha256(parent_file.read_bytes()).hexdigest(), rows=matched))

    retry_series(work, out)

    # Same round-trip computation as the archived auxiliary branch, with a pressure observer.
    pressure = []
    def observe(load, old):
        pc = old[351]["preconsolidation_pressure_kpa"]
        scaled = math.ldexp(pc, 448)
        pressure.append(dict(load=load, pc_kpa=pc, scaled_square=scaled**2,
                             scaled_square_subnormal=0 < scaled**2 < sys.float_info.min))
    source = inspect.getsource(nc_branch.run)
    for line in ('    checkpoint("fold")', '                accepted += 1'):
        assert source.count(line) == 1
        source = source.replace(line, line+"\n"+line[:len(line)-len(line.lstrip())]+"observe(lam, old)")
    namespace = dict(vars(nc_branch), observe=observe)
    exec(compile(source, nc_branch.__file__, "exec"), namespace)
    replay = work / "aux-replay"
    namespace["run"](work / "inputs", binaries / "p5_units", replay, 0.04)
    end = read(replay / "last.json.gz")
    assert end == read(P3 / "nc/roundtrip-0.04/last.json.gz")
    save(out / "pc-history.json", pressure)
    frozen = list(map(json.loads,(P4 / "nc/aux448-incoming-failures.jsonl").read_text().splitlines()))
    results = []
    for i, record in enumerate(frozen):
        assert record["point"] == 351 and record["state"] == end["active_material_states"][351]
        for mode in ("units448", "comparison480"):
            result = point_trace(binaries / "p5_point",record,mode,out / f"yield-{i+1}-{mode}.json")
            assert (result["reason"] is None) == (mode == "comparison480")
            results.append(result)
    save(out / "unit-comparison.json", results)
    save(out / "attempts.json",dict(extended_prefix=len(matched), homogeneous_point=1,
        auxiliary_including_fold=json.loads((replay/"summary.json").read_text())["attempts"]+1,
        auxiliary_frozen_local_calls=len(results), diagnostic_480_global_calls=0,
        budget=16520, default=2000))


if __name__ == "__main__":
    run(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve())
