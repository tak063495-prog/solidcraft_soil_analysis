"""Compare tensor, volume and shear MCC work on the same accepted coordinates."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

from cases import SIZES
from c03 import verify_axis
from c03diag_mcc_reference import STAGES
from c03diag_study import read, schedule, write

QUANTITIES = ("tensor", "volume", "shear")


def summarize(rows, force_tolerance):
    times = [(stage, i/32) for stage in STAGES for i in range(1, 33)]
    if [(r["stage"], r["progress"]) for r in rows] != times:
        raise ValueError("incomplete or unordered accepted history")
    if any(not all(math.isfinite(v) for v in (*r["work_kn_m"], r["relative_force_residual"],
               r["max_local_residual_kpa"], r["sum_abs_residual_kn_m"])) or
           len(r["work_kn_m"]) != 3 or r["relative_force_residual"] > force_tolerance for r in rows):
        raise ValueError("invalid work or failed original equilibrium gate")
    work = [math.fsum(r["work_kn_m"][k] for r in rows) for k in range(3)]
    return dict(work_kn_m=work, residual_kn_m=work[0]-(work[1]+work[2]),
                max_increment_residual_kn_m=max(abs(w-(v+s)) for w,v,s in (r["work_kn_m"] for r in rows)),
                max_local_residual_kpa=max(r["max_local_residual_kpa"] for r in rows),
                sum_abs_local_residual_kn_m=math.fsum(r["sum_abs_residual_kn_m"] for r in rows),
                substeps=sum(r["substeps"] for r in rows), plastic_substeps=sum(r["plastic_substeps"] for r in rows),
                max_force_residual=max(r["relative_force_residual"] for r in rows))


def audit_force(out):
    """Apply the survey's published residual criterion to the saved increments."""
    path = out / "increments.csv"
    by_route = {}
    with path.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            key = row["case"] + "-" + row["solver"]
            by_route.setdefault(key, []).append(float(row["relative_force_residual"]))
    if set(by_route) != {c + "-" + s for c in SIZES for s in ("native", "independent")}:
        raise ValueError("incomplete force audit")
    if any(len(values) != 224 or not all(math.isfinite(v) for v in values) for values in by_route.values()):
        raise ValueError("incomplete force history")
    limit = 1e-10
    exceed = {key: sum(v > limit for v in values) for key, values in by_route.items()}
    return dict(
        increments_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        criterion=limit,
        all_within_criterion=not any(exceed.values()),
        exceeded=exceed,
        maxima={key: max(values) for key, values in by_route.items()},
        scopes=[dict(id="force_all_points", status="未達" if any(exceed.values()) else "達成",
                     scope="原因調査の全受理点1e−10を各経路224座標で照合。")],
    )


def compare(root, out):
    out.mkdir()
    histories, totals, pairs, digests = {}, {}, {}, {}
    def fingerprint(p):
        with p.open("rb") as stream:
            digests[p.relative_to(root).as_posix()] = hashlib.file_digest(stream, "sha256").hexdigest()
    for p in (root/"basis/mesh.json", root/"ground.json"):
        fingerprint(p)
    with (out/"increments.csv").open("x", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["case","solver","stage","progress",*(f"{q}_kn_m" for q in QUANTITIES),
                         "residual_kn_m","max_local_residual_kpa","sum_abs_local_residual_kn_m",
                         "substeps","plastic_substeps","relative_force_residual"])
        for case in SIZES:
            folder = root/case
            cfg = read(folder/"stages.json")
            raw = read(folder/"raw/run.json")
            if raw["case"] != cfg or raw["ground"] != read(root/"ground.json"):
                raise ValueError("native input binding differs")
            schedule(raw, 32)
            for p in (folder/"stages.json", folder/"raw/run.json", folder/"raw/accepted.jsonl",
                      folder/"raw/constitutive.jsonl"):
                fingerprint(p)
            for side in ("native","independent"):
                path = folder/(side+"-work.jsonl")
                fingerprint(path)
                rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]
                key = case+"-"+side
                histories[key] = rows
                tolerance = cfg["settings"]["relative_force_tolerance"] if side=="native" else 1e-10
                totals[key] = summarize(rows,tolerance)
                for r in rows:
                    w,v,s = r["work_kn_m"]
                    writer.writerow([case,side,r["stage"],r["progress"],w,v,s,w-(v+s),
                                     r["max_local_residual_kpa"],r["sum_abs_residual_kn_m"],
                                     r["substeps"],r["plastic_substeps"],r["relative_force_residual"]])
    def pair(a,b,axis):
        verify_axis(read(root/a.split("-")[0]/"stages.json"),read(root/b.split("-")[0]/"stages.json"),axis)
        pairs[a+" minus "+b] = dict(axis=axis, signed_total_difference_kn_m=[
            x-y for x,y in zip(totals[a]["work_kn_m"], totals[b]["work_kn_m"], strict=True)],
            maximum_increment_difference_kn_m=[max(abs(x["work_kn_m"][k]-y["work_kn_m"][k])
             for x,y in zip(histories[a],histories[b],strict=True)) for k in range(3)])
    for side in ("native","independent"):
        for case,axis in (("m0c1","material"),("m1c0","contact")):
            pair("m1c1-"+side,case+"-"+side,axis)
    for case in SIZES:
        pair(case+"-native",case+"-independent","solver")
    write(out/"summary.json",dict(components=QUANTITIES,totals=totals,pairs=pairs,inputs_and_traces_sha256=digests))
    write(out/"force-audit.json", audit_force(out))
    # Each ledger row has a distinct scope. No physical acceptance inferred from executing a comparison.
    scopes = [
        ("work_evaluated","達成","本体・独立の6経路、324 GP・224受理座標。全材料内分割の終端状態を保存状態と完全一致照合。"),
        ("decomposition_evaluated","達成","独立に計算した体積・せん断の和との差を、材料内分割・受理増分・全履歴で計測。"),
        ("one_axis_comparisons","達成","同じ時刻の材料2対、接触2対、実装間3対の差を計測。"),
        ("quantitative_work_acceptance","未達","参照なし・未達。仕事量・分解残差・実装差の承認済み定量許容幅なし。"),
        ("general_c03","未達","一般材料・別荷重・開口・全場の刻み収束および元 C-03 全体。")]
    write(out/"validation.json",dict(scopes=[dict(id=k,status=s,scope=t) for k,s,t in scopes],
          measured_complete=True,quantitative_acceptance=False,c03_complete=False,p6_complete=False,
          work_tolerance=None,routes=6,comparisons=7,accepted_times=224,
          summary_sha256=hashlib.sha256((out/"summary.json").read_bytes()).hexdigest(),
          increments_sha256=hashlib.sha256((out/"increments.csv").read_bytes()).hexdigest()))


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root",type=Path)
    parser.add_argument("out",type=Path)
    args=parser.parse_args()
    compare(args.root,args.out)
