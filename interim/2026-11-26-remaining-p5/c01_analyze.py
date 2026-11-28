"""Decompose accepted slip changes on actual parents and score both archived full routes."""
import gzip
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
P4 = HERE.parent / "2026-11-25-remaining-p4/c01"
P3 = HERE.parent / "2026-11-25-remaining-p3/c01/fixed"
SIZES = (256, 512, 1024)
NAMES = ("geostatic", "shear", "open", "recontact", "reverse", "reload", "return")


def projection(a, b):
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    return float(np.vdot(a, b)/denom) if denom else None


def pair(values):
    ab, bc = values[0]-values[1], values[1]-values[2]
    return dict(norm=[np.linalg.norm(ab), np.linalg.norm(bc)], cosine=projection(ab, bc))


def run(out):
    report = dict(levels={}, common_anchors={}, routes={})
    snapshots = {}
    for n in SIZES:
        folder = out / f"n{n}"
        series = {}
        with gzip.open(folder / "trials.jsonl.gz", "rt") as src:
            rows = list(map(json.loads, src))
        # Telescope along the actually retained path. All three corrections compare
        # fixed-parent candidates; their sum never changes subsequent parents.
        cumulative = np.zeros((3, 32, 3))
        start_slip = np.array(rows[0]["parent"]["slip"])
        summaries, selected, closures = {}, {}, {}
        for idx, row in enumerate(rows):
            field = lambda key: np.array(row[key]["slip"])
            pre = field("before_c") - field("before_f")
            post = field("after_c") - field("after_f")
            assert abs(post).max() == row["post_max_slip"]
            end = row["end"] * 256
            region = "before_closure" if end <= 235 else "closure" if end <= 239 else "after_closure"
            a = summaries.setdefault(region, dict(trials=0, accepted=0, selected=[0,0,0],
                slip_over=0, dissipation_over=0, friction_event=0, changed_status_trials=0,
                before_max_m=0., after_max_m=0., accepted_before_max_m=0., accepted_after_max_m=0.,
                fixed_only_slip_over=0, fixed_only_slip_without_other_veto=0, local_direction_reversals=0))
            # Rescore only these candidates. The scaled trajectory has different parents.
            fixed_only = row["slip_budget"] < row["post_max_slip"] <= 1.e-4/n**2
            a["fixed_only_slip_over"] += fixed_only
            a["fixed_only_slip_without_other_veto"] += (fixed_only and row["event"] is None
                and row["post_max_dissipation"] <= row["dissipation_budget"])
            a["local_direction_reversals"] += int(np.vdot(pre, post) < 0)
            a["trials"] += 1
            a["accepted"] += row["accepted"]
            a["selected"] = [x+int(y is not None) for x,y in zip(a["selected"],row["selected_from"])]
            a["slip_over"] += row["post_max_slip"] > row["slip_budget"]
            a["dissipation_over"] += row["post_max_dissipation"] > row["dissipation_budget"]
            a["friction_event"] += row["event"] is not None
            switch = any(row[k]["status"] != row["parent"]["status"] for k in ("after_m","after_f","after_c"))
            a["changed_status_trials"] += switch
            for tag, diff in (("before",pre),("after",post)):
                val = abs(diff).max()
                a[tag+"_max_m"] = max(a[tag+"_max_m"],val)
                if row["accepted"]:
                    key = "accepted_"+tag+"_max_m"
                    a[key] = max(a[key],val)
            examples = {"max_pre":abs(pre).max(), "max_post":abs(post).max(),
                "max_correction":abs(post-pre).max()}
            if switch:
                examples["max_switch_post"] = abs(post).max()
            if row["accepted"]:
                examples["accepted_max_post"] = abs(post).max()
                examples["accepted_max_fine_correction"] = abs(field("after_f")-field("before_f")).max()
            for label, measure in examples.items():
                if label not in selected or measure > selected[label]["ranking_value"]:
                    gp, comp = np.unravel_index(abs(post).argmax(),post.shape)
                    selected[label] = dict(ranking_value=measure, trial_index=idx,
                        start=row["start"],end=row["end"],accepted=row["accepted"],
                        selected_from=row["selected_from"],gp=int(gp),component=int(comp),
                        old_mid_fine_status=[row[k]["status"][gp] for k in ("parent","after_m","after_f")],
                        before_error_m=pre[gp,comp],after_error_m=post[gp,comp],
                        before_max_m=abs(pre).max(),after_max_m=abs(post).max())
            if row["accepted"]:
                for gp, (left, middle, right) in enumerate(zip(row["parent"]["status"],
                        row["after_m"]["status"], row["after_f"]["status"])):
                    if left != middle or middle != right:
                        closures.setdefault(str(gp), []).append(dict(start=row["start"],
                            end=row["end"],old=left,mid=middle,fine=right))
                cumulative[0] += field("before_f") - field("parent")
                cumulative[1] += field("f_from_corrected_mid") - field("before_f")
                cumulative[2] += field("after_f") - field("f_from_corrected_mid")
                if end.is_integer():
                    np.testing.assert_allclose(start_slip+cumulative.sum(axis=0),field("after_f"),
                                               rtol=0,atol=5e-18)
                    series[int(end)] = (field("after_f"), cumulative.copy())
        snapshots[n] = series
        report["levels"][n] = dict(regions=summaries,examples=selected,contact_changes_by_gp=closures,
            slip_budget_m=rows[0]["slip_budget"],dissipation_budget=rows[0]["dissipation_budget"])
    for i in range(129,257):
        slips, corrections = zip(*(snapshots[n][i] for n in SIZES))
        cumulative = np.array(corrections)
        actual = np.array(slips)
        mid = cumulative[:,1]
        fine = cumulative[:,2]
        raw = actual - mid - fine
        report["common_anchors"][i] = dict(actual=pair(actual), subtract_local_corrections=pair(raw),
            accumulated_midpoint_effect=pair(mid), accumulated_fine_effect=pair(fine),
            correction_norm_by_level=[np.linalg.norm(x) for x in mid+fine],
            signed_correction_projection_fraction=[
                float(np.vdot(a-b, c-d) / np.vdot(a-b, a-b)) if np.any(a!=b) else None
                for a,b,c,d in zip(actual[:-1],actual[1:],(mid+fine)[:-1],(mid+fine)[1:])])
    for name,path in (("original",P4/"original-fixed-slip.npz"),("alternate",P3/"slip-common.npz")):
        with np.load(path) as archive:
            vals = np.array([archive[f"n{n}"] for n in SIZES])
        report["routes"][name] = {stage:pair(vals[:,i*256:(i+1)*256]) for i,stage in enumerate(NAMES)}
        report["routes"][name]["all"] = pair(vals)
        report["routes"][name]["endpoints"] = {stage:pair(vals[:,(i+1)*256-1]) for i,stage in enumerate(NAMES)}
    comparison = json.loads((P4/"summary.json").read_text())
    report["archived_budget_controls"] = {}
    for route in ("orig_fixed","orig_scaled","alt_fixed","alt_scaled"):
        data = comparison["routes"][route]
        report["archived_budget_controls"][route] = dict(overall=data["overall"], by_stage={
            key:{k:v[k] for k in ("cosine","coarse_norm","fine_norm","attempts","rebalances")}
            for key,v in data["by_stage"].items()})
    (out/"analysis.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")


if __name__ == "__main__":
    run(Path(sys.argv[1]))
