"""Resolve the selected projected enclosure into additive leaf terms (metres)."""
import gzip
import csv
import json
from pathlib import Path
import sys

import numpy as np


def score(folder, sign):
    prefix = f"{sign}-833"
    record = json.loads((folder/(sign+".json")).read_text())
    step = record["steps"][-1]
    assert step["step"] == 833 and not step["included"]
    inv = step["inverse_test"]
    def get(name):
        with np.load(folder/(prefix+"-"+name+".npz")) as z:
            return dict(z)
    leaves, pair, raw = get("leaves"), get("pair"), get("pair-supports")
    center = get("direct")["center_after_inverse_m"]
    norm = np.linalg.norm
    base = norm(center+pair["mid_action_m"])
    c = pair["component_radius_m"]
    assert norm(c) < pair["remainder_m"]
    e = raw["component_norm_share_m"]
    np.testing.assert_allclose(raw["summed_component_m"],c,rtol=1e-12,atol=1e-25)
    np.testing.assert_allclose(e.sum(),norm(c),rtol=1e-12)
    # Rationalize sqrt(base^2 + x) - base. This gives an exact additive
    # accounting for this evaluated bound, including its signed cross term.
    rad = norm(c)
    along = raw["directional_m"]
    x = rad**2+2*(norm(center)*along.sum()+norm(pair["mid_action_m"])*rad)
    hist = (e*(rad+2*norm(pair["mid_action_m"]))+2*norm(center)*along)/(np.sqrt(base**2+x)+base)
    hist *= 1+1e-9
    selected = leaves["selected_indices"]
    nom = leaves["tile_remainder_m"].sum(1)
    nom[selected] = leaves["child_tile_remainder_m"].sum(1)
    nom *= (1+1e-14)*(1+1e-9)
    key = lambda arr: [tuple(map(int,k)) for k in arr]
    by_nom = dict(zip(key(leaves["block_leaf"]),nom,strict=True))
    by_pair = {k:i for i,k in enumerate(key(raw["block_leaf"]))}
    refined = set(key(leaves["block_leaf"][selected]))
    rows = []
    for b,k in sorted(set(by_nom)|set(by_pair)):
        idx = by_pair.get((b,k))
        def support(a):
            return float(a[idx]*1e9) if idx is not None else 0.
        rows.append(dict(block=b,leaf=k,refined=(b,k) in refined,nominal_nm=1e9*by_nom.get((b,k),0.),
            norm_support_nm=support(raw["scalar_m"]),direction_support_nm=support(along),
            component_norm_share_nm=support(e),projected_history_share_nm=support(hist)))
    fixed = (base*(1+1e-9)+inv["force_guard_m"])*1e9+10*inv["contraction_upper"]
    def stats(rows):
        return dict(count=len(rows),**{k:sum(r[k] for r in rows) for k in
                    ("nominal_nm","norm_support_nm","direction_support_nm","projected_history_share_nm")})
    all_terms = stats(rows)
    image = inv["ball_image_upper_m"]*1e9
    np.testing.assert_allclose(fixed+all_terms["nominal_nm"]+all_terms["projected_history_share_nm"],image,
                               rtol=3e-13)
    assert inv["ball_image_upper_m"]==inv["projected_image_m"]
    for row in rows:
        row["total_nm"] = row["nominal_nm"]+row["projected_history_share_nm"]
    ranked = sorted(rows,key=lambda r:-r["total_nm"])
    cumsum = np.cumsum([r["total_nm"] for r in ranked])
    exceed = int(np.searchsorted(cumsum,10-fixed,side="right"))+1
    groups = {f"block{b}_{'refined' if ref else 'unrefined'}":
              stats([r for r in rows if r["block"]==b and r["refined"]==ref])
              for b in (0,1) for ref in (False,True)}
    with gzip.open(folder/(prefix+"-by-leaf.csv.gz"),"wt") as f:
        writer = csv.DictWriter(f,fieldnames=list(ranked[0]),lineterminator="\n")
        writer.writeheader()
        writer.writerows(ranked)
    return dict(image_nm=image,excess_nm=image-10,fixed_nm=fixed,all=all_terms,groups=groups,
                top_100=stats(ranked[:100]),ranked_prefix_exceeding_10_nm=exceed,
                image_before_last_rank_nm=fixed+cumsum[exceed-2],image_at_last_rank_nm=fixed+cumsum[exceed-1],
                largest_leaf=ranked[0],consecutive=record["consecutive"],
                full_history_relative_upper=record["full_history_relative_upper"])


if __name__ == "__main__":
    root = Path(sys.argv[1])
    results = {f"{label}-{s}":score(root/label,s) for label,s in
               (("baseline","minus"),("baseline","plus"),("tight","minus"))}
    def history(label):
        with gzip.open(root/label/"minus-833-by-leaf.csv.gz","rt") as f:
            return {(int(r["block"]),int(r["leaf"])):r for r in csv.DictReader(f)}
    before, after = history("baseline"),history("tight")
    with gzip.open(root/"radius-effect-by-leaf.csv.gz","wt") as f:
        writer = csv.writer(f,lineterminator="\n")
        cols = ("nominal_nm","norm_support_nm","direction_support_nm","projected_history_share_nm")
        writer.writerow(["block","leaf",*["decrease_"+k for k in cols]])
        for key in sorted(before.keys()|after.keys()):
            writer.writerow([*key,*[float(before.get(key,{}).get(k,0))-float(after.get(key,{}).get(k,0)) for k in cols]])
    (root/"leaf-groups.json").write_text(json.dumps(results,indent=2,allow_nan=False)+"\n")
