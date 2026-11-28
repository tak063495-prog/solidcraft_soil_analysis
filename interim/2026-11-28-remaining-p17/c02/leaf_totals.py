"""Collect the signed-path stopping metrics from the raw Bernstein leaf tables."""
import json
from pathlib import Path
import sys

import numpy as np


def run(out):
    rows = {}
    for sign in ("plus","minus"):
        history = json.loads((out/f"{sign}.json").read_text())
        for row in history["steps"]:
            step = row["step"]
            path = out/f"{sign}-{step}-leaves.npz"
            if not path.is_file():
                continue
            with np.load(path) as data:
                sums = data["tile_remainder_m"].sum(1)*1e9
                win = int(data["selected_index"])
                after = data["child_tile_remainder_m"].sum()*1e9
                before = sums[win]
                nom = row["center_only_inverse_test"]
                inv = row["inverse_test"]
                updated = sums.copy()
                updated[win] = after
                excess = inv["ball_image_upper_m"]*1e9-10
                # Diagnostic only. Hold the image direction, history term and qR fixed.
                ranked = np.cumsum(np.sort(updated)[::-1])
                rank = np.searchsorted(ranked,excess)
                rows[f"{sign}-{step}"] = dict(active=len(sums),total_leaves=int(data["roots_per_block"].sum()),
                    block_leaf=data["block_leaf"][win].tolist(),prior_remainder_nm=sums.sum(),
                    largest_leaf_nm=before,largest_share_percent=100*before/sums.sum(),split_remainder_nm=after,
                    new_remainder_nm=nom["direct_center_remainder_m"]*1e9,center_image_nm=nom["ball_image_upper_m"]*1e9,
                    history_image_nm=inv["ball_image_upper_m"]*1e9,contraction=inv["contraction_upper"],
                    consecutive=history["consecutive"],excess_nm=excess,
                    minimum_rank_count_if_only_remainder_removed_at_fixed_center_and_history=(
                        0 if excess<=0 else int(rank+1) if rank<len(ranked) else None))
    (out/"leaf-summary.json").write_text(json.dumps(rows,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    run(Path(sys.argv[1]))
