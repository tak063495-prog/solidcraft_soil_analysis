"""Locate the first negative signed-slip refinement product on the common load grid."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "tools"), str(ROOT / "interim/2026-11-24-c01-event4-retry")]
from c01_retry_fields import CONTACT, STAGES
from c01_retry_metrics import History, align
from joint_reference import System
from verify_joint_coupling import contact


def diagnose(folder, side, mode, anchor_index=None, point=None, save_differences=None):
    histories = []
    for n in (256, 512, 1024):
        key = f"{side}-{n}-{mode}"
        stages = json.loads((folder/key/"result.json").read_bytes())["stages"]
        histories.append(History(folder/(key+"-fields"), side, n, stages))
    a, b, c = align(histories, "plastic_slip", 256)
    ab, bc = a-b, b-c
    if save_differences is not None:
        with save_differences.open("xb") as sink:
            np.savez_compressed(sink, coarse_medium_m=ab, medium_fine_m=bc)
    products = np.einsum("tqc,tqc->tq", ab, bc)
    times = products.sum(axis=1)
    def coordinate(index):
        stage, step = divmod(int(index), 256)
        return [STAGES[stage], step+1, 256]
    negative = np.flatnonzero(times < 0)
    accumulated = np.flatnonzero(times.cumsum() < 0)
    result = dict(common_anchors=len(a), first_negative=None,
                  first_negative_cumulative=coordinate(accumulated[0]) if len(accumulated) else None)
    if len(negative):
        result["first_negative"] = coordinate(negative[0])
    if anchor_index is None and not len(negative):
        return result
    idx = int(negative[0] if anchor_index is None else anchor_index)
    if not 0 < idx < len(a):
        raise ValueError("diagnostic anchor must have a preceding common anchor")
    gp = int(products[idx].argmin() if point is None else point)
    result.update(measured_anchor=coordinate(idx), point=gp, anchor_product_m2=float(times[idx]),
                  point_product_m2=float(products[idx,gp]), coarse_medium_m=ab[idx,gp].tolist(),
                  medium_fine_m=bc[idx,gp].tolist(), grids={})
    mesh = json.loads((ROOT/"examples/joint-water-column-tet10-mesh.json").read_bytes())
    case = json.loads((folder/f"{side}-256-{mode}"/"stages.json").read_bytes())
    system = System(mesh, case["joint_system"])
    for history in histories:
        n = history.n
        stride = n//256
        # Include the preceding common anchor, then each intermediate native grid anchor.
        interval = slice(idx*stride-1, (idx+1)*stride)
        paths = {f: history.grid(f,n)[interval] for f in
                 ("jump","gap","traction","plastic_slip","contact_status")}
        result["grids"][str(n)] = dict(status=[CONTACT[i] for i in paths["contact_status"][:,gp]],
                gap_m=paths["gap"][:,gp].tolist(), traction_kpa=paths["traction"][:,gp].tolist(),
                jump_m=paths["jump"][:,gp].tolist(), slip_m=paths["plastic_slip"][:,gp].tolist(),
                increment_m=(paths["plastic_slip"][-1,gp]-paths["plastic_slip"][0,gp]).tolist())
        if n == 256:
            parent_jump = paths["jump"][0].copy()
            parent_tau = paths["traction"][0].copy()
        # Hold the initial state AND endpoint fixed; insert measured intermediate equilibria only.
        if n == 256:
            fixed_end = paths["jump"][-1].copy()
        targets = list(paths["jump"][1:-1]) + [fixed_end]
        jump, tau = parent_jump, parent_tau
        flow = np.zeros_like(jump)
        for target in targets:
            q, ds, _ = contact(system, jump, tau, target, 1)
            tau, jump = q[:,0], target
            flow += ds[:,0]
        result["grids"][str(n)]["fixed_parent_chord_replay_increment_m"] = flow[gp].tolist()
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    parser.add_argument("side", choices=("native","independent"))
    parser.add_argument("mode", choices=("friction-events","resolved-slip"))
    parser.add_argument("--anchor-index",type=int,help="zero-based anchor for a paired event replay")
    parser.add_argument("--point",type=int,help="zero-based interface quadrature point")
    parser.add_argument("--save-differences",type=Path,help="save both signed difference fields losslessly")
    args = parser.parse_args()
    print(json.dumps(diagnose(args.folder,args.side,args.mode,args.anchor_index,args.point,args.save_differences), indent=2, allow_nan=False))
