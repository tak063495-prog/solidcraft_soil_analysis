"""Bound the history difference at identical material positions through A D^T."""
import inspect
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
P15 = HERE.with_name("2026-11-28-remaining-p15")
sys.path.insert(0, str(P15))
import face_direct as p15_face

chain = p15_face.chain


def paired(red, refined, displacements, radii, leaf_box, kt, km, A, ix, direction, dest):
    mats = chain.loc.subdivision()
    restrict = (mats[:, None] @ mats[None]).reshape(16, 15, 15)
    norm_a = np.linalg.norm(A) * (1 + 1e-12)
    midpoint_force = np.zeros_like(displacements[-1]).reshape(-1, 3)
    nodal = np.zeros_like(midpoint_force)
    remainder = 0.
    component = np.zeros(len(A))
    along = 0.
    tiles = 0
    for block, mesh in zip(red.blocks, refined, strict=True):
        sl, _, ids, *_ = block
        _, weights, _, sign, areas, gain2 = mesh
        active = np.flatnonzero(np.any((leaf_box[0][sl] != 0) | (leaf_box[1][sl] != 0), axis=1))
        weights = weights.reshape(-1, 15, 15)
        jumps = [(u.reshape(-1, 3)[ids[1]] - u.reshape(-1, 3)[ids[0]]) * [1, 1, sign]
                 for u in displacements]
        columns = [(A[:, ix[3*ids[1]+d]] - A[:, ix[3*ids[0]+d]]).T for d in (0, 1)]
        for k in active:
            cells = restrict @ weights[k]
            size = areas[k] / len(cells)
            gain = np.sqrt(gain2[k]) * (1 + 1e-12)
            z = np.zeros((len(cells), 2))
            nominal, actual = (z, z), (z, z)
            distance = np.zeros(len(cells))
            for jump, u, radius in zip(jumps, displacements, [*radii, 0.], strict=True):
                coeff = cells @ jump
                pad = 1e-11 * abs(u).max() + 1e-90
                box = (coeff.min(1) - pad, coeff.max(1) + pad)
                uncertainty = gain*radius
                enlarged = box[0]-uncertainty, box[1]+uncertainty
                assert enlarged[1][:, 2].max() < 0
                nominal, stick_n = chain.project(box, nominal, kt, km)
                actual, stick_a = chain.project(enlarged, actual, kt, km)
                # At a common x the disk projection is nonexpansive in history.
                # If both stick, the old difference is copied exactly.
                distance += np.where(stick_n & stick_a, 0., (kt+km)*uncertainty)
                distance *= 1 + 1e-12
            pad = 1e-12 * (np.maximum(abs(actual[0]), abs(actual[1])) +
                           np.maximum(abs(nominal[0]), abs(nominal[1]))) + 1e-80
            lo = np.maximum.reduce([np.broadcast_to(leaf_box[0][sl][k], z.shape),
                                    actual[0]-nominal[1]-pad, np.broadcast_to(-distance[:, None], z.shape)])
            hi = np.minimum.reduce([np.broadcast_to(leaf_box[1][sl][k], z.shape),
                                    actual[1]-nominal[0]+pad, np.broadcast_to(distance[:, None], z.shape)])
            assert (lo <= hi).all()
            mid = (hi+lo)/2
            width = np.maximum(hi-mid, mid-lo)*(1+1e-12) + 1e-80
            use_ball = distance < np.linalg.norm(width, axis=1)
            mid[use_ball] = 0.
            rad = np.minimum(distance, np.linalg.norm(width, axis=1))*(1+1e-12)
            width[use_ball] = np.minimum(distance[use_ball, None], np.maximum(abs(lo),abs(hi))[use_ball])
            a, b = (cells @ v for v in columns)
            norms = np.hypot(a,b)
            amp_guard = 1e-10*norm_a*gain
            amps = np.sqrt((norms**2).sum(2).max(1)) + amp_guard
            remainder += size * (amps @ rad)
            za, zb = (cells @ (v @ direction) for v in columns)
            directional = np.minimum((np.hypot(za,zb)*rad[:,None]).max(1),
                (abs(za)*width[:,0,None]+abs(zb)*width[:,1,None]).max(1))
            along += size * (directional + amp_guard*rad).sum()
            # Bound each component of the image on the very same subtriangle.
            # Convexity bounds each directional support by its largest coefficient.
            sphere = norms * rad[:, None, None]
            rect = abs(a)*width[:,0,None,None] + abs(b)*width[:,1,None,None]
            component += size*(np.minimum(sphere.max(1),rect.max(1)).sum(0) +
                               amp_guard*rad.sum())
            m = size * cells.mean(1).T @ mid
            e = size * np.max(abs(cells), axis=1).T @ width
            for nodes, side in zip(ids, (-1,1), strict=True):
                midpoint_force[nodes,:2] += side*m
                nodal[nodes,:2] += e
            tiles += len(cells)
    action = A @ midpoint_force.ravel()[red.free]
    nodal_action = abs(A) @ nodal.ravel()[red.free]
    np.savez_compressed(dest, mid_action_m=action, nodal_action_m=nodal_action,
                        remainder_m=remainder, component_radius_m=component, center_direction_support_m=along, tiles=tiles)
    return action, nodal_action, remainder, component, along, tiles


def run(sign, out, mode):
    source = inspect.getsource(chain.run)
    key = '        upper = (fc[1]+fp[1])'
    assert source.count(key) == 1
    source = source.replace(key, '''        if idx+1 >= 831:
            us = data["u"][p13["step"]-1:idx+1]
            row.update(nodal_center=center_only,nodal_full_history=inv)
            dest = out/f"{sign}-{idx+1}-direct.npz"
            center_only, inv, tiles = p15_face.integrate(red,refined,center_def,us,kn,kt,mu,
                f,pf,pe,guard,center_only,inv,dest)
            row["integration_tiles"] = tiles
            if idx+1 >= 832:
                row["p15_nodal_past"] = inv.copy()
                A, ix = p15_face.fixed_inverse(red,kn,kt)
                hist = [p13["certified_radius_m"]] + [r["certified_radius_m"] for r in steps]
                with np.load(dest) as saved:
                    center_action = saved["center_after_inverse_m"]
                c_norm = np.linalg.norm(center_action)
                direction = center_action/c_norm if c_norm else np.zeros_like(center_action)
                action, nodal, rem, component, along, ntiles = paired(red,refined,us,hist,(pl,ph),kt,km,A,ix,
                    direction,out/f"{sign}-{idx+1}-pair.npz")
                base = center_action + action
                roundoff = center_only["force_guard_m"] + np.linalg.norm(A)*1e-12*fp[1]
                nominal_rem = center_only["direct_center_remainder_m"]
                q = inv["contraction_upper"]
                scalar = np.linalg.norm(base)+rem
                vector = np.linalg.norm(abs(base)+component)
                error_norm = min(rem,np.linalg.norm(component))
                # ||base+e||^2 includes the signed cross term. Bound its projection
                # in the nominal center direction directly inside each surface tile.
                projected = np.sqrt(np.linalg.norm(base)**2 + error_norm**2 +
                    2*(c_norm*along + np.linalg.norm(action)*error_norm))
                upper_norm = dict(scalar=scalar,component=min(scalar,vector),
                                  projected=min(scalar,vector,projected))[mode]
                eta = upper_norm*(1+1e-9)+nominal_rem+roundoff
                image = eta + q*R
                ok = q < 1 and image < R
                inv = dict(center_after_inverse_upper_m=eta,ball_image_upper_m=image,
                    contraction_upper=q,included=bool(ok),root_radius_upper_m=eta/(1-q)*(1+1e-12) if ok else None,
                    tiles=ntiles,paired_remainder_m=rem,scalar_image_m=scalar*(1+1e-9)+nominal_rem+roundoff+q*R,
                    projected_image_m=projected*(1+1e-9)+nominal_rem+roundoff+q*R,
                    component_image_m=vector*(1+1e-9)+nominal_rem+roundoff+q*R,past_midpoint_norm_m=np.linalg.norm(action),
                    force_guard_m=roundoff, same_paired_integral_nodal_image_m=
                    np.linalg.norm(abs(base)+nodal)*(1+1e-9)+nominal_rem+roundoff+q*R)
''' + key)
    scope = dict(vars(chain), p15_face=p15_face, paired=paired, mode=mode)
    exec(compile(source, str(HERE/"paired_face.py"), "exec"), scope)
    scope["run"](sign, out)


if __name__ == "__main__":
    run(sys.argv[1], Path(sys.argv[2]), sys.argv[3] if len(sys.argv)>3 else "projected")
