"""Nested local triangle quadrature on the unchanged native clock and P4 volume."""
import inspect
import json
from math import factorial, fsum
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
P5 = HERE.parent / "2026-11-26-remaining-p5"
sys.path.insert(0, str(P5))
# Import the archived reference by path: its own dependencies also use c02_reference.
import importlib.util
spec = importlib.util.spec_from_file_location("p5_face_driver", P5 / "c02_reference.py")
driver = importlib.util.module_from_spec(spec)
spec.loader.exec_module(driver)
face = driver.c02_face


def split(tile):
    a, b, c = tile
    u, v, w = (a+b)/2, (b+c)/2, (c+a)/2
    return [np.array(t) for t in ((a,u,w), (u,b,v), (w,v,c), (u,v,w))]


def run(out, level, focused=False):
    assert level in (1, 2) and (not focused or level == 2)
    native = np.load(driver.P3 / "principal.npz")
    archived = np.load(P5 / "c02/d8/principal.npz")
    positions = native["position"][:, :2]
    parents = native["parent"]
    # Integrated squared field discrepancy at each passive point in final return.
    indicator = ((native["slip"][-64:] - archived["slip"][-64:])**2).sum(axis=(0, 2))
    is_inside = np.minimum(positions, 1-positions).min(axis=1) >= .1
    limits = np.quantile(indicator[is_inside], [.8, .95])
    bary, unit_weights = face.triangle_rule(1)
    edge_cut = np.median(indicator[~is_inside])
    records = []

    def rule(parent, corners):
        leaves, depths = [], []
        def visit(tile, depth):
            center = tile.mean(axis=0)
            point = np.r_[1-center.sum(), center] @ corners[:, :2]
            edge_distance = np.minimum(point, 1-point).min()
            ids = np.flatnonzero(parents == parent)
            near = ids[np.linalg.norm(positions[ids]-point, axis=1).argmin()]
            refine = False
            if depth < level:
                refine = (edge_distance < (.1, .05)[depth] or
                          (is_inside[near] and indicator[near] >= limits[depth]))
                if depth == 1 and focused:
                    refine = ((edge_distance < .025 and indicator[near] >= edge_cut) or
                              (is_inside[near] and indicator[near] >= np.quantile(indicator[is_inside], .99)))
            if refine:
                for child in split(tile):
                    visit(child, depth+1)
            else:
                leaves.append(tile)
                depths.append(depth)
        for i in range(8):
            for j in range(8-i):
                visit(np.array([(i,j), (i+1,j), (i,j+1)])/8, 0)
                if i+j+1 < 8:
                    visit(np.array([(i+1,j), (i+1,j+1), (i,j+1)])/8, 0)
        tiles = np.array(leaves)
        xy = (bary @ tiles).reshape(-1, 2)
        areas = np.abs(np.linalg.det(tiles[:, 1:] - tiles[:, :1]))
        weights = np.outer(areas, unit_weights).ravel()
        assert min(weights) > 0 and np.any(np.array(depths) == 0)
        for i in range(6):
            for j in range(6-i):
                exact = 2*factorial(i)*factorial(j)/factorial(i+j+2)
                np.testing.assert_allclose(weights @ (xy[:,0]**i * xy[:,1]**j), exact,
                                           atol=3e-15, rtol=0)
        records.append(dict(parent=int(parent), tiles_by_depth=np.bincount(depths,minlength=3).tolist(),
                            corners=corners.tolist(), leaves=tiles.tolist(), points=len(weights)))
        return np.column_stack((1-xy.sum(axis=1),xy)), weights

    source = inspect.getsource(face.surface)
    get_rule = "    bary, weights = triangle_rule(divisions)\n"
    corners = "        corners = xyz[left]"
    assert source.count(get_rule) == source.count(corners) == 1
    source = source.replace(get_rule, "").replace(corners, corners+"\n        bary, weights = rule(i, corners)")
    assert source.count("assert_allclose(sum(areas)") == 1
    source = source.replace("assert_allclose(sum(areas)", "assert_allclose(fsum(areas)")
    namespace = dict(vars(face), rule=rule, fsum=fsum)
    exec(compile(source, face.__file__, "exec"), namespace)
    face.surface = namespace["surface"]
    # Reuse D8 invocation for volume degree, time grid and observation only.
    driver.run(out, 8)
    cost = json.loads((out / "cost.json").read_text())
    del cost["subdivisions_per_edge"]
    cost.update(base_subdivisions=8, level=level)
    assert cost["area_points"] == sum(r["points"] for r in records)
    if focused:
        cost["focused_second_pass"] = True
        assert cost["area_points"] < 2 * 16**2 * 7
    (out/"cost.json").write_text(json.dumps(cost, indent=2)+"\n")
    (out/"tiles.json").write_text(json.dumps(dict(return_indicator=indicator.tolist(),
        interior_thresholds=limits.tolist(),faces=records))+"\n")


if __name__ == "__main__":
    run(Path(sys.argv[1]), int(sys.argv[2]), "--focused" in sys.argv[3:])
