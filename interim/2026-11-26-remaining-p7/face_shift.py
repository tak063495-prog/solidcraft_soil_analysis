"""Translate both P6 marking boundaries, with the number of splits fixed per face and depth."""
import hashlib
import inspect
import json
from math import factorial, fsum
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
P6 = HERE.parent / "2026-11-26-remaining-p6"
sys.path.insert(0, str(P6))
import face as p6

SHIFT = np.array([1/16, 1/32])  # physical surface coordinates, metres


def roots():
    for i in range(8):
        for j in range(8-i):
            yield np.array([(i,j), (i+1,j), (i,j+1)]) / 8
            if i+j+1 < 8:
                yield np.array([(i+1,j), (i+1,j+1), (i,j+1)]) / 8


def key(tile):
    return tuple(sorted(map(tuple, tile)))


def marked_leaves(old):
    leaves = {key(t) for t in old["leaves"]}
    corners = np.array(old["corners"])[:, :2]

    def center(t):
        x, y = t.mean(axis=0)
        return np.array([1-x-y, x, y]) @ corners

    original = list(roots())
    candidates = list(roots())
    new_leaves = []
    marks = []
    for depth in range(2):
        old_split = [t for t in original if key(t) not in leaves]
        targets = np.array([center(t) + SHIFT for t in old_split])
        # Stable ties follow the barycentric traversal of the base triangles.
        def distance(i):
            return np.min(np.sum((targets-center(candidates[i]))**2, axis=1))
        ranked = sorted(range(len(candidates)), key=distance) if len(targets) else []
        selected = set(ranked[:len(targets)])
        picked = [t for i,t in enumerate(candidates) if i in selected]
        marks.append(dict(depth=depth, count=len(picked),
                          old_centers=[center(t).tolist() for t in old_split],
                          new_centers=[center(t).tolist() for t in picked],
                          retained=len({key(t) for t in old_split} & {key(t) for t in picked})))
        new_leaves.extend(t for i,t in enumerate(candidates) if i not in selected)
        candidates = [child for t in picked for child in p6.split(t)]
        original = [child for t in old_split for child in p6.split(t)]
    assert all(key(t) in leaves for t in original)
    new_leaves.extend(candidates)
    assert len(new_leaves) == len(leaves)
    return np.array(new_leaves), marks


def run(out):
    old = json.loads((P6/"c02/focused/tiles.json").read_text())["faces"]
    rules, records = {}, []
    bary, w = p6.face.triangle_rule(1)
    for original in old:
        tiles, marks = marked_leaves(original)
        xy = (bary @ tiles).reshape(-1,2)
        weights = np.outer(abs(np.linalg.det(tiles[:,1:]-tiles[:,:1])), w).ravel()
        assert weights.min() > 0
        for i in range(6):
            for j in range(6-i):
                np.testing.assert_allclose(weights @ (xy[:,0]**i * xy[:,1]**j),
                    2*factorial(i)*factorial(j)/factorial(i+j+2), rtol=0, atol=3e-15)
        rules[original["parent"]] = np.column_stack((1-xy.sum(axis=1),xy)), weights
        records.append(dict(parent=original["parent"], marks=marks, leaves=tiles.tolist()))
    source = inspect.getsource(p6.face.surface)
    for a,b in [('    bary, weights = triangle_rule(divisions)\n',''),
                ('        corners = xyz[left]','        corners = xyz[left]\n        bary, weights = rules[i]'),
                ('assert_allclose(sum(areas)', 'assert_allclose(fsum(areas)')]:
        assert source.count(a) == 1
        source = source.replace(a,b)
    ns = dict(vars(p6.face), rules=rules, fsum=fsum)
    exec(compile(source,p6.face.__file__,"exec"),ns)
    p6.face.surface = ns["surface"]
    p6.driver.run(out, 8)
    cost = json.loads((out/"cost.json").read_text())
    assert cost.pop("subdivisions_per_edge") == 8 and cost["area_points"] == 2660
    cost.update(base_subdivisions=8, shift_m=SHIFT.tolist())
    (out/"cost.json").write_text(json.dumps(cost,indent=2)+"\n")
    (out/"tiles.json").write_text(json.dumps(dict(shift_m=SHIFT.tolist(), faces=records))+"\n")


def score(out):
    paths = {"native":p6.driver.P3, **{k:P6/"c02"/k for k in ("local1","focused","local2")},
             "shifted":out}
    data = {k:np.load(v/"principal.npz") for k,v in paths.items()}
    report = {}
    for left,right in [("native",k) for k in list(paths)[1:]] + [(k,"shifted") for k in ("local1","focused","local2")]:
        a,b = data[left],data[right]
        for name in ("parent","position","weight"):
            np.testing.assert_array_equal(a[name],b[name])
        pair = {}
        for field in ("slip","signed_area_integral"):
            weight = np.sqrt(b["weight"])[None,:,None] if field=="slip" else 1.
            assert a[field].shape == b[field].shape and len(b[field]) == 448
            delta = (a[field]-b[field])*weight
            norm = np.linalg.norm(delta)
            denominator = np.linalg.norm(b[field]*weight)
            pair[field] = dict(numerator=norm, denominator=denominator, relative=norm/denominator,
                               within_width=bool(norm/denominator <= .001),
                               stage_squared_share=(np.square(delta).reshape(7,64,-1).sum(axis=(1,2))/np.square(delta).sum()).tolist())
            if field=="slip":
                xy = b["position"][:,:2]
                edge = np.minimum(xy,1-xy).min(axis=1)<.1
                pair[field]["edge_squared_share"] = np.square(delta[:,edge]).sum()/np.square(delta).sum()
        report[left+"-"+right] = pair
    prior = json.loads((P6/"c02/summary.json").read_text())
    hashes = {k:hashlib.sha256((v/"principal.npz").read_bytes()).hexdigest() for k,v in paths.items()}
    root = HERE.parents[1]
    for k,v in paths.items():
        if k!="shifted":
            assert hashes[k] == prior["sha256"][str((v/"principal.npz").relative_to(root))]
    clock = p6.driver.P3/"increments.npz"
    assert hashlib.sha256(clock.read_bytes()).hexdigest()==prior["sha256"][str(clock.relative_to(root))]
    old_comparisons = {}
    for pair in ("native-local1","native-focused","native-local2"):
        values = {}
        for k,v in report[pair].items():
            old = prior["pairs"][pair][k]["relative"]
            values[k] = dict(archived=old, independently_rescored=v["relative"],
                             bitwise_equal=bool(old==v["relative"]),
                             ulps=float((v["relative"]-old)/np.spacing(old)))
        old_comparisons[pair] = values
    (out/"comparison.json").write_text(json.dumps(dict(width=.001,pairs=report,sha256=hashes,
        historical_norm_check=old_comparisons),indent=2)+"\n")


def compare_repeat(out,repeat):
    old = P6/"c02/focused"
    report = {}
    with np.load(old/"principal.npz") as a,np.load(repeat/"principal.npz") as b:
        for k in ("weight","position","parent"):
            np.testing.assert_array_equal(a[k],b[k])
        for k in ("slip","signed_area_integral"):
            w = np.sqrt(a["weight"])[None,:,None] if k=="slip" else 1.
            report[k] = dict(relative=np.linalg.norm((a[k]-b[k])*w)/np.linalg.norm(a[k]*w),
                             max_absolute=float(abs(a[k]-b[k]).max()))
    # Mark selection on the repeated old case should match triangle-for-triangle.
    assert json.loads((repeat/"tiles.json").read_text())==json.loads((old/"tiles.json").read_text())
    (out/"repeat-comparison.json").write_text(json.dumps(report,indent=2)+"\n")


if __name__ == "__main__":
    out = Path(sys.argv[1]).resolve()
    if "--repeat" in sys.argv:
        compare_repeat(out,Path(sys.argv[-1]))
    else:
        (score if "--score" in sys.argv else run)(out)
