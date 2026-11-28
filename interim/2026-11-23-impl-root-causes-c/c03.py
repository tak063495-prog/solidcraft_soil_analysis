"""Join same-time classification disagreements to material, contact and nodal fields."""
import argparse
import copy
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "2026-12-18-c03-plastic-bulk-cause-diagnosis"))
from c03diag_metrics import CATEGORIES, load, metric
from c03diag_study import Model, read, write
from c03diag_reference import CONTACT, RETURNS, STAGES


def verify_axis(a, b, axis):
    left, right = copy.deepcopy(a), copy.deepcopy(b)
    def controls(case):
        return {"material": case["settings"].pop("material_maximum_substep"),
                "contact": case["joint_system"].pop("maximum_increment_m")}
    first, second = controls(left), controls(right)
    changed = {key for key in first if first[key] != second[key]}
    expected = set() if axis == "solver" else {axis}
    if left != right or changed != expected:
        raise ValueError("inputs differ outside the selected axis, or the axis did not change")
    return [first, second]


def mismatch_rows(rows, left, right, mesh, materials, contacts):
    """Interpolate u at its own GP. Stress and traction stay on their native supports."""
    records = []
    elements = mesh.get("tetrahedra10") or mesh["tetrahedra"]
    for label, names in (("contact_status", CONTACT), ("material_return", RETURNS)):
        for values in (left[label], right[label]):
            if not np.issubdtype(values.dtype, np.integer) or np.any((values < 0) | (values >= len(names))):
                raise ValueError("unknown classification code")
        for t, gp in zip(*np.nonzero(left[label] != right[label]), strict=True):
            a, b = int(left[label][t,gp]), int(right[label][t,gp])
            du = left["u"][t] - right["u"][t]
            ds = left["stress"][t] - right["stress"][t]
            dq = left["traction"][t] - right["traction"][t]
            record = dict(stage=rows[t]["stage"], progress=rows[t]["progress"], field=label,
                          gp=int(gp), left=names[a], right=names[b])
            if label == "contact_status":
                p = contacts[gp]
                owners = (p["owner_a"], p["owner_b"])
                record.update(joint=p["joint_id"], face=p["face_index"], position_m=p["position_m"],
                              delta_traction_kpa=dq[gp].tolist(),
                              delta_u_a_m=(p["shape"] @ du[p["nodes_a"]]).tolist(),
                              delta_u_b_m=(p["shape"] @ du[p["nodes_b"]]).tolist(),
                              adjacent_stress=[dict(gp=i, element=m["element_index"],
                                                    local=m["local_index"], delta_kpa=ds[i].tolist())
                                               for i,m in enumerate(materials) if m["element_index"] in owners])
            else:
                p = materials[gp]
                e = p["element_index"]
                record.update(element=e, local=p["local_index"], position_m=p["position_m"],
                              delta_stress_kpa=ds[gp].tolist(),
                              delta_u_m=(p["shape_functions"] @ du[elements[e]]).tolist(),
                              adjacent_traction=[dict(gp=i, joint=q["joint_id"], face=q["face_index"],
                                                      delta_kpa=dq[i].tolist()) for i,q in enumerate(contacts)
                                                 if e in (q["owner_a"], q["owner_b"])])
            records.append(record)
    return records


def equilibrium(folder, expected):
    native = folder / "schedule.json"
    independent = folder / "result.json"
    if native.is_file() == independent.is_file():
        raise ValueError("one solver's accepted schedule is required")
    if native.is_file():
        schedule = read(native)
        if [(r["stage"], r["accepted"]) for r in schedule] != [
                (stage, sum(key[0] == stage for key in expected)) for stage in STAGES]:
            raise ValueError("native schedule differs from field history")
        return "native", max(r["maximum_force"] for r in schedule)
    result = read(independent)
    accepted = [(r["stage"],r["progress"]) for r in result["rows"]]
    anchors = set(expected)
    if not result["complete"] or [pair for pair in accepted if pair in anchors] != expected:
        raise ValueError("independent route incomplete")
    return "independent", max(r["relative_force_residual"] for r in result["rows"])


def measure(matrix, left_dir, right_dir, native_map, axis):
    a, b = [read(folder / "stages.json") for folder in (left_dir, right_dir)]
    controls = verify_axis(a, b, axis)
    mesh, ground = [read(matrix / f"{name}.json") for name in ("mesh", "ground")]
    if len(mesh["tetrahedra"]) != 12:
        raise ValueError("only the 12-element case is supported")
    n = round(1 / a["settings"]["initial_increment"])
    expected = [(s, i/n) for s in STAGES for i in range(1,n+1)]
    solvers, forces = zip(*(equilibrium(folder, expected) for folder in (left_dir, right_dir)), strict=True)
    if (solvers[0] != solvers[1]) != (axis == "solver"):
        raise ValueError("solver identity inconsistent with comparison axis")
    (rows, left), (other, right) = [load(folder / "fields.jsonl") for folder in (left_dir, right_dir)]
    for history in (rows, other):
        if [(r["stage"], r["progress"]) for r in history] != expected:
            raise ValueError("missing, duplicate or unordered exact accepted times")
    if left.keys() != right.keys() or any(left[k].shape != right[k].shape for k in left):
        raise ValueError("full-field inventory or dimensions differ")
    raw = read(native_map)
    materials, contacts = raw["integration_points"], raw["joint_integration_points"]
    model = Model(mesh, ground, a)
    if raw["ground"] != ground or not verify_map_case(raw["case"], a):
        raise ValueError("GP descriptors from a different physical case")
    if (len(materials) != left["stress"].shape[1] or len(contacts) != left["traction"].shape[1]
            or left["u"].shape[1] != len(mesh["nodes"])):
        raise ValueError("GP or node count differs")
    np.testing.assert_allclose([p["position_m"] for p in materials], model.positions, rtol=0, atol=1e-12)
    expected_materials = [(element, local, shape) for element, (op, _) in enumerate(model.ops)
                          for local, shape in enumerate(op["shape"])]
    for p, (e, local, shape) in zip(materials, expected_materials, strict=True):
        if (p["element_index"], p["local_index"]) != (e, local):
            raise ValueError("material GP ownership differs")
        np.testing.assert_allclose(p["shape_functions"], shape, rtol=0, atol=1e-12)
    for p,q in zip(contacts, model.system.points, strict=True):
        if any(p[k] != q[k] for k in ("joint_id","face_index","owner_a","owner_b","nodes_a","nodes_b")):
            raise ValueError("contact ownership differs")
        for k in ("position_m","shape","weight_m2","normal"):
            np.testing.assert_allclose(p[k], q[k], rtol=0, atol=1e-12)
    mismatches = mismatch_rows(rows, left, right, mesh, materials, contacts)
    fields = {k: (dict(mismatches=int(np.count_nonzero(left[k] != right[k]))) if k in CATEGORIES
                  else metric(left[k],right[k])) for k in sorted(left)}
    return dict(axis=axis, left=left_dir.name, right=right_dir.name, controls=controls,
                accepted=len(rows), maximum_force_residual=forces, mesh_nodes=len(mesh["nodes"]), full_fields=fields,
                mismatches=mismatches, mcc_plastic_work=None, gci=None, passed=False)


def verify_map_case(raw, case):
    def without_sizes(value):
        result = copy.deepcopy(value)
        result["settings"].pop("material_maximum_substep")
        result["joint_system"].pop("maximum_increment_m")
        return result
    return without_sizes(raw) == without_sizes(case)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matrix", type=Path)
    parser.add_argument("left")
    parser.add_argument("right")
    parser.add_argument("--native-map", type=Path, required=True)
    parser.add_argument("--axis", choices=("material", "contact", "solver"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    write(args.output, measure(args.matrix, args.matrix/args.left, args.matrix/args.right,
                               args.native_map, args.axis))
