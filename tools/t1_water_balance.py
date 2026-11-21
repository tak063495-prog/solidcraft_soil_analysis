"""Audit accepted ONE_WAY states and instantaneous construction transfers in m³."""
import argparse
import json
import math
from collections import defaultdict
from pathlib import Path

LIMITS = {"verification": (1e-7, 1e-8, 1e-12), "practical": (1e-3, 5e-4, 1e-9)}


class Account:
    def __init__(self, initial):
        self.initial = initial
        self.incoming = []
        self.outgoing = []
        self.denominators = []
        self.defects = []

    def record(self, before, after, incoming, outgoing, profile):
        if not all(map(math.isfinite, (before, after, incoming, outgoing))) or min(incoming, outgoing) < 0:
            raise ValueError("invalid inventory or transfer")
        change = after - before
        defect = math.fsum([change, outgoing, -incoming])
        denominator = (incoming + outgoing + abs(change)) / 2
        self.incoming.append(incoming)
        self.outgoing.append(outgoing)
        self.denominators.append(denominator)
        self.defects.append(abs(defect))
        cumulative = math.fsum([after, -self.initial, *self.outgoing, *(-x for x in self.incoming)])
        unsigned = math.fsum(self.defects)
        total_d = math.fsum(self.denominators)
        step_r, cum_r, floor = LIMITS[profile]
        limit = max(floor, step_r * denominator)
        cum_limit = max(floor, cum_r * total_d)
        return dict(delta_m=change, P=incoming, N=outgoing, e=defect, D=denominator,
                    E=cumulative, U=unsigned, sum_D=total_d, step_limit=limit,
                    cumulative_limit=cum_limit, ratio=abs(defect)/denominator if denominator else None,
                    passed=abs(defect) <= limit and max(abs(cumulative), unsigned) <= cum_limit)


def split(values):
    return (math.fsum(-v for v in values if v < 0), math.fsum(v for v in values if v > 0))


class MissingQuantity(ValueError):
    """An active account lacks a quantity needed to judge it."""


def required(mapping, key):
    value = mapping.get(key) if mapping is not None else None
    if value is None:
        raise MissingQuantity(str(key))
    return value


def inventory(state, expected_points=None):
    cells = defaultdict(list)
    volumes = defaultdict(list)
    local_points = defaultdict(set)
    for point in required(required(state, "fluid"), "hydraulic_points"):
        i = required(point, "element_index")
        if expected_points is not None:
            local = required(point, "local_index")
            if local in local_points[i]:
                raise ValueError("duplicate hydraulic integration point")
            local_points[i].add(local)
        w = required(point, "weight_m3")
        cells[i].append(w * required(point, "water_equivalent_mass_per_volume"))
        volumes[i].append(w)
    hydraulic = required(state, "hydraulic")
    heads = required(hydraulic, "relative_cell_head_m")
    storage = required(hydraulic, "cell_storage_offset_m3")
    active = {i for i, h in enumerate(heads) if h is not None}
    if len(storage) != len(heads) or any(storage[i] is None for i in active):
        raise MissingQuantity("active cell storage")
    if expected_points is not None and any(local_points[i] != expected_points for i in active):
        raise MissingQuantity("incomplete hydraulic quadrature")
    if cells.keys() != active:
        raise MissingQuantity("required hydraulic integration points missing")
    return ({i: math.fsum(v) for i, v in cells.items()},
            {i: math.fsum(v) for i, v in volumes.items()})


def audit(mesh, report, profile):
    if profile not in LIMITS:
        raise ValueError("unknown water-balance profile")
    try:
        if "stages" in report.get("case", {}):
            raise MissingQuantity("steady or prescribed pressure has no elapsed flow volume or domain inventory")
        return _audit(mesh, report, profile)
    except MissingQuantity as error:
        return dict(profile=profile, unit="m3", complete=report["complete"],
                    interface_storage_m3=None, passed=None, rows=[], missing=str(error))


def _audit(mesh, report, profile):
    case = report["case"].get("one_way", report["case"])
    # An omitted optional source in the input schema denotes a physical zero.
    sources = {m["id"]: m.get("volume_source_per_s", 0.)
               for m in required(required(case, "flow"), "materials")}
    construction = "initialization" in report
    if construction:
        prior = required(report["initialization"], "last_committed")
        frames = required(report, "frames")
    else:
        prior = required(report, "seed")
        frames = required(report, "frames")
    rows, accounts, generations, events = [], {}, defaultdict(int), set()

    def add(kind, key, start, end, p, n):
        label = f"{kind}/{key}"
        if label not in accounts:
            accounts[label] = Account(start)
        value = accounts[label].record(start, end, p, n, profile)
        rows.append(dict(kind=kind, id=label, time_s=time, **value))

    expected_points = set(range((mesh.get("quadrature_order") or 3)**3)) if "tetrahedra" in mesh else None
    for frame in frames:
        state = required(frame, "state")
        hydraulic = required(state, "hydraulic")
        time = required(hydraulic, "time_s")
        dt = time - required(required(prior, "hydraulic"), "time_s")
        if dt < 0:
            raise ValueError("accepted time reversed")
        old, _ = inventory(prior, expected_points)
        new, volumes = inventory(state, expected_points)
        if construction and "hydraulic_cell_element_indices" in frame:
            ids = required(frame, "hydraulic_cell_element_indices")
            if set(ids) != new.keys():
                raise MissingQuantity("active hydraulic cell head or inventory")
        elif not construction and "tetrahedra" in mesh:
            structure_names = set((case["mechanics"].get("joint_system") or {}).get("structure_bodies", []))
            excluded = {r["id"] for r in mesh["regions"] if r["body"] in structure_names}
            if new.keys() != {i for i, rid in enumerate(mesh["region_ids"]) if rid not in excluded}:
                raise MissingQuantity("active hydraulic cell head or inventory")
        event = frame.get("phase") in ("excavation_release", "material_transition")
        if dt == 0 and (not event or frame["name"] in events):
            if old != new:
                raise ValueError("unbooked inventory change at fixed time")
            prior = state
            continue
        inward, outward = [], []
        if event:
            if dt != 0:
                raise ValueError("event advanced physical time")
            events.add(frame["name"])
            records = [r for r in report["excavations"] + report.get("transitions", [])
                       if r["name"] == frame["name"] and r["time_s"] == time]
            if not records:
                raise MissingQuantity("construction transfer record")
            if len(records) != 1:
                raise ValueError("duplicate construction transfer record")
            record = records[0]
            moves = {}
            if "removal" in record:
                for cell in required(record["removal"], "removed_cells"):
                    if cell["element_index"] in moves:
                        raise ValueError("duplicate removed cell")
                    moves[cell["element_index"]] = -required(cell, "water_equivalent_volume_m3")
                if moves.keys() != old.keys() - new.keys():
                    raise MissingQuantity("removed cell inventory")
                p, n = 0., required(record["removal"], "water_removed_m3")
            else:
                water = required(record, "water")
                internal = water.get("internal_transfers", [])
                for cell in required(water, "changed_cells"):
                    index = cell["element_index"]
                    if index in moves:
                        raise ValueError("duplicate changed cell")
                    previous = cell.get("before")
                    if index in old and previous is None:
                        raise MissingQuantity(f"cell {index} pre-event inventory")
                    if index not in old and previous is not None:
                        raise ValueError("inactive cell has live pre-event inventory")
                    moves[cell["element_index"]] = (required(required(cell, "after"), "water_equivalent_volume_m3")
                        - (required(previous, "water_equivalent_volume_m3") if previous is not None else 0.))
                if not (new.keys() - old.keys()) <= moves.keys():
                    raise MissingQuantity("new cell inventory")
                p, n = required(water, "water_added_m3"), required(water, "water_removed_m3")
                if internal:
                    if p or n or record.get("water_source") is not None:
                        raise ValueError("internal and external supply counted together")
                    net = defaultdict(list)
                    receivers = set()
                    for transfer in internal:
                        donor, receiver = transfer["donor_element"], transfer["receiver_element"]
                        amount = required(transfer, "volume_m3")
                        if (donor not in old or receiver in old or receiver not in new or
                                receiver in receivers or amount < 0 or not math.isfinite(amount)):
                            raise ValueError("invalid internal construction transfer")
                        receivers.add(receiver)
                        net[donor].append(-amount)
                        net[receiver].append(amount)
                    if receivers != new.keys() - old.keys() or net.keys() != moves.keys():
                        raise ValueError("incomplete construction transfer pair")
                    moves = {i: math.fsum(v) for i, v in net.items()}
            inward.append(p)
            outward.append(n)
            for i in sorted(old.keys() | new.keys()):
                if i not in old:
                    generations[i] += 1
                v = moves.get(i, 0.)
                m0, m1 = old.get(i, 0.), new.get(i, 0.)
                add("cell", (i,generations[i]), m0, m1, max(v,0.), max(-v,0.))
                if i in moves:
                    add("event_pair", (frame["name"],i,generations[i]), m0, m1, max(v,0.), max(-v,0.))
        else:
            active = sorted(new)
            if old.keys() != new.keys():
                raise ValueError("unbooked activity change")
            faces = required(frame if construction else report, "hydraulic_faces")
            ids = state["mechanical"].get("element_material_ids") or mesh["material_ids"]
            for i in active:
                source = required(sources, ids[i]) * volumes[i] * dt
                inward.append(max(source, 0.))
                outward.append(max(-source, 0.))
                all_outflow = required(hydraulic, "cell_outflow_m3_s")
                q_face = all_outflow[i] if i < len(all_outflow) else None
                if q_face is None or len(q_face) != 4 or any(q is None for q in q_face):
                    raise MissingQuantity(f"cell {i} outward flow")
                flows = [q * dt for q in q_face] + [-source]
                add("cell", (i,generations[i]), old[i], new[i], *split(flows))
            wall_faces = {}
            for point in state["fluid"].get("joint_traces", []):
                for side in ("side_a", "side_b"):
                    trace = point[side]
                    if trace is not None:
                        wall_faces[trace["face"]] = (point["joint_id"],point["face_index"])
            for fi, face in enumerate(faces):
                owners = [(active[c], j) for c, j in face["owners"]]
                values = [hydraulic["cell_outflow_m3_s"][i][j] * dt for i,j in owners]
                if len(values) == 2:
                    key = tuple((i,j,generations[i]) for i,j in owners)
                    add("face_pair", key, 0., 0., *split(values))
                elif len(values) == 1:
                    p, n = split(values)
                    inward.append(p)
                    outward.append(n)
                    if fi in wall_faces:
                        i, j = owners[0]
                        definition = next(d for d in case["mechanics"]["joint_system"]["joints"]
                                          if d["id"] == wall_faces[fi][0])
                        if definition["water"]["kind"] == "soil":
                            add("impermeable_wall",(i,j,generations[i]),0.,0.,p,n)
                        else:
                            rows.append(dict(kind="reservoir",id=[i,j,generations[i]],time_s=time,P=p,N=n))
                else:
                    raise ValueError("invalid RT0 owner count")
        add("global", "domain", math.fsum(old.values()), math.fsum(new.values()),
            math.fsum(inward), math.fsum(outward))
        prior = state
    judgments = [r["passed"] for r in rows if "passed" in r]
    return dict(profile=profile, unit="m3", complete=report["complete"], interface_storage_m3=None,
                passed=bool(report["complete"] and judgments and all(judgments)), rows=rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mesh", type=Path)
    parser.add_argument("report", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--profile", choices=LIMITS, required=True)
    args = parser.parse_args()
    load = lambda path: json.loads(path.read_text(encoding="utf-8-sig"))
    result = audit(load(args.mesh), load(args.report), args.profile)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
