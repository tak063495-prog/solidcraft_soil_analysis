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


def inventory(state):
    cells = defaultdict(list)
    volumes = defaultdict(list)
    for point in state["fluid"]["hydraulic_points"]:
        i = point["element_index"]
        w = point["weight_m3"]
        cells[i].append(w * point["water_equivalent_mass_per_volume"])
        volumes[i].append(w)
    active = {i for i, h in enumerate(state["hydraulic"]["relative_cell_head_m"]) if h is not None}
    if cells.keys() != active:
        raise ValueError("required hydraulic integration points missing")
    return ({i: math.fsum(v) for i, v in cells.items()},
            {i: math.fsum(v) for i, v in volumes.items()})


def audit(mesh, report, profile):
    case = report["case"].get("one_way", report["case"])
    sources = {m["id"]: m.get("volume_source_per_s", 0.) for m in case["flow"]["materials"]}
    construction = "initialization" in report
    if construction:
        prior = report["initialization"]["last_committed"]
        frames = report["frames"]
    else:
        prior = report["seed"]
        frames = report["frames"]
    rows, accounts, generations, events = [], {}, defaultdict(int), set()

    def add(kind, key, start, end, p, n):
        label = f"{kind}/{key}"
        if label not in accounts:
            accounts[label] = Account(start)
        value = accounts[label].record(start, end, p, n, profile)
        rows.append(dict(kind=kind, id=label, time_s=time, **value))

    for frame in frames:
        state = frame["state"]
        hydraulic = state["hydraulic"]
        time = hydraulic["time_s"]
        dt = time - prior["hydraulic"]["time_s"]
        if dt < 0:
            raise ValueError("accepted time reversed")
        old, _ = inventory(prior)
        new, volumes = inventory(state)
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
            record = next(r for r in report["excavations"] + report.get("transitions", [])
                          if r["name"] == frame["name"] and r["time_s"] == time)
            moves = {}
            if "removal" in record:
                for cell in record["removal"]["removed_cells"]:
                    moves[cell["element_index"]] = -cell["water_equivalent_volume_m3"]
                p, n = 0., record["removal"]["water_removed_m3"]
            else:
                water = record["water"]
                for cell in water["changed_cells"]:
                    previous = cell["before"]
                    moves[cell["element_index"]] = (cell["after"]["water_equivalent_volume_m3"]
                        - (previous["water_equivalent_volume_m3"] if previous else 0.))
                p, n = water["water_added_m3"], water["water_removed_m3"]
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
            faces = frame["hydraulic_faces"] if construction else report["hydraulic_faces"]
            ids = state["mechanical"].get("element_material_ids") or mesh["material_ids"]
            for i in active:
                source = sources[ids[i]] * volumes[i] * dt
                inward.append(max(source, 0.))
                outward.append(max(-source, 0.))
                flows = [q * dt for q in hydraulic["cell_outflow_m3_s"][i]] + [-source]
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
