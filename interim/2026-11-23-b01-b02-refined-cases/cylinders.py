"""Build the approved SSLS104A curved-shell inputs in m, kN and kPa."""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "2026-11-21-impl-t5-curved-tri6-shell"))
import standard

GRID_SHA256 = "609a8d0944fe8df99103650bf85a987c32fd199af61334627d444cc4188e9fcd"


def grid():
    raw = (HERE / "ssls104a-grid.json").read_bytes()
    if hashlib.sha256(raw).hexdigest() != GRID_SHA256:
        raise ValueError("SSLS104A input grid fingerprint changed")
    return json.loads(raw)


def edf_geometry():
    data = grid()
    points = [np.array([x, 300 - z, y]) for x, y, z in data["vertices"]]
    directors = [np.array([x, 0, z]) / np.hypot(x, z) for x, _, z in points]
    edges = {}
    faces = []
    for a, b, c in data["triangles"]:
        if np.cross(points[b] - points[a], points[c] - points[a]) @ directors[a] < 0:
            b, c = c, b
        face = [a, b, c]
        for i, j in ((a, b), (b, c), (c, a)):
            key = tuple(sorted((i, j)))
            if key not in edges:
                direction = directors[i] + directors[j]
                direction /= np.linalg.norm(direction)
                edges[key] = len(points)
                directors.append(direction)
                points.append(300 * direction + [0, (points[i][1] + points[j][1]) / 2, 0])
            face.append(edges[key])
        faces.append(dict(id=f"f{len(faces):04}", group="shell", nodes=face, offset_m=0.0))
    return np.array(points).tolist(), np.array(directors).tolist(), faces


def cylinder_input(case_id):
    card = next(c for c in json.loads(standard.CARDS.read_bytes())["benchmarks"]
                if c["id"] == "pinched_cylinder")
    case = standard.build(card, 16, False)
    if case_id == "b01-regular16":
        free = set(case["steps"][1]["free_dofs"])
        # EDF excludes endpoints from the generator groups. Remove the two
        # additional drilling restraints from the T5 builder at B and A.
        corner_b = 32 * 33
        free.update((6 * corner_b + 5, 6 * (corner_b + 32) + 3))
    elif case_id == "b01-edf364":
        points, directions, faces = edf_geometry()
        case["faces"] = faces
        case["nodes"] = points
        case["directors"] = [dict(group="shell", node=i, direction=v)
                             for i, v in enumerate(directions)]
        fixed = set()
        for i, (x, y, z) in enumerate(points):
            if y == 300:
                axes = [0, 2, 4]
            else:
                axes = [1, 3, 5] if y == 0 else []
                if x == 0:
                    axes.extend([0, 4, 5])
                if z == 0:
                    axes.extend([2, 3, 4])
            fixed.update(6 * i + j for j in axes)
        free = set(range(6 * len(points))) - fixed
    else:
        raise ValueError(f"unknown cylinder case: {case_id}")
    probe, = [i for i, p in enumerate(case["nodes"]) if p == [0, 0, 300]]
    size = 6 * len(case["nodes"])
    load = [0.0] * size
    load[6 * probe + 2] = -0.00025
    case["steps"][1].update(free_dofs=sorted(free), prescribed_configuration=[0.0] * size,
                            nodal_load=load)
    case["section"]["young_kpa"] = 3000.0
    return case, probe
