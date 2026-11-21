"""A stepped pair tied to the top and bottom faces of the same dry Tet10 patch."""

import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import fem_reference as ref


def attached():
    mesh = json.loads((ROOT / "examples/plate6-offset-mesh.json").read_bytes())
    xyz = np.array(mesh["nodes"])
    topology = ref.topology(mesh)
    groups = ["lower", "upper"]
    faces = []
    owners = {}
    nodes = []
    wings = []
    for height, group in enumerate(groups):
        for tri, cells in topology.items():
            if len(cells) != 1 or not np.all(xyz[list(tri), 2] == height):
                continue
            cell = cells[0]
            ids = ref.face_nodes(mesh, cell, list(tri)).tolist()
            a, b, c = xyz[ids[:3]]
            if np.cross(b - a, c - a)[2] * (2 * height - 1) < 0:
                ids = [ids[k] for k in (0, 2, 1, 5, 4, 3)]
            for k in ids:
                if k not in nodes:
                    nodes.append(k)
            fid = f"face-{len(faces)}"
            owners[fid] = dict(cell=int(cell), nodes=ids)
            faces.append(
                dict(
                    id=fid,
                    group=group,
                    offset_m=-0.05,
                    nodes=[nodes.index(k) for k in ids],
                )
            )
        edge = sorted(
            (nodes.index(k) for k in nodes if xyz[k, 2] == height and xyz[k, 1] == 0),
            key=lambda n: xyz[nodes[n], 0],
        )
        wings.append(dict(group=group, edge_nodes=edge))
    total = 6 * len(nodes) + 3 * (len(xyz) - len(nodes))
    # Fix the soil and shell away from the paired edge. All line variables
    # can respond to the applied partial line traction.
    free = [6 * n + k for w in wings for n in w["edge_nodes"] for k in range(6)]
    q = np.zeros(total)
    q[: 6 * len(nodes)].reshape(-1, 6)[:, :3] = [0.001, 0.0, 0.0]
    q[6 * len(nodes) :].reshape(-1, 3)[:] = [0.001, 0.0, 0.0]

    def load(free):
        return dict(
            event="load",
            free_dofs=free,
            prescribed_configuration=q.tolist(),
            nodal_load=[0.0] * total,
            plate_traction_kpa=[0.0] * 3,
        )

    install = dict(event="install", holder="soil", contacts={})
    shell = dict(
        nodes=xyz[nodes].tolist(),
        faces=faces,
        directors=[],
        section=dict(
            id="liner",
            young_kpa=2e5,
            poisson=0.25,
            thickness_m=0.1,
            unit_weight_kn_m3=0.0,
            drilling_factor=1e-4,
        ),
        soil=dict(
            mesh=mesh,
            owners=owners,
            materials={
                str(m): dict(young_kpa=12000.0 * m, poisson=0.25)
                for m in mesh["material_ids"]
            },
            exposed_side={f: 1 for f in owners},
            tie_at_physical_surface=True,
        ),
        steps=[
            load([]),
            install,
            load(free),
            dict(event="remove", holder="soil"),
            load([]),
            install,
            load(free),
        ],
    )
    entry = dict(
        id="working",
        wing="upper",
        interval_m=[0.25, 1.0],
        force_kn_m=[0.1, 0.25, -0.5],
        couple_kn_m_per_m=[0.0, 0.0, 0.0],
    )
    return dict(
        shell=shell,
        joint=dict(
            id="step",
            wings=wings,
            rotation=dict(kind="spring", stiffness_kn_m_per_rad_m=20.0),
            line_m=[[i / 2, 0.0, 0.5] for i in range(3)],
            separators=[
                dict(
                    negative="lower",
                    positive="upper",
                    point_m=[0.0, 0.0, 0.5],
                    normal=[0.0, 0.0, 1.0],
                )
            ],
            loads={3: [entry], 7: [entry]},
        ),
    ), nodes
