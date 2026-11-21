"""Two separately trimmed, affine Tri6 wings and an explicitly owned line."""

import numpy as np


def patch(n=2, law="rigid", alpha=1e-4, angle=90.0, axes=None, distorted=False):
    axes = np.eye(3) if axes is None else np.asarray(axes)
    x = axes[:, 0]
    y = axes[:, 1]
    z = axes[:, 2]
    radians = np.deg2rad(angle)
    # The separating plane is y=0, with both finite thicknesses on their
    # declared sides (including a parallel step when angle=0).
    sides = [y, -np.cos(radians) * y + np.sin(radians) * z]
    origins = [0.14 * y + 0.03 * z, -0.14 * y + 0.12 * z]
    if sides[1] @ y > 0:
        origins[1] -= (np.dot(sides[1], y) + 0.04) * y
    nodes, faces, wings, fixed = [], [], [], []
    length = 1.0
    line = [x * (i / (2 * n)) for i in range(2 * n + 1)]
    for s, (direction, origin) in enumerate(zip(sides, origins)):
        group = f"wing-{s}"
        start = len(nodes)

        def node(i, j):
            return start + (2 * n + 1) * j + i

        # Distort the coarse vertices, then construct true edge midpoints.
        vertices = {}
        for i in range(n + 1):
            for j in range(n + 1):
                p = np.array([i / n, j / n])
                if distorted and 0 < i < n and 0 < j < n:
                    p += np.array([0.13, -0.08]) * np.sin(i + j) / n
                vertices[i, j] = p

        def point(i, j):
            a, b = i // 2, j // 2
            if i % 2 and j % 2:
                return (vertices[a + 1, b] + vertices[a, b + 1]) / 2
            return (vertices[a, b] + vertices[a + i % 2, b + j % 2]) / 2

        for j in range(2 * n + 1):
            for i in range(2 * n + 1):
                xx, yy = point(i, j)
                nodes.append((origin + x * xx + direction * yy).tolist())
                if j == 2 * n:
                    fixed.extend(range(6 * node(i, j), 6 * node(i, j) + 6))
        for j in range(n):
            for i in range(n):
                a, b = 2 * i, 2 * j
                for tri in [
                    [(a, b), (a + 2, b), (a, b + 2)],
                    [(a + 2, b + 2), (a, b + 2), (a + 2, b)],
                ]:
                    ids = tri + [
                        ((u[0] + v[0]) // 2, (u[1] + v[1]) // 2)
                        for u, v in zip(tri, tri[1:] + tri[:1])
                    ]
                    faces.append(
                        dict(
                            id=f"{group}-{len(faces)}",
                            group=group,
                            nodes=[node(*p) for p in ids],
                            offset_m=[0.022, -0.014][s],
                        )
                    )
        wings.append(
            dict(group=group, edge_nodes=[node(i, 0) for i in range(2 * n + 1)])
        )
    total = 6 * len(nodes)
    free = sorted(set(range(total)) - set(fixed))
    section = dict(
        id="dry",
        young_kpa=2e5,
        poisson=0.25,
        thickness_m=0.08,
        unit_weight_kn_m3=0.0,
        drilling_factor=alpha,
    )
    shell = dict(
        nodes=nodes,
        faces=faces,
        directors=[],
        section=section,
        steps=[
            dict(event="install", holder="assembly", contacts={}),
            dict(
                event="load",
                prescribed_configuration=[0.0] * total,
                free_dofs=free,
                nodal_load=[0.0] * total,
                plate_traction_kpa=[0.0] * 3,
            ),
        ],
    )
    loads = [
        dict(
            id="edge-traction",
            wing="wing-0",
            interval_m=[0.25, length],
            force_kn_m=(axes @ np.array([0.2, -0.7, -1.0])).tolist(),
            couple_kn_m_per_m=(axes @ np.array([0.03, 0.01, -0.02])).tolist(),
        )
    ]
    rotation = dict(kind=law)
    if law == "spring":
        rotation["stiffness_kn_m_per_rad_m"] = 25.0
    return dict(
        shell=shell,
        joint=dict(
            id="fold",
            rotation=rotation,
            line_m=np.asarray(line).tolist(),
            wings=wings,
            loads={2: loads},
            separators=[
                dict(
                    negative="wing-1",
                    positive="wing-0",
                    point_m=[0.0] * 3,
                    normal=y.tolist(),
                )
            ],
        ),
    )
