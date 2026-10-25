"""External Gmsh preprocessing process; no Gmsh linking in the Rust solver."""
import argparse
import hashlib
import json
from pathlib import Path

import gmsh
import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("step", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--mesh-size-m", type=float, default=2.5)
    parser.add_argument("--crs", default="LOCAL")
    args = parser.parse_args()
    if not np.isfinite(args.mesh_size_m) or args.mesh_size_m <= 0:
        parser.error("mesh size must be finite and positive")
    if args.output.exists():
        parser.error("output already exists")
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.logger.start()
        gmsh.option.setNumber("General.NumThreads", 1)
        # STEP exports from SolidCraft are mm. Explicitly import in mm then convert all nodes.
        gmsh.option.setString("Geometry.OCCTargetUnit", "MM")
        gmsh.model.occ.importShapes(str(args.step.resolve()))
        gmsh.model.occ.synchronize()
        if len(gmsh.model.getEntities(3)) != 1:
            raise ValueError("bootstrap requires exactly one closed volume; body/material mapping is planned")
        h = args.mesh_size_m * 1000.0
        gmsh.option.setNumber("Mesh.MeshSizeMin", h)
        gmsh.option.setNumber("Mesh.MeshSizeMax", h)
        gmsh.option.setNumber("Mesh.ElementOrder", 1)
        gmsh.model.mesh.generate(3)
        tags, coords, _ = gmsh.model.mesh.getNodes()
        coordinates = np.asarray(coords).reshape(-1, 3) * 0.001
        by_tag = {int(tag): coordinates[i].tolist() for i, tag in enumerate(tags)}
        kinds, _, blocks = gmsh.model.mesh.getElements(3)
        tets = []
        for kind, block in zip(kinds, blocks):
            if kind != 4:
                raise ValueError(f"unsupported volume element type {kind}, expected linear tetrahedron")
            tets.extend(np.asarray(block).reshape(-1, 4).tolist())
        if not tets:
            raise ValueError("no volume mesh")
        used = sorted({int(tag) for tet in tets for tag in tet})
        index = {tag: i for i, tag in enumerate(used)}
        nodes = [by_tag[tag] for tag in used]
        tet_indices = [[index[int(tag)] for tag in tet] for tet in tets]
        for tet in tet_indices:
            p = np.asarray([nodes[i] for i in tet])
            det = np.linalg.det((p[1:] - p[0]).T)
            if det < 0:
                tet[1], tet[2] = tet[2], tet[1]
        minimum_z = min(p[2] for p in nodes)
        tolerance = max(np.ptp(np.asarray(nodes), axis=0).max(), 1.0) * 1e-9
        bottom = [i for i, p in enumerate(nodes) if abs(p[2] - minimum_z) < tolerance]
        if len(bottom) < 3:
            raise ValueError("no flat bottom node set: configure physical boundaries before solving")
        mesh = {
            "schema": "soil-fem-mesh/1", "units": "m-kN-kPa", "axis": "z_up",
            "crs": args.crs, "origin_m": [0.0, 0.0, 0.0], "nodes": nodes,
            "tetrahedra": tet_indices, "material_ids": [1] * len(tets),
            "node_sets": {"bottom": bottom},
        }
        with args.output.open("x", encoding="utf-8") as f:
            json.dump(mesh, f, allow_nan=False)
        gmsh.write(str(args.output.with_suffix(".msh")))
        provenance = {"source": str(args.step.resolve()),
                      "source_sha256": hashlib.sha256(args.step.read_bytes()).hexdigest(),
                      "gmsh": gmsh.__version__, "step_import_unit": "mm", "solver_unit": "m",
                      "mesh_size_m": args.mesh_size_m, "nodes": len(nodes), "tetrahedra": len(tets)}
        args.output.with_suffix(".provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
        print(f"Gmsh {gmsh.__version__}: {len(nodes)} nodes, {len(tets)} tetrahedra; mm -> m")
    finally:
        args.output.with_suffix(".gmsh.log").write_text("\n".join(gmsh.logger.get()), encoding="utf-8")
        gmsh.logger.stop()
        gmsh.finalize()


if __name__ == "__main__":
    main()
