import argparse
import json
from pathlib import Path
import meshio
import numpy as np


def export(mesh_path, state_path, output):
    if output.exists():
        raise ValueError("VTU output already exists")
    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    state = json.loads(state_path.read_text(encoding="utf-8"))
    if mesh.get('tetrahedra10') or state.get('schema')=='soil-fem-initial-state/2':
        from high_order_io import initial
        return initial(mesh,state,output)
    points = state["point_states"]
    if len(points) != len(mesh["tetrahedra"]):
        raise ValueError("mesh/state length mismatch")
    cell = {"material_id":[np.asarray(mesh["material_ids"])],"region_id":[np.asarray(mesh["region_ids"])]}
    for name in ("pore_pressure_kpa","porosity","mcc_pc_over_mean_effective","normalized_yield",
                 "bulk_modulus_kpa","shear_modulus_kpa","total_stress_compression_positive_kpa"):
        cell[name] = [np.asarray([p[name] for p in points])]
    for name in ("void_ratio","preconsolidation_pressure_kpa","effective_stress_compression_positive_kpa"):
        cell[name] = [np.asarray([p["mcc_state"][name] for p in points])]
    meshio.write_points_cells(output,mesh["nodes"],[("tetra",mesh["tetrahedra"])],cell_data=cell)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mesh",type=Path)
    parser.add_argument("state",type=Path)
    parser.add_argument("output",type=Path)
    a = parser.parse_args()
    export(a.mesh,a.state,a.output)
