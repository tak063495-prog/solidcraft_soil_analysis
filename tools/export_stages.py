"""Export completed staged FEM states: active-only VTU, PVD and CSV diagnostics."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import meshio
import numpy as np


def export(mesh_file,result_file,output):
    result=json.loads(result_file.read_text(encoding="utf-8"))
    mesh=json.loads(mesh_file.read_text(encoding="utf-8"))
    if mesh.get('tetrahedra10') or result.get('schema') in ('soil-fem-stage-result/3','soil-fem-stage-result/4','soil-fem-stage-result/5','soil-fem-stage-result/6','soil-fem-stage-result/7','soil-fem-stage-result/8','soil-fem-stage-result/9','soil-fem-stage-result/10','soil-fem-stage-result/11','soil-fem-stage-result/12','soil-fem-stage-result/13'):
        from high_order_io import stages
        return stages(mesh_file,result_file,output)
    if result.get("schema")!="soil-fem-stage-result/1" or not result["complete"]:
        raise ValueError("only completed staged runs are exported as accepted results")
    output.mkdir(exist_ok=False,parents=True)
    nodes=np.asarray(mesh["nodes"]);tets=np.asarray(mesh["tetrahedra"])
    collection=ET.Element("VTKFile",type="Collection",version="0.1",byte_order="LittleEndian")
    data=ET.SubElement(collection,"Collection")
    rows=[]
    for i,frame in enumerate(result["frames"]):
        state=frame["state"];active=np.asarray(state["active_elements"],dtype=bool)
        if len(active)!=len(tets) or len(state["displacement_m"])!=len(nodes):
            raise ValueError("mesh/result size mismatch")
        ids=np.flatnonzero(active);used=np.unique(tets[active]);mapping=np.full(len(nodes),-1,dtype=int);mapping[used]=np.arange(len(used))
        u=np.asarray(state["displacement_m"])
        previous=np.asarray(result["frames"][max(i-1,0)]["state"]["displacement_m"])
        cell={"original_element_id":[ids],"material_id":[np.asarray(mesh["material_ids"])[active]],"region_id":[np.asarray(mesh["region_ids"])[active]]}
        for name in ("effective_stress_compression_positive_kpa","void_ratio","preconsolidation_pressure_kpa","plastic_strain_compression_positive"):
            cell[name]=[np.asarray([state["material_states"][j][name] for j in ids])]
        filename=f"stage-{i:03}.vtu"
        meshio.write_points_cells(output/filename,nodes[used],[("tetra",mapping[tets[active]])],
            point_data={"original_node_id":used,"displacement_m":u[used],"stage_displacement_increment_m":(u-previous)[used]},cell_data=cell)
        checked=meshio.read(output/filename)
        assert np.array_equal(checked.cells_dict["tetra"],mapping[tets[active]])
        assert np.allclose(checked.point_data["displacement_m"],u[used])
        ET.SubElement(data,"DataSet",timestep=str(i),group=frame["name"],part="0",file=filename)
        row={"stage_index":i,"stage":frame["name"],"active_elements":int(active.sum()),"relative_force_residual":frame["relative_force_residual"],
             "maximum_normalized_yield":frame["maximum_normalized_yield"]}
        for name,m in frame["monitors"].items():
            row[name+"_u_z_m"]=m["mean_displacement_z_m"];row[name+"_increment_z_m"]=m["stage_increment_z_m"]
            row[name+"_minimum_u_z_m"]=m.get("minimum_displacement_z_m")
            row[name+"_maximum_u_z_m"]=m.get("maximum_displacement_z_m")
        rows.append(row)
    ET.ElementTree(collection).write(output/"stages.pvd",encoding="utf-8",xml_declaration=True)
    with (output/"stages.csv").open("w",encoding="utf-8-sig",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    with (output/"attempts.csv").open("w",encoding="utf-8-sig",newline="") as f:
        fields=list(result["attempts"][0]) if result["attempts"] else ["stage","accepted"]
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(result["attempts"])
    package={"schema":"soil-fem-stage-export/1","result_sha256":hashlib.sha256(result_file.read_bytes()).hexdigest(),
             "mesh_sha256":hashlib.sha256(mesh_file.read_bytes()).hexdigest(),"vtu_roundtrip":"passed","frame_count":len(rows)}
    (output/"export-verification.json").write_text(json.dumps(package,indent=2),encoding="utf-8")
    print(f"Exported {len(rows)} active-domain frames, PVD and CSV: {output}")


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("mesh",type=Path);p.add_argument("result",type=Path);p.add_argument("output",type=Path)
    a=p.parse_args();export(a.mesh,a.result,a.output)
