"""Conformal multi-body volume mesh using Gmsh's fragment ancestry map."""
import argparse
import hashlib
import json
from pathlib import Path
import gmsh
import numpy as np


EDGES=((0,1),(1,2),(2,0),(0,3),(1,3),(2,3))
FACES=((1,2,3),(0,3,2),(0,1,3),(0,2,1))

def generate(manifest_path, output, size, order=1, quadrature_order=3):
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != "soil-fem-cad-bodies/1" or manifest.get("units") != "m-kN-kPa":
        raise ValueError("invalid CAD manifest")
    if output.exists() or not np.isfinite(size) or size <= 0:
        raise ValueError("output exists or invalid mesh size")
    if order not in (1,2) or quadrature_order not in (3,4,5):
        raise ValueError('element order must be 1/2; volume quadrature order must be 3/4/5')
    regions = manifest["regions"]
    tools = manifest.get("excavation_tools", [])
    if not regions or len({r["body"] for r in regions}) != len(regions):
        raise ValueError("regions missing or duplicated")
    all_names=[r['body'] for r in regions+tools]
    if len(set(all_names))!=len(all_names) or any(not n for n in all_names):
        raise ValueError('material and excavation body names must be distinct and nonempty')
    if len({r['id'] for r in regions})!=len(regions):
        raise ValueError('duplicate material region IDs')
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal",0)
        gmsh.logger.start()
        gmsh.option.setNumber("General.NumThreads",1)
        gmsh.option.setString("Geometry.OCCTargetUnit","MM")
        imported = []
        for r in regions + tools:
            file = manifest_path.parent / r["step"]
            if hashlib.sha256(file.read_bytes()).hexdigest() != r["source_sha256"]:
                raise ValueError("body STEP hash mismatch")
            entities = gmsh.model.occ.importShapes(str(file.resolve()))
            volumes = [x for x in entities if x[0] == 3]
            if len(volumes) != 1:
                raise ValueError(f"body {r['body']} must contain one closed solid")
            volume = gmsh.model.occ.getMass(*volumes[0])*1e-9
            if not np.isclose(volume,r["volume_m3"],rtol=1e-8):
                raise ValueError("CAD/STEP volume mismatch")
            imported.extend(volumes)
        if len(imported) > 1:
            _, ancestry = gmsh.model.occ.fragment(imported[:1],imported[1:])
        else:
            ancestry = [imported]
        owner = {}
        for r, children in zip(regions,ancestry[:len(regions)]):
            for dim, tag in children:
                if dim != 3:
                    continue
                if tag in owner:
                    raise ValueError("CAD solids overlap in positive volume; repair geometry")
                owner[tag] = r
        selections = {r["body"]: {tag for dim,tag in children if dim == 3} & set(owner)
                      for r,children in zip(tools,ancestry[len(regions):])}
        if any(not tags for tags in selections.values()):
            raise ValueError("excavation tool has no positive-volume intersection with soil")
        # Tool-only volume above/outside the ground is a selector, not additional soil.
        extra = [(dim,tag) for dim,tag in gmsh.model.occ.getEntities(3) if tag not in owner]
        if extra:
            gmsh.model.occ.remove(extra, recursive=True)
        gmsh.model.occ.synchronize()
        if set(owner) != {tag for dim,tag in gmsh.model.getEntities(3)}:
            raise ValueError("lost body ancestry during fragmentation")
        for r in regions:
            tags = [tag for tag, value in owner.items() if value["id"] == r["id"]]
            if not tags:
                raise ValueError("body disappeared during fragmentation")
            volume = sum(gmsh.model.occ.getMass(3,tag) for tag in tags)*1e-9
            if not np.isclose(volume,r["volume_m3"],rtol=1e-8):
                raise ValueError("fragment volume changed")
            group = gmsh.model.addPhysicalGroup(3,tags)
            gmsh.model.setPhysicalName(3,group,r["body"])
        gmsh.option.setNumber("Mesh.MeshSizeMin",size*1000)
        gmsh.option.setNumber("Mesh.MeshSizeMax",size*1000)
        gmsh.option.setNumber("Mesh.ElementOrder",order)
        gmsh.option.setNumber("Mesh.SecondOrderLinear",0)
        if order==2:gmsh.option.setNumber("Mesh.HighOrderOptimize",1)
        gmsh.model.mesh.generate(3)
        ordering=list(range(4))
        if order==2:
            local=np.array(gmsh.model.mesh.getElementProperties(11)[4]).reshape(-1,3)
            ordering=list(range(4))
            for a,b in EDGES:
                distances=np.linalg.norm(local-(local[a]+local[b])/2,axis=1)
                middle=int(np.argmin(distances))
                if distances[middle]>1e-12:raise ValueError('unknown Gmsh Tet10 edge ordering')
                ordering.append(middle)
        tags, coords, _ = gmsh.model.mesh.getNodes()
        by_tag = {int(tag):p.tolist() for tag,p in zip(tags,np.asarray(coords).reshape(-1,3)*0.001)}
        tets, materials, region_ids = [], [], []
        excavation_sets = {name: [] for name in selections}
        for tag,r in sorted(owner.items()):
            kinds, _, blocks = gmsh.model.mesh.getElements(3,tag)
            for kind,block in zip(kinds,blocks):
                if kind != (4 if order==1 else 11):
                    raise ValueError("mesher returned an unexpected tetrahedral order")
                batch = np.asarray(block).reshape(-1,len(ordering))[:,ordering].tolist()
                for name, tags in selections.items():
                    if tag in tags:
                        excavation_sets[name].extend(range(len(tets), len(tets)+len(batch)))
                tets.extend(batch)
                materials.extend([r["material_id"]]*len(batch))
                region_ids.extend([r["id"]]*len(batch))
        used = sorted({int(tag) for tet in tets for tag in tet})
        indices = {tag:i for i,tag in enumerate(used)}
        nodes = np.asarray([by_tag[tag] for tag in used])
        if len(np.unique(np.round(nodes,10),axis=0)) != len(nodes):
            raise ValueError("coincident disconnected mesh nodes")
        full_elements = [[indices[int(tag)] for tag in tet] for tet in tets]
        for tet in full_elements:
            p = nodes[tet[:4]]
            if np.linalg.det((p[1:]-p[0]).T) < 0:
                old=tet.copy();permutation=[0,2,1,3];tet[:4]=[old[i] for i in permutation]
                if order==2:
                    for k,(a,b) in enumerate(EDGES):
                        pair={permutation[a],permutation[b]}
                        old_edge=next(j for j,e in enumerate(EDGES) if set(e)==pair)
                        tet[k+4]=old[old_edge+4]
        elements=[t[:4] for t in full_elements]
        def face_nodes(face,owner):
            cell=full_elements[owner]
            for f in FACES:
                if set(cell[a] for a in f)!=set(face):continue
                ids=[cell[a] for a in f]
                if order==2:
                    for a,b in ((0,1),(1,2),(2,0)):
                        edge=next(i for i,e in enumerate(EDGES) if set(e)=={f[a],f[b]})
                        ids.append(cell[edge+4])
                return ids
            raise ValueError('lost owner face')
        bounds = [nodes.min(axis=0),nodes.max(axis=0)]
        sets = {}
        for axis,names in enumerate((("xmin","xmax"),("ymin","ymax"),("bottom","top"))):
            for side,name in enumerate(names):
                sets[name] = np.flatnonzero(np.isclose(nodes[:,axis],bounds[side][axis],atol=1e-8,rtol=0)).tolist()
        mesh = {"schema":"soil-fem-mesh/2" if order==1 else "soil-fem-mesh/3","units":"m-kN-kPa","axis":"z_up","crs":"LOCAL",
                "origin_m":[0.,0.,0.],"nodes":nodes.tolist(),"tetrahedra":elements,
                "material_ids":materials,"region_ids":region_ids,"node_sets":sets,
                "regions":[{k:r[k] for k in ('id','body','material_id','bounds_m','volume_m3','source_sha256')} for r in regions]}
        if order==2:mesh.update(tetrahedra10=full_elements,quadrature_order=quadrature_order)
        faces = {}
        for ei,tet in enumerate(elements):
            for f in FACES:
                key=tuple(sorted(tet[a] for a in f));faces.setdefault(key,[]).append(ei)
        face_sets={}
        for axis,names in enumerate((('xmin','xmax'),('ymin','ymax'),('bottom','top'))):
            for side,name in enumerate(names):
                face_sets[name]=[list(f) for f,owners in faces.items() if len(owners)==1
                    and np.allclose(nodes[face_nodes(f,owners[0]),axis],bounds[side][axis],atol=1e-8,rtol=0)]
        mesh['face_sets']=face_sets
        region_names={r['id']:r['body'] for r in regions}
        for face,owners in faces.items():
            if len(owners)==2 and region_ids[owners[0]]!=region_ids[owners[1]]:
                names=sorted(region_names[region_ids[e]] for e in owners)
                face_sets.setdefault('interface:'+':'.join(names),[]).append(list(face))
        exterior=[list(f) for f,owners in faces.items() if len(owners)==1]
        free=[]
        for f in exterior:
            owner=faces[tuple(f)][0];p=nodes[f];normal=np.cross(p[1]-p[0],p[2]-p[0]);center=p.mean(axis=0)
            if normal@(center-nodes[elements[owner]].mean(axis=0))<0:normal=-normal
            if normal[2]>1e-8*np.linalg.norm(normal):free.append(f)
        face_sets['exterior']=exterior;face_sets['free_surface']=free
        free_keys={tuple(f) for f in free}
        sets['all']=list(range(len(nodes)))
        sets['free_surface']=sorted({n for f in free for n in face_nodes(f,faces[tuple(f)][0])})
        if excavation_sets:
            mesh["excavation_sets"] = excavation_sets
            # Surface monitors are tied to actual cut facets, not tool bounding boxes.
            for name, selection in excavation_sets.items():
                selected=set(selection);floor=set();wall=set();surface=set()
                floor_faces=[];wall_faces=[];surface_faces=[]
                for face, owners in faces.items():
                    p=nodes[list(face)];center=p.mean(axis=0)
                    if len(owners)==2 and sum(i in selected for i in owners)==1:
                        retained=next(i for i in owners if i not in selected)
                        n=np.cross(p[1]-p[0],p[2]-p[0])
                        if n@(center-nodes[elements[retained]].mean(axis=0))<0:n=-n
                        is_floor=n[2]>0.5*np.linalg.norm(n)
                        (floor if is_floor else wall).update(face_nodes(face,retained))
                        (floor_faces if is_floor else wall_faces).append(list(face))
                    if len(owners)==1 and owners[0] not in selected and face in free_keys:
                        surface.update(face_nodes(face,owners[0]));surface_faces.append(list(face))
                sets[f"{name}:floor"]=sorted(floor);sets[f"{name}:wall"]=sorted(wall)
                sets[f"{name}:surrounding_surface"]=sorted(surface)
                face_sets[f'{name}:floor']=floor_faces;face_sets[f'{name}:wall']=wall_faces
                face_sets[f'{name}:surrounding_surface']=surface_faces
        output.write_text(json.dumps(mesh,allow_nan=False),encoding="utf-8")
        gmsh.write(str(output.with_suffix(".msh")))
        output.with_suffix(".provenance.json").write_text(json.dumps({
            "gmsh":gmsh.__version__,"cad_manifest_sha256":hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
            "ground_sha256":manifest["ground_sha256"],"mesh_size_m":size,"element_order":order,"quadrature_order":quadrature_order if order==2 else None,
            "excavation_boolean_volume_m3":{name:sum(gmsh.model.occ.getMass(3,tag) for tag in tags)*1e-9 for name,tags in selections.items()}},indent=2),encoding="utf-8")
        print(f"Conformal mesh: {len(nodes)} nodes, {len(tets)} tetrahedra, {len(regions)} bodies")
    finally:
        output.with_suffix(".gmsh.log").write_text("\n".join(gmsh.logger.get()),encoding="utf-8")
        gmsh.logger.stop()
        gmsh.finalize()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest",type=Path)
    parser.add_argument("output",type=Path)
    parser.add_argument("--mesh-size-m",type=float,default=2.)
    parser.add_argument("--order",type=int,default=1,choices=(1,2))
    parser.add_argument("--quadrature-order",type=int,default=3,choices=(3,4,5))
    a = parser.parse_args()
    generate(a.manifest,a.output,a.mesh_size_m,a.order,a.quadrature_order)
