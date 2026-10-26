"""Bind named, straight exterior or declared shared-interface faces to plates.

Tri6 is split into four Tri3 plates using all six soil nodes. This preserves
the source geometry and exact nodal transfer, but does not make the quadratic
soil trace equal to the piecewise-linear plate trace between the shared nodes.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import re
import numpy as np
import fem_reference as ref


def attach(mesh, case, specification):
    version=specification.get('schema')
    if version not in ('soil-fem-plate-surfaces/1','soil-fem-plate-surfaces/2'):
        raise ValueError('plate surface specification schema required')
    surfaces=specification.get('surfaces',[])
    if not surfaces or len(surfaces)>20:raise ValueError('require 1 to20 plate surfaces')
    xyz=np.asarray(mesh['nodes'],float)
    if xyz.ndim!=2 or xyz.shape[1]!=3 or not np.isfinite(xyz).all():raise ValueError('invalid mesh positions')
    result_mesh=copy.deepcopy(mesh);result_case=copy.deepcopy(case)
    if result_case.get('beam_system') is None:
        result_case['beam_system']=dict(schema='soil-fem-structures/1',nodes=[],sections=[],members=[],
            plate_sections=[],plates=[],rotation_length_scale_m=1.,maximum_condensed_dofs=3000)
    cfg=result_case['beam_system']
    if cfg['schema'] not in ('soil-fem-structures/1','soil-fem-beams/1'):raise ValueError('unsupported structure schema')
    cfg['schema']='soil-fem-structures/1';cfg.setdefault('plate_sections',[]);cfg.setdefault('plates',[])
    section_ids={s['id'] for s in cfg['plate_sections']}
    for section in specification.get('sections',[]):
        if section['id'] in section_ids:raise ValueError('plate section ID collision')
        cfg['plate_sections'].append(copy.deepcopy(section));section_ids.add(section['id'])
    node_ids={n['id'] for n in cfg['nodes']};member_ids={m['id'] for m in cfg['members']+cfg['plates']}
    if len(node_ids)!=len(cfg['nodes']) or len(member_ids)!=len(cfg['members'])+len(cfg['plates']):raise ValueError('duplicate structure identifiers')
    attached={}
    for node in cfg['nodes']:
        soil=None
        if node.get('soil_node_set') is not None:
            group=mesh['node_sets'].get(node['soil_node_set'],[])
            if len(group)!=1:raise ValueError('existing structural attachment must identify one soil node')
            soil=group[0]
        elif node.get('match_soil_node'):
            found=np.flatnonzero(np.linalg.norm(xyz-np.asarray(node['position_m']),axis=1)<=1e-8)
            if len(found)!=1:raise ValueError('ambiguous existing coordinate attachment')
            soil=int(found[0])
        if soil is not None:
            if soil in attached:raise ValueError('multiple existing rotations at one soil node; select their connection explicitly')
            if np.linalg.norm(xyz[soil]-node['position_m'])>1e-8:raise ValueError('existing attachment coordinate mismatch')
            attached[soil]=node['id']
    topology=ref.topology(mesh);regions={r['id']:r['body'] for r in mesh['regions']}
    selected=set();groups={};records=[];generated=[]

    def node_id(soil):
        soil=int(soil)
        if soil in attached:return attached[soil]
        name=f'plate_mesh_node_{soil}'
        if name in node_ids:raise ValueError('generated structural node ID collision')
        group=f'plate:mesh-node:{soil}'
        if group in result_mesh.setdefault('node_sets',{}):raise ValueError('generated attachment node set collision')
        result_mesh['node_sets'][group]=[soil]
        cfg['nodes'].append(dict(id=name,position_m=xyz[soil].tolist(),soil_node_set=group))
        node_ids.add(name);attached[soil]=name;generated.append(soil)
        return name

    for surface in surfaces:
        name=surface['id']
        if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_-]{0,79}',name) or name in groups or name in member_ids:
            raise ValueError('unique safe surface ID required')
        if surface['section_id'] not in section_ids:raise ValueError('unknown plate surface section')
        if not isinstance(surface['initially_active'],bool):raise ValueError('explicit initial plate activity required')
        attachment=surface.get('attachment','exterior')
        if attachment not in ('exterior','shared_interface'):
            raise ValueError('unknown plate surface attachment')
        if attachment=='shared_interface' and version!='soil-fem-plate-surfaces/2':
            raise ValueError('shared-interface plate attachment requires surface schema /2')
        opposite=surface.get('opposite_body')
        if attachment=='shared_interface' and (not isinstance(opposite,str) or not opposite or opposite==surface['owner_body']):
            raise ValueError('shared interface requires a distinct opposite body')
        if attachment=='exterior' and opposite is not None:
            raise ValueError('an exterior plate surface cannot declare an opposite body')
        faces=mesh.get('face_sets',{}).get(surface['face_set'],[])
        if not faces:raise ValueError('unknown or empty plate face set')
        cells=[];total_area=0.;source_area=0.;groups[name]=[];owner_pairs=[]
        for face_index,key in enumerate(sorted(tuple(sorted(f)) for f in faces)):
            if key in selected:raise ValueError('overlapping/duplicate plate surface faces')
            selected.add(key);owners=topology.get(key,[])
            owned=[e for e in owners if regions[mesh['region_ids'][e]]==surface['owner_body']]
            other=[e for e in owners if regions[mesh['region_ids'][e]]==opposite]
            if attachment=='exterior':
                if len(owners)!=1 or len(owned)!=1:
                    raise ValueError('plate surface must be exterior to its declared owner body')
            elif len(owners)!=2 or len(owned)!=1 or len(other)!=1:
                raise ValueError('plate interface must have exactly its two declared body owners')
            q=ref.face(mesh,owned[0],key);nodes=q['nodes'];points=xyz[nodes]
            if attachment=='shared_interface':
                opposing=ref.face(mesh,other[0],key)
                if set(map(int,nodes))!=set(map(int,opposing['nodes'])):
                    raise ValueError('shared plate interface has nonconforming trace nodes')
                normal=q['area_vector'].sum(axis=0);opposing_normal=opposing['area_vector'].sum(axis=0)
                if np.linalg.norm(normal+opposing_normal)>1e-10*max(np.linalg.norm(normal),1.):
                    raise ValueError('shared plate interface owner orientations disagree')
                owner_pairs.append(dict(owner_element=owned[0],opposite_element=other[0],face=list(key)))
            scale=max(np.linalg.norm(points.max(axis=0)-points.min(axis=0)),1.)
            if len(nodes)==6:
                expected=np.array([(points[i]+points[j])/2 for i,j in ((0,1),(1,2),(2,0))])
                if np.max(np.abs(points[3:]-expected))>1e-10*scale:
                    raise ValueError('curved Tri6 trace is not silently flattened; straight faces required')
                parts=((0,3,5),(3,1,4),(5,4,2),(3,4,5))
            elif len(nodes)==3:parts=((0,1,2),)
            else:raise ValueError('Tri3 or Tri6 trace required')
            source_area+=float(np.linalg.norm(q['area_vector'],axis=1).sum())
            outward=q['area_vector'].sum(axis=0)
            for part_index,part in enumerate(parts):
                cell=[int(nodes[i]) for i in part]
                cross=np.cross(xyz[cell[1]]-xyz[cell[0]],xyz[cell[2]]-xyz[cell[0]])
                if np.linalg.norm(cross)<1e-12*scale*scale:raise ValueError('degenerate plate trace triangle')
                if cross@outward<=0:raise ValueError('plate trace orientation is inconsistent with owner exterior')
                total_area+=float(np.linalg.norm(cross)/2);cells.append(cell)
                member=f'{name}_face{face_index}_part{part_index}'
                if member in member_ids:raise ValueError('generated plate member ID collision')
                member_ids.add(member);groups[name].append(member)
                plate=dict(id=member,nodes=[node_id(n) for n in cell],section_id=surface['section_id'],initially_active=surface['initially_active'])
                for field in ('constant_traction_kpa','boundary_traction_kpa'):
                    value=np.asarray(surface.get(field,[0.]*3),float)
                    if value.shape!=(3,) or not np.isfinite(value).all():raise ValueError('finite global plate traction required')
                    plate[field]=value.tolist()
                cfg['plates'].append(plate)
        if abs(total_area-source_area)>1e-10*max(source_area,1.):raise ValueError('plate and source trace areas differ')
        used=sorted({n for cell in cells for n in cell})
        records.append(dict(id=name,face_set=surface['face_set'],owner_body=surface['owner_body'],source_faces=len(faces),
            plate_triangles=len(cells),shared_soil_nodes=len(used),area_m2=total_area,source_area_m2=source_area,
            outward_orientation_checked=True,all_quadratic_trace_nodes_retained=True))
        if version=='soil-fem-plate-surfaces/2':
            records[-1].update(attachment=attachment,opposite_body=opposite,owner_pairs=owner_pairs)
    if len(cfg['nodes'])>500 or len(cfg['members'])+len(cfg['plates'])>2000:
        raise ValueError('surface plate system exceeds declared Rust structural budget; no truncation')
    for stage in result_case['stages']:
        if stage['kind'] in ('install_structures','remove_structures','install_beams','remove_beams'):
            stage['members']=[member for name in stage['members'] for member in groups.get(name,[name])]
            if len(set(stage['members']))!=len(stage['members']):raise ValueError('overlapping structural stage groups')
    binding=dict(schema='soil-fem-plate-surface-binding/1',surfaces=records,member_groups=groups,
        soil_node_to_structure_id={str(n):attached[n] for n in sorted(attached)},generated_soil_attachment_nodes=generated,
        geometry_and_connectivity_changed=False,soil_plate_coupling='Exact shared nodal translations; structural rotations shared by node ID',
        interpolation_limitation='Tri6 soil trace and four Tri3 plate traces agree at all six nodes, not pointwise between nodes',
        geometry_scope='Straight exterior Tri3/Tri6 faces; no persistent CAD face ID claim',goal_complete=False,p6_complete=False)
    if version=='soil-fem-plate-surfaces/2':
        binding.update(schema='soil-fem-plate-surface-binding/2',
            geometry_scope='Straight exterior or explicitly declared conforming shared-interface faces; no persistent CAD face ID claim',
            interface_limitation='Shared soil translations on both sides before excavation; no opening, slip or JOINT introduced. Owner-side orientation is explicit; retaining-wall mechanics and pointwise trace identity require separate verification.')
    return result_mesh,result_case,binding


def main():
    parser=argparse.ArgumentParser()
    for name in ('mesh','case','specification','output'):parser.add_argument(name,type=Path)
    args=parser.parse_args();inputs=[args.mesh,args.case,args.specification]
    data=[json.loads(p.read_text(encoding='utf-8-sig')) for p in inputs]
    mesh,case,binding=attach(*data)
    args.output.mkdir(parents=True,exist_ok=False)
    for name,value in [('mesh.json',mesh),('stages.json',case)]:
        (args.output/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    binding['source_sha256']={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs+[Path(__file__),Path(ref.__file__)]}
    binding['output_sha256']={name:hashlib.sha256((args.output/name).read_bytes()).hexdigest() for name in ('mesh.json','stages.json')}
    (args.output/'plate-surface-binding.json').write_text(json.dumps(binding,indent=2),encoding='utf-8')
    print(json.dumps(binding['surfaces']))


if __name__=='__main__':main()
