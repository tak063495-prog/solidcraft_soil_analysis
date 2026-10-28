"""Version-bound, lazy-frame data for the forthcoming native SolidCraft viewer.

This is an external FEM postprocessor. It does not execute a CAD design, smooth
stress to nodes, or certify a run as mechanically valid. A completed, audited
native-capture job is required; the existing raw result remains authoritative.
"""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from high_order_io import state_fields
from import_solidcraft_capture import equal,validate
from verify_acceptance_packet import sha

FACES=((1,2,3),(0,3,2),(0,1,3),(0,2,1))
EDGES=((0,1),(1,2),(2,0),(0,3),(1,3),(2,3))
MAX_FILE=256*1024**2
MAX_TOTAL=512*1024**2


def read(path):
    if path.stat().st_size>MAX_FILE:raise ValueError('viewer input exceeds 256 MiB limit')
    def reject(value):raise ValueError('non-finite JSON: '+value)
    def real(value):
        number=float(value)
        if not math.isfinite(number):reject(value)
        return number
    return json.loads(path.read_text(encoding='utf-8'),parse_constant=reject,parse_float=real)


def clean(value):
    # Material-specific unavailable fields (MCC pc on MC, for example) remain
    # null. They must never appear as zero-valued physical measurements.
    if isinstance(value,np.ndarray):return clean(value.tolist())
    if isinstance(value,list):return [clean(v) for v in value]
    if isinstance(value,dict):return {k:clean(v) for k,v in value.items()}
    if isinstance(value,float) and not np.isfinite(value):return None
    return value


def boundary(cells,active,xyz):
    owners={}
    for ei in np.flatnonzero(active):
        cell=cells[ei]
        for local in FACES:
            face=[int(cell[a]) for a in local]
            if len(cell)==10:
                face += [int(cell[4+next(i for i,e in enumerate(EDGES) if set(e)=={a,b})])
                         for a,b in ((local[0],local[1]),(local[1],local[2]),(local[2],local[0]))]
            key=tuple(sorted(face[:3]))
            records=owners.setdefault(key,[])
            if records and set(records[0][1])!=set(face):raise ValueError('nonconforming shared quadratic face')
            records.append((int(ei),face))
            if len(records)>2:raise ValueError('nonmanifold active mesh face')
    result=[]
    for records in owners.values():
        if len(records)!=1:continue
        ei,face=records[0]
        corners=xyz[face[:3]];normal=np.cross(corners[1]-corners[0],corners[2]-corners[0])
        if np.dot(normal,corners.mean(axis=0)-xyz[cells[ei,:4]].mean(axis=0))<0:
            face=[face[0],face[2],face[1],face[5],face[4],face[3]] if len(face)==6 else [face[0],face[2],face[1]]
        result.append(dict(element_index=ei,nodes=face))
    return result


def scene_frame(mesh,run,index):
    return _scene_frame(mesh,run,index)


def _scene_frame(mesh,run,index,*,validated_mixed=False):
    mixed=None
    if run.get('schema') in ('soil-fem-mixed-mcc-stage-result/1','soil-fem-projected-mcc-stage-result/1','soil-fem-mixed-mcc-stage-result/2','soil-fem-projected-mcc-stage-result/2'):
        # These fields change the display contract. Validate the whole path,
        # including DEATH offsets and frozen histories, before extracting one
        # frame. Old native readers must reject the version-2 frame explicitly.
        import mixed_mcc_io
        if not validated_mixed:mixed_mcc_io.validate(mesh,run)
        mixed=mixed_mcc_io.frame_fields(mesh,run,run['frames'][index])
    elif any(f.get('state',{}).get('mixed_mcc') is not None for f in run.get('frames',[])):
        raise ValueError('mixed primary state requires its own result schema')
    xyz=np.asarray(mesh['nodes'],dtype=float)
    cells=np.asarray(mesh.get('tetrahedra10') or mesh['tetrahedra'],dtype=int)
    frame=run['frames'][index];state=frame['state'];layout=run.get('integration_points')
    if xyz.ndim!=2 or xyz.shape[1]!=3 or not np.isfinite(xyz).all():raise ValueError('invalid mesh positions')
    if cells.ndim!=2 or cells.shape[1] not in (4,10) or cells.min()<0 or cells.max()>=len(xyz):raise ValueError('invalid mesh connectivity')
    if layout is None:
        if run.get('schema')!='soil-fem-stage-result/1' or cells.shape[1]!=4:raise ValueError('missing integration layout')
        layout=[]
        for ei,cell in enumerate(cells):
            p=xyz[cell]
            layout.append(dict(element_index=ei,position_m=p.mean(axis=0).tolist(),
                weight_m3=float(np.linalg.det((p[1:]-p[0]).T)/6)))
    active=np.asarray(state['active_elements'])
    if active.dtype!=bool or active.shape!=(len(cells),):raise ValueError('invalid active element mask')
    if not active.any():raise ValueError('no active continuum cells to display')
    u=np.asarray(state['displacement_m'],dtype=float)
    initial=np.asarray(run.get('reference_displacement_m',np.zeros_like(xyz)),dtype=float)
    reference=np.asarray(state.get('display_reference_displacement_m') or initial,dtype=float)
    previous=np.asarray(run['frames'][max(0,index-1)]['state']['displacement_m'],dtype=float)
    if index and run['case']['stages'][index-1]['kind']=='birth':
        previous=np.asarray(next(a for a in run['transition_audits'] if a['stage']==frame['name'])['prepared_state']['displacement_m'],dtype=float)
    if any(a.shape!=xyz.shape or not np.isfinite(a).all() for a in (u,initial,reference,previous)):raise ValueError('invalid nodal displacement field')
    used=np.unique(cells[active]);ids=np.flatnonzero(active)
    weights=np.asarray([p['weight_m3'] for p in layout],dtype=float)
    element_ids=np.asarray([p['element_index'] for p in layout],dtype=int)
    positions=np.asarray([p['position_m'] for p in layout],dtype=float)
    if weights.ndim!=1 or not np.isfinite(weights).all() or np.any(weights<=0):raise ValueError('invalid integration weights')
    if positions.shape!=(len(weights),3) or not np.isfinite(positions).all():raise ValueError('invalid integration positions')
    if element_ids.shape!=weights.shape or element_ids.min()<0 or element_ids.max()>=len(cells):raise ValueError('invalid integration element index')
    grouped=[[] for _ in cells]
    for qi,ei in enumerate(element_ids):grouped[ei].append(qi)
    groups=[np.asarray(grouped[ei],dtype=int) for ei in ids]
    if any(not len(g) for g in groups):raise ValueError('active cell has no integration points')
    states=state['material_states']
    if len(states)!=len(layout):raise ValueError('integration history length mismatch')
    raw=state_fields(states)
    if mixed is not None:raw.update(mixed['raw'])
    for name in raw:
        if name=='material_model':continue
        if any(name in s and not np.isfinite(np.asarray(s[name],dtype=float)).all() for s in states):raise ValueError('non-finite stored material field')
    # Missing model fields from state_fields are represented with NaN internally
    # and null in JSON. Required state data must always be finite.
    for name in ('effective_stress_compression_positive_kpa','plastic_strain_compression_positive','void_ratio'):
        if not np.isfinite(raw[name]).all():raise ValueError('non-finite required material field')
    raw_models=raw.pop('material_model')
    if any(len(set(raw_models[g].tolist()))!=1 for g in groups):raise ValueError('mixed constitutive models within one cell')
    materials=np.asarray(state.get('element_material_ids') or mesh['material_ids'])
    regions=np.asarray(mesh['region_ids'])
    if materials.shape!=(len(cells),) or regions.shape!=(len(cells),):raise ValueError('cell identifier length mismatch')
    qp=np.flatnonzero(active[element_ids])
    cell_fields={name:np.asarray([np.average(values[g],axis=0,weights=weights[g]) for g in groups]) for name,values in raw.items()}
    cell_fields.update(material_id=materials[ids],region_id=regions[ids],material_model=np.asarray([raw_models[g[0]] for g in groups]))
    node_fields=dict(displacement_m=(u-reference)[used],stage_displacement_increment_m=(u-previous)[used],
                     total_displacement_from_seed_m=u[used],display_reference_displacement_m=reference[used])
    point_fields={name:values[qp] for name,values in raw.items()}
    point_fields['material_model']=raw_models[qp]
    point_fields['material_id']=materials[element_ids[qp]]
    # All connectivity uses original node ids. No interpolation or smoothing is
    # hidden here; a viewport can tessellate the supplied Tri6 geometry explicitly.
    surface=boundary(cells,active,xyz)
    view=dict(schema='soil-fem-view-frame/1',index=index,name=frame['name'],
        active_element_ids=ids,active_node_ids=used,boundary_faces=surface,
        node_fields=node_fields,cell_fields=cell_fields,
        monitors=frame['monitors'],relative_force_residual=frame['relative_force_residual'],
        maximum_normalized_yield=frame['maximum_normalized_yield'])
    points=dict(schema='soil-fem-view-integration-points/1',frame_index=index,
        integration_point_ids=qp,element_ids=element_ids[qp],positions_m=positions[qp],weights_m3=weights[qp],fields=point_fields)
    if mixed is not None:
        view['schema']='soil-fem-view-frame/2'
        points['schema']='soil-fem-view-integration-points/2'
        keys=mixed['keys']
        primary=dict(original_node_ids=[node for node,material in keys],
            material_ids=[material for node,material in keys],
            positions_m=xyz[[node for node,material in keys]],
            mean_effective_stress_compression_positive_kpa=mixed['primary'])
        if 'projected_volume_dofs' in mixed:
            primary.update(cumulative_material_volume_strain_compression_positive=mixed['projected_volume_dofs'],
                frozen_compatibility_offset_m3=mixed['projected_compatibility_offsets'])
        view['mixed_primary']=primary
        view['field_semantics']=dict(
            source_result_schema=run['schema'],
            material_stress='effective_stress_compression_positive_kpa; raw GP constitutive history',
            primary_pressure='material-split P1 mean effective pressure; not pore pressure or a displacement-node average',
            equilibrium_stress='work_conjugate_effective_stress_compression_positive_kpa' if 'projected_volume_dofs' in mixed else 'effective_stress_compression_positive_kpa',
            volume='material and kinematic cumulative volume strains remain distinct',
            primary_identity='(original_node_ids, material_ids); coincident material-boundary values remain separate',
            inactive='primary fields contain active corner/material pairs only; removed GP history stays in the authoritative result')
        view['mixed_equation_audit']=mixed['metrics']
    from solidcraft_structure_fields import fields as structure_fields
    structures=structure_fields(mesh,run,index)
    if structures is not None:
        view['schema']='soil-fem-view-frame/3';points['schema']='soil-fem-view-integration-points/3'
        view['structures']=structures
    return clean(view),clean(points)


def export(job,output):
    job=job.resolve();output=output.resolve();analysis=job/'analysis';prepared=analysis/'ground';conformal=job/'conformal'
    capture,ground,captured_stages,_=validate(conformal,manifest_name='native-capture.json')
    link=read(conformal/'native-capture-import.json')
    if link.get('passed') is not True or link.get('input_sha256')!=capture['input_sha256'] or link.get('source_capture_sha256')!=sha(conformal/'native-capture.json'):
        raise ValueError('native capture import linkage mismatch')
    result_path=analysis/'stage-results.json';mesh_path=prepared/'mesh.json';accepted_path=analysis/'accepted-stage-run.json'
    accepted=read(accepted_path);run=read(result_path);mesh=read(mesh_path)
    imported=read(prepared/'model-input-package.json');original=imported.get('original_inventory',{})
    if imported.get('schema')!='soil-fem-prepared-model-import/1' or original.get('schema')!='soil-fem-prepared-model/1':raise ValueError('prepared geometry import record required')
    inventory=original.get('files_sha256',{})
    if not {'mesh.json','bodies.json','mesh.provenance.json','native-capture.json'}<=set(inventory):raise ValueError('prepared geometry inventory is incomplete')
    geometry_sources=[]
    for name,digest in inventory.items():
        if not isinstance(name,str) or '/' in name or '\\' in name or name in ('.','..'):raise ValueError('invalid prepared inventory filename')
        path=prepared/('prepared-ground.json' if name=='ground.json' else name)
        if not path.resolve().is_relative_to(prepared) or sha(path)!=digest:raise ValueError('prepared geometry file hash mismatch')
        geometry_sources.append(path)
    if imported.get('ground_sha256')!=sha(prepared/'ground.json') or sha(prepared/'native-capture.json')!=sha(conformal/'native-capture.json'):
        raise ValueError('prepared geometry refers to a different captured input')
    provenance=read(prepared/'mesh.provenance.json');bodies=read(prepared/'bodies.json')
    if provenance.get('cad_manifest_sha256')!=sha(prepared/'bodies.json'):raise ValueError('mesh/CAD manifest fingerprint mismatch')
    body_records=bodies.get('regions',[])+bodies.get('excavation_tools',[])
    if {r['body']:(r['step'],r['source_sha256']) for r in body_records}!={r['analysis_body']:(r['step'],r['source_sha256']) for r in capture['bodies']}:
        raise ValueError('meshed bodies differ from native captured geometry')
    inputs={'mesh.json':mesh_path,'ground.json':prepared/'ground.json','stages.json':analysis/'stages.json','stage-results.json':result_path}
    if accepted.get('schema')!='soil-fem-accepted-staged-run/1':raise ValueError('accepted staged job required')
    recorded=accepted.get('inputs_and_results_sha256',{})
    names=[name.replace('\\','/').rsplit('/',1)[-1] for name in recorded]
    if len(names)!=4 or set(names)!=set(inputs):raise ValueError('invalid accepted job file inventory')
    for name,digest in recorded.items():
        if sha(inputs[name.replace('\\','/').rsplit('/',1)[-1]])!=digest:raise ValueError('accepted job input/result hash mismatch')
    if not equal(read(inputs['ground.json']),ground):raise ValueError('job ground differs from native capture')
    from prepare_native_plates import verify as verify_plate_expansion
    stages=read(inputs['stages.json'])
    verify_plate_expansion(conformal, prepared, mesh, stages, captured_stages)
    solver=analysis/'solver/soil-fem-3d.exe';solver_hash=sha(solver)
    if accepted.get('solver_binary_sha256')!=solver_hash or run.get('provenance',{}).get('solver_binary_sha256')!=solver_hash:raise ValueError('job solver fingerprint mismatch')
    if run.get('complete') is not True or run.get('failure_reason') is not None:raise ValueError('incomplete analysis cannot be displayed as accepted')
    projected=run.get('schema') in ('soil-fem-projected-mcc-stage-result/1','soil-fem-projected-mcc-stage-result/2')
    structural=isinstance(run.get('case',{}).get('beam_system'),dict)
    if not projected and run.get('schema') not in {f'soil-fem-stage-result/{i}' for i in range(1,15)}:raise ValueError('unsupported result schema for native viewer')
    if run['schema']=='soil-fem-projected-mcc-stage-result/2' and not structural:raise ValueError('structural projected source lacks members')
    if projected:
        from verify_projected_stage_package import audit
        audit(mesh,ground,stages,run)
    if not 1<=len(run.get('frames',[]))<=257:raise ValueError('viewer requires 1..257 frames')
    if not equal(run['last_committed'],run['frames'][-1]['state']):raise ValueError('final frame differs from committed state')
    verification=read(analysis/'stage-verification.json');exports=read(analysis/'results/export-verification.json')
    if verification.get('vtu_roundtrip')!='passed' or (run['schema']!='soil-fem-stage-result/1' and verification.get('raw_integration_point_roundtrip')!='passed'):
        raise ValueError('independent force and raw field audit required')
    if len(verification.get('frames',[]))!=len(run['frames']):raise ValueError('independent audit frame count mismatch')
    if exports.get('result_sha256')!=sha(result_path) or exports.get('mesh_sha256')!=sha(mesh_path):raise ValueError('export audit input/result mismatch')
    if projected and (verification.get('schema')!='soil-fem-stage-verification/3' or verification.get('independent_load_and_domain_audit')!='passed'):
        raise ValueError('projected load/domain and work-conjugate audit required')
    for frame in verification['frames']:
        r=frame['independent_force_residual']
        if not np.isfinite(r) or r>stages['settings']['relative_force_tolerance']*1.01:raise ValueError('independent residual exceeds job tolerance')
    source_paths=[*inputs.values(),*geometry_sources,prepared/'model-input-package.json',accepted_path,solver,analysis/'stage-verification.json',analysis/'results/export-verification.json',conformal/'native-capture.json',conformal/'native-capture-import.json']
    hashes={str(p.relative_to(job)).replace('\\','/'):sha(p) for p in source_paths}
    output.mkdir(parents=True,exist_ok=False);files={};total=0
    def write(name,value,*,track=True):
        nonlocal total
        data=json.dumps(value,ensure_ascii=False,allow_nan=False,separators=(',',':')).encode('utf-8')
        total+=len(data)
        if len(data)>MAX_FILE or total>MAX_TOTAL:raise ValueError('view packet exceeds bounded size; incomplete output retained')
        with (output/name).open('xb') as f:f.write(data)
        if track:files[name]=sha(output/name)
    geometry=dict(schema='soil-fem-view-geometry/1',length_unit='m',axis='Z-up',nodes_m=mesh['nodes'],
        cells=mesh.get('tetrahedra10') or mesh['tetrahedra'],cell_kind='tetra10' if mesh.get('tetrahedra10') else 'tetra4',
        regions=mesh['regions'],node_sets=mesh.get('node_sets',{}),crs=mesh.get('crs'),origin_m=mesh.get('origin_m',[0.,0.,0.]))
    write('geometry.json',geometry);frames=[]
    for i in range(len(run['frames'])):
        # The projected audit above has checked every frame, frozen history and
        # domain transition. Do not repeat the whole path for each extracted
        # frame (quadratic work in the number of construction stages).
        view,points=_scene_frame(mesh,run,i,validated_mixed=projected);name=f'frame-{i:03}.json';point_name=f'points-{i:03}.json'
        write(name,view);write(point_name,points)
        frames.append(dict(index=i,name=view['name'],frame=name,integration_points=point_name))
    if any(sha(job/p)!=digest for p,digest in hashes.items()):raise ValueError('job changed while exporting view packet')
    # Check the STEP/input archive again before publishing a completion marker.
    checked,_,_,_=validate(conformal,manifest_name='native-capture.json')
    if not equal(checked,capture):raise ValueError('native capture changed while exporting view packet')
    manifest=dict(schema='soil-fem-solidcraft-view/1',complete=True,input_sha256=capture['input_sha256'],
        cad_sha256=read(conformal/'snapshot-canonical.json')['cad_sha256'],
        exporter_binary_sha256=capture['exporter_binary_sha256'],solver_binary_sha256=solver_hash,
        geometry='geometry.json',frames=frames,files_sha256=files,job_sources_sha256=hashes,
        semantics=dict(length_unit='m',stress_unit='kPa',z_positive='up; settlement is negative',
            stress_sign='compression positive',stress_components=['xx','yy','zz','xy','yz','xz'],
            displacement_m='total displacement minus the explicit display reference',
            stage_displacement_increment_m='increment from previous frame, or prepared BIRTH state',
            total_displacement_from_seed_m='unreset cumulative displacement from the original seed',
            cell_fields='integration-weighted volume averages; constant per cell; not smoothed to nodes',
            integration_points='raw stored material values; never averaged',
            unavailable_value='null; constitutive-model-specific unavailable quantity',
            boundary_faces='original node ids; Tri3 or Tri6; tessellation is a viewer approximation'),
        omissions=['beam/anchor/contact result layers','water fields','principal-direction glyphs'],
        scope='Native-viewer data foundation only; no native viewport implementation or P5 completion is implied.')
    if projected:
        manifest['schema']='soil-fem-solidcraft-view/2'
        manifest['semantics'].update(source_result_schema=run['schema'],
            material_stress='effective_stress_compression_positive_kpa; raw constitutive GP history',
            equilibrium_stress='work_conjugate_effective_stress_compression_positive_kpa',
            primary_pressure='material-split P1 mean effective stress; distinct from pore pressure',
            mixed_primary='active corner/material pairs; no cross-material averaging',
            volume='cumulative material volume and kinematic divergence remain distinct')
    if structural:
        manifest['schema']='soil-fem-solidcraft-view/3'
        manifest['semantics'].update(source_result_schema=run['schema'],
            structures='raw beam section forces and DKT plate samples in their reported local axes; archived inactive state retained',
            structural_displacement='from CAD; independent of soil display RESET',
            structural_units='beam kN and kN m; plate membrane kN/m and bending kN m/m; rotations rad')
        manifest['omissions']=['anchor/contact result layers','water fields']
        manifest['scope']='Versioned continuum and linear frame/plate display data; native rendering and mechanical/spatial acceptance require separate evidence.'
    write('view-manifest.json',manifest,track=False)
    print('Native viewer data:',len(frames),'frames; explicit nodal, cell-average and raw integration-point fields')
    return manifest


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('job',type=Path);p.add_argument('output',type=Path);a=p.parse_args();export(a.job,a.output)
