"""Validate and export unsmoothed primary mixed-MCC fields.

Mean effective pressure is distinct from pore pressure; material volume is
distinct from the kinematic volume strain. The latter agree weakly, not at
every integration point. Validation precedes creation of the output directory.
"""
import csv
import xml.etree.ElementTree as ET
import meshio
import numpy as np
import fem_reference as ref

SCHEMA='soil-fem-mixed-mcc-stage-result/1'
PROJECTED_SCHEMA='soil-fem-projected-mcc-stage-result/1'
PROJECTED_FORMULATION='tet10_u_material_p1_p_p1_volume_frozen_compatibility'
LEGACY_TANGENT='finite_difference_fixed_partition'
TANGENT_METHODS={LEGACY_TANGENT, 'implicit_branch_fixed_partition'}


def array(value,shape,name):
    result=np.asarray(value,dtype=float)
    if result.shape!=shape or not np.isfinite(result).all():raise ValueError('invalid mixed '+name)
    return result


def frame_fields(mesh,result,frame,check_geometry=False):
    layout=result['integration_points'];ng=len(layout);ne=len(mesh['tetrahedra'])
    state=frame['state'];mixed=state.get('mixed_mcc')
    if any(state.get(k) is not None for k in ('pore_pressure_kpa','hydrostatic_water_level_m','joint_water','joints','anchors','beams')):
        raise ValueError('mixed wet/member/interface export needs explicit supported coupling')
    if not isinstance(mixed,dict):raise ValueError('missing primary mixed MCC state')
    projected=result.get('schema')==PROJECTED_SCHEMA
    projection=mixed.get('projected_volume')
    if projected != isinstance(projection,dict):raise ValueError('mixed volume formulation/state mismatch')
    if projected and projection.get('formulation')!=PROJECTED_FORMULATION:raise ValueError('unknown projected MCC formulation')
    if mixed.get('tangent_method',LEGACY_TANGENT) not in TANGENT_METHODS:
        raise ValueError('unknown mixed MCC tangent method')
    active=np.asarray(state['active_elements'])
    if active.shape!=(ne,) or active.dtype!=bool or not active.any():raise ValueError('invalid mixed active elements')
    material=state.get('element_material_ids') or mesh['material_ids']
    if material!=mesh['material_ids']:raise ValueError('mixed material changes need an explicit supported policy')
    u=array(state['displacement_m'],(len(mesh['nodes']),3),'displacement')
    theta=array(mixed['cumulative_volume_strain_compression_positive'],(ng,),'cumulative volume')
    reference=array(mixed['reference_void_ratio'],(ng,),'reference void ratio')
    partitions=np.asarray(mixed['last_increment_partitions'])
    if partitions.shape!=(ng,) or partitions.dtype.kind not in 'iu' or np.any((partitions<1)|(partitions>100000)):
        raise ValueError('invalid mixed material partitions')
    states=state['material_states']
    if len(states)!=ng or any('preconsolidation_pressure_kpa' not in s for s in states):raise ValueError('MCC histories required')
    stress=array([s['effective_stress_compression_positive_kpa'] for s in states],(ng,6),'stress')
    void=array([s['void_ratio'] for s in states],(ng,),'void ratio')
    pc=array([s['preconsolidation_pressure_kpa'] for s in states],(ng,),'pc')
    array([s['plastic_strain_compression_positive'] for s in states],(ng,6),'plastic strain')
    if np.any(reference<=0) or np.any(void<=0) or np.any(pc<=0):raise ValueError('nonpositive mixed MCC state')
    with np.errstate(over='ignore',invalid='ignore'):
        expected_void=(1+reference)*np.exp(-theta)-1
    if not np.all(np.isfinite(expected_void)) or np.any(abs(expected_void-void)>1e-10*(1+abs(expected_void))):
        raise ValueError('mixed volume and void ratio are inconsistent')
    keys=[(d['node_index'],d['material_id']) for d in mixed['pressure_dofs']]
    expected={(node,material[e]) for e,tet in enumerate(mesh['tetrahedra']) if active[e] for node in tet}
    if len(keys)!=len(set(keys)) or set(keys)!=expected:raise ValueError('mixed pressure inventory differs from active corner/material pairs')
    primary=array([d['compression_positive_kpa'] for d in mixed['pressure_dofs']],(len(keys),),'primary pressure')
    lookup=dict(zip(keys,range(len(keys))))
    nodal_theta=offsets=None
    if projected:
        def nodal(name,field):
            rows=projection[name];inventory=[(d['node_index'],d['material_id']) for d in rows]
            if len(inventory)!=len(keys) or len(set(inventory))!=len(keys) or set(inventory)!=set(keys):
                raise ValueError('projected '+name+' inventory differs')
            values={k:d[field] for k,d in zip(inventory,rows)}
            return array([values[k] for k in keys],(len(keys),),name)
        nodal_theta=nodal('volume_dofs','compression_positive')
        offsets=nodal('compatibility_offsets','integral_m3')
    gp_pressure=np.full(ng,np.nan);kinematic=np.empty((ng,6));volume=np.zeros(len(keys));scale=np.zeros(len(keys));measure=0.
    weak_pressure=np.zeros(len(keys));pressure_scale=np.zeros(len(keys));mean=stress[:,:3].mean(axis=1)
    owner=np.empty(ng,dtype=int);cursor=0
    for e,tet in enumerate(mesh['tetrahedra']):
        q=ref.volume(mesh,e);nq=len(q['weight']);indices=np.arange(cursor,cursor+nq)
        if cursor+nq>ng:raise ValueError('mixed quadrature inventory too short')
        owner[indices]=e
        for local,i in enumerate(indices):
            p=layout[i]
            if p['element_index']!=e or p['local_index']!=local:raise ValueError('mixed GP owner/order mismatch')
            if check_geometry and (not np.allclose(p['position_m'],q['position'][local],rtol=0.,atol=1e-8)
                or not np.allclose(p['natural_coordinates'],q['natural'][local],rtol=0.,atol=1e-12)
                or not np.isclose(p['weight_m3'],q['weight'][local],rtol=1e-10,atol=1e-12)
                or not np.allclose(p['shape_functions'],q['shape'][local],rtol=0.,atol=1e-11)):
                raise ValueError('mixed integration geometry differs from the mesh')
        kinematic[indices]=-np.einsum('qai,i->qa',q['b'],u[q['nodes']].reshape(-1))
        if active[e]:
            pi=np.array([lookup[(node,material[e])] for node in tet]);l=np.column_stack((1-q['natural'].sum(axis=1),q['natural']))
            gp_pressure[indices]=l@primary[pi]
            div=-kinematic[indices,:3].sum(axis=1)
            weighted=q['weight'][:,None]*l
            volume[pi]-=weighted.T@(div+theta[indices])
            scale[pi]+=weighted.T@(abs(div)+abs(theta[indices]))
            if projected:
                if np.any(abs(l@nodal_theta[pi]-theta[indices])>1e-11*(1+abs(theta[indices]))):
                    raise ValueError('projected nodal and GP volume histories disagree')
                weak_pressure[pi]+=weighted.T@(mean[indices]-gp_pressure[indices])
                pressure_scale[pi]+=weighted.T@np.maximum(np.maximum(abs(mean[indices]),abs(gp_pressure[indices])),1.)
            measure+=q['weight'].sum()
        cursor+=nq
    if cursor!=ng:raise ValueError('mixed quadrature inventory too long')
    mask=active[owner];mean=stress[:,:3].mean(axis=1)
    parameters={m['id']:m.get('mcc') for m in result['ground']['materials']}
    if any(p is None for p in parameters.values()):raise ValueError('mixed MCC material model mismatch')
    critical=np.array([parameters[material[e]]['critical_state_m'] for e in owner])
    dev=stress.copy();dev[:,:3]-=mean[:,None]
    q2=1.5*np.square(dev[:,:3]).sum(axis=1)+3*np.square(dev[:,3:]).sum(axis=1)
    if np.any(mean<=0) or not np.isfinite(critical).all() or np.any(critical<=0):raise ValueError('inadmissible mixed mean pressure or critical-state ratio')
    normalized_yield=q2/(critical*pc)**2+(mean/pc)*(mean/pc-1)
    if np.any(normalized_yield>1e-8):raise ValueError('mixed stress lies outside the MCC yield surface')
    pressure=float(np.max(abs(gp_pressure[mask]-mean[mask])/np.maximum(np.maximum(abs(gp_pressure[mask]),abs(mean[mask])),1.)))
    if projected:
        gap=float(abs(gp_pressure[mask]-mean[mask]).max())
        cached=projection['maximum_gp_pressure_projection_difference_kpa']
        if not np.isfinite(cached) or cached<0 or abs(cached-gap)>1e-11*(1+gap):raise ValueError('projected pressure-gap diagnostic differs')
        pressure=float(np.linalg.norm(weak_pressure)/max(np.linalg.norm(pressure_scale),measure*1e-12))
        volume+=offsets
    volume_relative=float(np.linalg.norm(volume)/max(np.linalg.norm(scale),measure*1e-12))
    tolerance=result['case']['settings']['relative_force_tolerance']
    if not np.isfinite(tolerance) or not 1e-11<=tolerance<=1e-5:raise ValueError('invalid mixed convergence tolerance')
    if not np.isfinite(pressure+volume_relative) or pressure>tolerance or volume_relative>tolerance:
        raise ValueError('mixed pressure consistency or weak volume equation failed')
    for name,actual in [('relative_volume_equation_residual',volume_relative),('relative_pressure_consistency_residual',pressure)]:
        stored=mixed[name]
        if not np.isfinite(stored) or stored<0 or abs(stored-actual)>1e-11:raise ValueError('mixed stored residual disagrees: '+name)
    material_strain=kinematic.copy();difference=theta-kinematic[:,:3].sum(axis=1)
    material_strain[:,:3]+=difference[:,None]/3
    maximum_material=float(abs(material_strain[mask]).max())
    cached=frame.get('maximum_absolute_material_engineering_strain_from_seed')
    if cached is not None and (not isinstance(cached,(int,float)) or isinstance(cached,bool)
        or not np.isfinite(cached) or abs(cached-maximum_material)>1e-11):
        raise ValueError('mixed cached maximum material strain disagrees with independent reconstruction')
    raw=dict(primary_mean_effective_stress_compression_positive_kpa=gp_pressure,
        material_mean_effective_stress_compression_positive_kpa=mean,
        cumulative_material_volume_strain_compression_positive=theta,
        cumulative_kinematic_volume_strain_compression_positive=kinematic[:,:3].sum(axis=1),
        material_minus_kinematic_volume_strain=difference,
        cumulative_material_engineering_strain_compression_positive=material_strain,
        cumulative_kinematic_engineering_strain_compression_positive=kinematic,
        reference_void_ratio=reference,last_material_partition_count=partitions)
    if projected:
        work_stress=stress.copy();work_stress[:,:3]+=(gp_pressure-mean)[:,None]
        raw['work_conjugate_effective_stress_compression_positive_kpa']=work_stress
    out=dict(raw=raw,keys=keys,primary=primary,owner=owner,active=active,
                metrics=dict(relative_volume_equation_residual=volume_relative,relative_pressure_consistency_residual=pressure,
                    maximum_material_engineering_strain_from_seed=maximum_material))
    if projected:
        out.update(projected_volume_dofs=nodal_theta,projected_compatibility_offsets=offsets)
        out['metrics'].update(projected_force_audit(mesh,result,frame,work_stress))
    return out


def projected_force_audit(mesh,result,frame,work_stress):
    state=frame['state'];active=np.asarray(state['active_elements']);internal=np.zeros((len(mesh['nodes']),3));cursor=0
    for e in range(len(mesh['tetrahedra'])):
        q=ref.volume(mesh,e);nq=len(q['weight'])
        if active[e]:
            local=-np.einsum('qai,qa,q->i',q['b'],work_stress[cursor:cursor+nq],q['weight']).reshape(-1,3)
            np.add.at(internal,q['nodes'],local)
        cursor+=nq
    target=array(state['applied_nodal_load_kn'],(len(mesh['nodes'])*3,),'applied nodal force').reshape(-1,3)
    case=result['case']
    if case['settings']['kinematics'] not in ('explicit','side_rollers'):
        raise ValueError('projected force export audit currently requires independent-node kinematics')
    case={**case,'boundary':case.get('boundary') or {'constraints':[]}}
    bases,error=ref.staged_constraint_space(mesh,case,state)
    residual=internal-target;free=[];support=np.zeros(3)
    for node,basis in bases.items():
        reduced=basis.T@residual[node];free.extend(reduced);support+=residual[node]-basis@reduced
    actual=float(np.linalg.norm(free)/max(np.linalg.norm(target),1.))
    tolerance=case['settings']['relative_force_tolerance']
    if actual>tolerance or error>1e-10 or abs(actual-frame['relative_force_residual'])>1e-11:
        raise ValueError('projected work-conjugate force or cached residual disagrees')
    if not np.allclose(support,frame['support_reaction_kn'],rtol=1e-11,atol=1e-8):
        raise ValueError('projected work-conjugate support reactions disagree')
    if not np.allclose(target.sum(axis=0),frame['applied_force_kn'],rtol=1e-11,atol=1e-8):
        raise ValueError('projected applied force summary disagrees')
    return dict(independent_work_conjugate_force_residual=actual,
        support_reaction_maximum_error_kn=float(abs(support-frame['support_reaction_kn']).max()))


def projected_offsets(mesh,state,active,keys):
    lookup=dict(zip(keys,range(len(keys))));out=np.zeros(len(keys));cursor=0
    u=np.asarray(state['displacement_m']);theta=np.asarray(state['mixed_mcc']['cumulative_volume_strain_compression_positive'])
    for e,tet in enumerate(mesh['tetrahedra']):
        q=ref.volume(mesh,e);nq=len(q['weight'])
        if active[e]:
            div=np.einsum('qai,i->qa',q['b'],u[q['nodes']].reshape(-1))[:,:3].sum(axis=1)
            l=np.column_stack((1-q['natural'].sum(axis=1),q['natural']))
            pi=[lookup[(n,mesh['material_ids'][e])] for n in tet]
            out[pi]+=(q['weight'][:,None]*l).T@(div+theta[cursor:cursor+nq])
        cursor+=nq
    return out


def validate(mesh,result):
    projected=result.get('schema')==PROJECTED_SCHEMA
    if result.get('schema') not in (SCHEMA,PROJECTED_SCHEMA) or not result.get('complete') or not result.get('frames') or result.get('failure_reason') is not None:
        raise ValueError('completed mixed MCC result required')
    if result['case']['schema']!=('soil-fem-projected-mcc-stages/1' if projected else 'soil-fem-mixed-mcc-stages/1') or not mesh.get('tetrahedra10'):
        raise ValueError('mixed MCC Tet10 input/result schema mismatch')
    if len(result['frames'])!=len(result['case']['stages'])+1 or any(s['kind'] not in ('boundary','pressure','excavate','reset_display') for s in result['case']['stages']):
        raise ValueError('unsupported or incomplete mixed construction stages')
    if mesh['units']!='m-kN-kPa' or result['ground']['units']!='m-kN-kPa' or result['ground']['initialization']['water_level_m'] is not None:
        raise ValueError('mixed dry material/mesh units or water policy mismatch')
    cells=np.asarray(mesh['tetrahedra10']);tet=np.asarray(mesh['tetrahedra'])
    if cells.shape!=(len(tet),10) or not np.array_equal(cells[:,:4],tet):raise ValueError('mixed Tet10 corners differ')
    previous=None
    method=result['frames'][0]['state'].get('mixed_mcc',{}).get('tangent_method',LEGACY_TANGENT)
    for i,frame in enumerate(result['frames']):
        fields=frame_fields(mesh,result,frame,check_geometry=i==0)
        if projected:
            if previous is None:
                expected=np.zeros(len(fields['keys']))
            elif previous['active_elements']!=frame['state']['active_elements']:
                expected=projected_offsets(mesh,previous,fields['active'],fields['keys'])
            else:
                old={(d['node_index'],d['material_id']):d['integral_m3'] for d in previous['mixed_mcc']['projected_volume']['compatibility_offsets']}
                expected=np.asarray([old[k] for k in fields['keys']])
            if not np.allclose(expected,fields['projected_compatibility_offsets'],rtol=1e-11,atol=1e-12):
                raise ValueError('frozen projected compatibility offset differs from stage origin')
        cfg=result['ground']['initialization'].get('geostatic')
        if i==0 and cfg is not None:
            maximum=max(fields['metrics']['maximum_material_engineering_strain_from_seed'],float(abs(fields['raw']['cumulative_kinematic_engineering_strain_compression_positive'][fields['active'][fields['owner']]]).max()))
            if maximum>cfg['maximum_initial_engineering_strain']:
                raise ValueError('mixed initial material/kinematic strain exceeds the configured domain')
        if frame['state']['mixed_mcc'].get('tangent_method',LEGACY_TANGENT)!=method:
            raise ValueError('mixed MCC tangent method changed between frames')
        if previous is not None:
            for point,e in enumerate(fields['owner']):
                if not fields['active'][e]:
                    for name in ('material_states',):
                        if frame['state'][name][point]!=previous[name][point]:raise ValueError('removed mixed material history changed')
                    for name in ('cumulative_volume_strain_compression_positive','reference_void_ratio','last_increment_partitions'):
                        if frame['state']['mixed_mcc'][name][point]!=previous['mixed_mcc'][name][point]:
                            raise ValueError('removed mixed primary history changed')
        previous=frame['state']
    if result.get('last_committed')!=previous:
        raise ValueError('mixed final committed state differs from the final frame')


def write_primary(output,index,mesh,fields):
    node=np.array([n for n,m in fields['keys']],dtype=int);material=np.array([m for n,m in fields['keys']],dtype=int)
    coordinates=np.asarray(mesh['nodes'])[node];p=fields['primary']
    filename=f'mean-pressure-dofs-{index:03}.vtu'
    values=dict(original_node_id=node,material_id=material,mean_effective_stress_compression_positive_kpa=p)
    extra=[]
    if 'projected_volume_dofs' in fields:
        values.update(cumulative_material_volume_strain_compression_positive=fields['projected_volume_dofs'],
            frozen_compatibility_offset_m3=fields['projected_compatibility_offsets'])
        extra=['cumulative_material_volume_strain_compression_positive','frozen_compatibility_offset_m3']
    meshio.write_points_cells(output/filename,coordinates,[('vertex',np.arange(len(node))[:,None])],point_data=values)
    checked=meshio.read(output/filename)
    assert np.array_equal(checked.points,coordinates)
    for k,v in values.items():assert np.array_equal(checked.point_data[k],v)
    with (output/f'mean-pressure-dofs-{index:03}.csv').open('x',encoding='utf-8-sig',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['original_node_id','material_id','x_m','y_m','z_m','mean_effective_stress_compression_positive_kpa']+extra)
        writer.writerows([[int(node[i]),int(material[i]),*coordinates[i],p[i],*[values[name][i] for name in extra]] for i in range(len(p))])
    with (output/f'mean-pressure-dofs-{index:03}.csv').open(encoding='utf-8-sig',newline='') as stream:
        rows=list(csv.DictReader(stream))
    assert [(int(r['original_node_id']),int(r['material_id'])) for r in rows]==fields['keys']
    assert np.array_equal([float(r['mean_effective_stress_compression_positive_kpa']) for r in rows],p)
    for name in extra:assert np.array_equal([float(r[name]) for r in rows],values[name])
    return filename


def write_history(output,index,frame,fields):
    filename=f'material-history-{index:03}.csv'
    state=frame['state'];raw=fields['raw'];columns=['integration_point_id','original_element_id','active',
        'primary_mean_effective_stress_kpa','material_mean_effective_stress_kpa','reference_void_ratio','void_ratio','pc_kpa',
        'cumulative_material_volume_strain','cumulative_kinematic_volume_strain','last_material_partition_count']
    components=('xx','yy','zz','xy','yz','zx')
    columns.extend('effective_stress_'+c+'_kpa' for c in components)
    projected='projected_volume_dofs' in fields
    if projected:columns.extend('work_conjugate_effective_stress_'+c+'_kpa' for c in components)
    columns.extend('plastic_engineering_strain_'+c for c in components)
    with (output/filename).open('x',encoding='utf-8-sig',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(columns)
        for i,e in enumerate(fields['owner']):
            work=[]
            if projected:work=raw['work_conjugate_effective_stress_compression_positive_kpa'][i] if fields['active'][e] else ['']*6
            writer.writerow([i,int(e),int(fields['active'][e]),
                raw['primary_mean_effective_stress_compression_positive_kpa'][i] if fields['active'][e] else '',
                raw['material_mean_effective_stress_compression_positive_kpa'][i],raw['reference_void_ratio'][i],state['material_states'][i]['void_ratio'],
                state['material_states'][i]['preconsolidation_pressure_kpa'],raw['cumulative_material_volume_strain_compression_positive'][i],
                raw['cumulative_kinematic_volume_strain_compression_positive'][i],raw['last_material_partition_count'][i],
                *state['material_states'][i]['effective_stress_compression_positive_kpa'],
                *work,
                *state['material_states'][i]['plastic_strain_compression_positive']])
    with (output/filename).open(encoding='utf-8-sig',newline='') as stream:rows=list(csv.DictReader(stream))
    assert len(rows)==len(fields['owner'])
    assert np.array_equal([float(r['cumulative_material_volume_strain']) for r in rows],raw['cumulative_material_volume_strain_compression_positive'])
    assert np.array_equal([float(r['pc_kpa']) for r in rows],[s['preconsolidation_pressure_kpa'] for s in state['material_states']])
    assert np.array_equal([float(r['reference_void_ratio']) for r in rows],raw['reference_void_ratio'])
    assert np.array_equal([[float(r['effective_stress_'+c+'_kpa']) for c in components] for r in rows],
        [s['effective_stress_compression_positive_kpa'] for s in state['material_states']])
    assert np.array_equal([[float(r['plastic_engineering_strain_'+c]) for c in components] for r in rows],
        [s['plastic_strain_compression_positive'] for s in state['material_states']])


def write_collection(output,frames):
    root=ET.Element('VTKFile',type='Collection',version='0.1',byte_order='LittleEndian');data=ET.SubElement(root,'Collection')
    for i,frame in enumerate(frames):ET.SubElement(data,'DataSet',timestep=str(i),group=frame['name'],part='0',file=f'mean-pressure-dofs-{i:03}.vtu')
    ET.ElementTree(root).write(output/'mean-pressure-dofs.pvd',encoding='utf-8',xml_declaration=True)
