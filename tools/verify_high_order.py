"""Independent NumPy Tet10 geostatic and accepted-stage force verification."""
import hashlib
import json
from pathlib import Path
import meshio
import numpy as np
import fem_reference as ref
from high_order_io import point_state


def hashfile(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def prepare(mesh,layout=None):
    q=[ref.volume(mesh,ei) for ei in range(len(mesh['tetrahedra']))]
    offset=np.cumsum([0]+[len(v['weight']) for v in q])
    if layout is not None:
        assert len(layout)==offset[-1]
        for ei,v in enumerate(q):
            for qi in range(len(v['weight'])):
                p=layout[offset[ei]+qi]
                assert p['element_index']==ei and p['local_index']==qi
                assert np.allclose(p['position_m'],v['position'][qi],rtol=1e-12,atol=1e-12)
                assert np.allclose(p['shape_functions'],v['shape'][qi],rtol=1e-12,atol=1e-12)
                assert np.isclose(p['weight_m3'],v['weight'][qi],rtol=1e-11,atol=1e-13)
    return q,offset


def internal(mesh,ground,states,active,q,offset,material_ids=None,pore_pressure_kpa=None):
    out=np.zeros((len(mesh['nodes']),3));max_yield=-np.inf
    materials={m['id']:m for m in ground['materials']}
    for ei in np.flatnonzero(active):
        qp=q[ei];points=states[offset[ei]:offset[ei+1]]
        stress=np.array([p['effective_stress_compression_positive_kpa'] for p in points])
        e=np.array([p['void_ratio'] for p in points])
        mean=stress[:,:3].mean(axis=1);dev=stress.copy();dev[:,:3]-=mean[:,None]
        q2=1.5*(np.square(dev[:,:3]).sum(axis=1)+2*np.square(dev[:,3:]).sum(axis=1))
        material=materials[(material_ids or mesh['material_ids'])[ei]]
        if 'mcc' in material:
            pc=np.array([p['preconsolidation_pressure_kpa'] for p in points]);m2=material['mcc']['critical_state_m']**2
            yf=(q2+m2*mean*(mean-pc))/(m2*pc**2);assert np.all(mean>0) and np.all(pc>0)
        elif 'invariant_double_hardening' in material:
            from hardening_reference import audit
            yf=audit(material['invariant_double_hardening'],points)
        else:
            from mc_reference import normalized_yield
            yf=normalized_yield(material['mohr_coulomb'],stress)
            assert all(p['equivalent_plastic_strain']>=0 and p['plastic_work_kpa']>=-1e-9 for p in points)
        assert np.isfinite(stress).all() and np.all(e>0) and np.all(yf<=1e-8)
        max_yield=max(max_yield,float(yf.max()))
        total=stress.copy();h=ground['initialization']['water_level_m']
        if pore_pressure_kpa is not None:total[:,:3]+=np.asarray(pore_pressure_kpa)[offset[ei]:offset[ei+1],None]
        elif h is not None:total[:,:3]+=np.maximum(h-qp['position'][:,2],0)[:,None]*ground['initialization']['water_unit_weight_kn_m3']
        force=-np.einsum('qai,qa,q->i',qp['b'],total,qp['weight'])
        out[qp['nodes']]+=force.reshape(-1,3)
    return out,max_yield


def physical(mesh,ground,active,pressure,original_only,material_ids=None,joint_system=None):
    out=np.zeros((len(mesh['nodes']),3));materials={m['id']:m for m in ground['materials']}
    for ei in np.flatnonzero(active):
        # Different quadrature than the Rust solver, exact for body loading.
        qp=ref.volume(mesh,ei,independent=True)
        out[qp['nodes'],2]-=qp['weight']@qp['shape']*materials[(material_ids or mesh['material_ids'])[ei]]['unit_weight_kn_m3']
    cut_area=0.;loaded_area=0.;cut_water=np.zeros(3);top=np.max(np.asarray(mesh['nodes'])[:,2])
    cfg=ground.get('initialization',{}).get('geostatic');hydro=set()
    if cfg:
        hydro={tuple(sorted(f)) for name in cfg.get('hydrostatic_face_sets',[]) for f in mesh['face_sets'][name]}
    for face,owners in ref.topology(mesh).items():
        live=[i for i in owners if active[i]]
        if len(live)!=1:continue
        opposite=None if joint_system is None else joint_system.counterparts.get(face)
        if opposite is not None and active[opposite]:continue
        qp=ref.face(mesh,live[0],face);normal=qp['area_vector'];area=np.linalg.norm(normal,axis=1)
        up=(normal[:,2]>0)&(np.abs(normal[:,0])+np.abs(normal[:,1])<1e-8*area)
        cut=len(owners)==2 or opposite is not None
        if cut:cut_area+=area.sum()
        h=ground['initialization']['water_level_m'];water_face=np.full(len(area),cut or face in hydro) if cfg else (cut|up)
        extra=up & (np.abs(qp['position'][:,2]-top)<1e-8 if original_only else np.ones(len(area),bool))
        force=-(extra*pressure)[:,None]*normal
        out[qp['nodes']]+=qp['shape'].T@force
        if h is not None and np.any(water_face):
            assert np.all(water_face),'water face must be a consistent physical boundary'
            from hydrostatic_reference import face_force
            wet_force,_=face_force(np.asarray(mesh['nodes'])[qp['nodes']],h,ground['initialization']['water_unit_weight_kn_m3'])
            out[qp['nodes']]+=wet_force
            if cut:cut_water+=wet_force.sum(axis=0)
        loaded_area+=area[extra].sum()
    if cfg:out+=ref.conditions_load(mesh,cfg.get('base_loads',{}),active,skip_removed_surfaces=True)
    return out,{'cut_surface_area_m2':float(cut_area),'loaded_surface_area_m2':float(loaded_area),'cut_water_force_kn':cut_water.tolist()}


def constraints(mesh,active,kinematics):
    xyz=np.asarray(mesh['nodes']);used=np.unique(ref.cells(mesh)[active]);fixed=np.ones_like(xyz,dtype=bool);fixed[used]=False
    if kinematics=='oedometer':fixed[used,:2]=True
    else:
        for a in (0,1):
            fixed[used,a]=np.isclose(xyz[used,a],xyz[:,a].min(),atol=1e-8,rtol=0)|np.isclose(xyz[used,a],xyz[:,a].max(),atol=1e-8,rtol=0)
    fixed[used,2]=np.isclose(xyz[used,2],xyz[:,2].min(),atol=1e-8,rtol=0)
    return used,fixed


def reduced_residual(xyz,residual,fixed,used,kinematics):
    if kinematics!='oedometer':return residual[~fixed]
    heights=np.copysign(np.floor(np.abs(xyz[:,2])*1e8+.5),xyz[:,2]);values={}
    for n in used:
        if not fixed[n,2]:values[heights[n]]=values.get(heights[n],0)+residual[n,2]
    return np.asarray(list(values.values()))


def verify_ground(run,solver=None):
    mesh=json.loads((run/'mesh.json').read_text());ground=json.loads((run/'ground.json').read_text());state=json.loads((run/'initial-state.json').read_text())
    manifest=json.loads((run/'bodies.json').read_text());provenance=json.loads((run/'mesh.provenance.json').read_text())
    assert hashfile(run/'ground.json')==manifest['ground_sha256']==provenance['ground_sha256']
    assert hashfile(run/'bodies.json')==provenance['cad_manifest_sha256']
    for r in manifest['regions']+manifest.get('excavation_tools',[]):assert hashfile(run/r['step'])==r['source_sha256']
    assert state['ground']==ground
    q,offset=prepare(mesh);points=state['point_states'];assert len(points)==offset[-1]
    materials={m['id']:m for m in ground['materials']};volumes=np.array([v['weight'].sum() for v in q]);maximum=0.
    for ei,qp in enumerate(q):
        expected,pore,material,layer=ref.initial_fields(mesh,ground,ei,qp['position'])
        for qi in range(len(qp['weight'])):
            p=points[offset[ei]+qi];history=point_state(p);s=np.array(history['effective_stress_compression_positive_kpa'])
            assert p['element_index']==ei and p['local_index']==qi and p['material_id']==mesh['material_ids'][ei]
            assert np.allclose(p['centroid_m'],qp['position'][qi],rtol=1e-12,atol=1e-12)
            assert np.isclose(p['weight_m3'],qp['weight'][qi],rtol=1e-11,atol=1e-13)
            maximum=max(maximum,float(np.max(np.abs(s-expected[qi]))));assert np.allclose(s,expected[qi],rtol=1e-12,atol=1e-10)
            assert np.isclose(p['pore_pressure_kpa'],pore[qi],rtol=1e-12,atol=1e-10)
            if 'mcc' in material:
                mean=s[:3].mean();pc=history['preconsolidation_pressure_kpa'];params=material['mcc']
                e=params['specific_volume_ncl_at_reference']-params['lambda']*np.log(pc/params['reference_pressure_kpa'])+params['kappa']*np.log(pc/mean)-1
                assert np.isclose(pc/mean,p['mcc_pc_over_mean_effective'],rtol=1e-12)
            elif 'invariant_double_hardening' in material:
                params=material['invariant_double_hardening'];e=layer['void_ratio']['value'];mean=s[:3].mean()
                young=params['unloading_young_ref_kpa']*(history['cap_pressure_kpa']/params['reference_pressure_kpa'])**params['stiffness_exponent']
                assert np.isclose(p['bulk_modulus_kpa'],young/(3*(1-2*params['poisson'])))
                assert np.isclose(p['shear_modulus_kpa'],young/(2*(1+params['poisson'])))
                consolidation=layer['consolidation']
                cap=consolidation['ratio']*mean if consolidation['definition']=='hardening_cap_over_mean_effective' else consolidation['cap_pressure_kpa']
                assert np.isclose(cap,history['cap_pressure_kpa']) and history['initial_cap_pressure_kpa']==history['cap_pressure_kpa']
                assert history['last_return']=='elastic' and history['shear_plastic_strain']==history['cap_plastic_volumetric_strain']==history['plastic_work_kpa']==0
                assert p['mcc_pc_over_mean_effective'] is None and p['vertical_ocr'] is None
            else:
                params=material['mohr_coulomb'];e=layer['void_ratio']['value']
                assert p['mcc_pc_over_mean_effective'] is None and p['vertical_ocr'] is None
                assert np.isclose(p['bulk_modulus_kpa'],params['young_kpa']/(3*(1-2*params['poisson'])))
                assert np.isclose(p['shear_modulus_kpa'],params['young_kpa']/(2*(1+params['poisson'])))
                assert history['last_return']=='elastic' and history['plastic_work_kpa']==history['equivalent_plastic_strain']==0
            assert np.isclose(e,history['void_ratio'],rtol=1e-12) and not np.array(history['plastic_strain_compression_positive']).any()
    weight=0.
    for r in mesh['regions']:
        v=volumes[np.array(mesh['region_ids'])==r['id']].sum();assert np.isclose(v,r['volume_m3'],rtol=1e-10,atol=1e-10)
        weight+=v*materials[r['material_id']]['unit_weight_kn_m3']
    active=np.ones(len(q),bool)
    internal_force,_=internal(mesh,ground,[point_state(p) for p in points],active,q,offset)
    f,_=physical(mesh,ground,active,ground['initialization']['surface_surcharge_kpa'],True)
    xyz=np.array(mesh['nodes']);used,fixed=constraints(mesh,active,'side_rollers');relative=np.linalg.norm((f-internal_force)[~fixed])/max(np.linalg.norm(f),1)
    assert relative<1e-10
    reaction=np.where(fixed,internal_force-f,0);assert np.allclose(f.sum(axis=0),state['applied_force_kn'],rtol=1e-12,atol=1e-8)
    assert np.allclose(reaction.sum(axis=0),state['support_reaction_kn'],rtol=1e-12,atol=1e-8)
    interfaces=sum(len(owners)==2 and mesh['region_ids'][owners[0]]!=mesh['region_ids'][owners[1]] for owners in ref.topology(mesh).values())
    assert interfaces==state['shared_interface_triangles']
    vtk=meshio.read(run/'initial-state.vtu');raw=meshio.read(run/'initial-state-points.vtu')
    quadratic=bool(mesh.get('tetrahedra10'));kind='tetra10' if quadratic else 'tetra'
    assert np.array_equal(vtk.cells_dict[kind],ref.cells(mesh))
    assert np.allclose(raw.point_data['void_ratio'],[point_state(p)['void_ratio'] for p in points])
    report={'nodes':len(xyz),'tetrahedra':len(q),'integration_points':len(points),'shared_interface_triangles':interfaces,
            'body_volumes_m3':state['body_volumes_m3'],'total_weight_kn':float(weight),'maximum_analytic_stress_error_kpa':maximum,
            'true_relative_free_force_residual':float(relative),'relative_global_force_balance_error':float(np.linalg.norm((f+reaction).sum(axis=0))/weight),
            'vtu_roundtrip':'passed','raw_integration_point_roundtrip':'passed'}
    (run/'ground-verification.json').write_text(json.dumps(report,indent=2))
    root=Path(__file__).resolve().parents[1];files=sorted(p for p in run.iterdir() if p.is_file() and p.name!='accepted-run.json')
    package={'schema':'soil-fem-accepted-ground-run/1','files_sha256':{p.name:hashfile(p) for p in files},
             'implementation_sha256':{str(p.resolve()):hashfile(p) for p in [(solver or root/'target/release/soil-fem-3d.exe'),root/'Cargo.lock',root/'rust-toolchain.toml']}}
    (run/'accepted-run.json').write_text(json.dumps(package,indent=2));print(json.dumps(report,indent=2))


def verify_stages(mesh_path,ground_path,result_path,exports,solver=None):
    mesh=json.loads(mesh_path.read_text());ground=ref.normalize_ground(json.loads(ground_path.read_text()));run=json.loads(result_path.read_text())
    if run['schema']=='soil-fem-stage-result/13':
        from verify_joint_water import verify_package
        return verify_package(mesh_path,ground_path,result_path,exports,solver)
    assert run['complete'] and run['schema'] in ('soil-fem-stage-result/1','soil-fem-stage-result/2','soil-fem-stage-result/3','soil-fem-stage-result/4','soil-fem-stage-result/5','soil-fem-stage-result/6','soil-fem-stage-result/7','soil-fem-stage-result/8','soil-fem-stage-result/9','soil-fem-stage-result/10','soil-fem-stage-result/11','soil-fem-stage-result/12','soil-fem-stage-result/13','soil-fem-stage-result/14') and run['ground']==ground
    case=ref.normalize_stage_case(json.loads((result_path.parent/'stages.json').read_text()))
    assert case==run['case']
    q,offset=prepare(mesh,run['integration_points']);xyz=np.array(mesh['nodes']);reports=[]
    beam_system=None
    if case.get('beam_system'):
        from beam_reference import System,coupled_forces
        beam_system=System(mesh,case['beam_system'])
    anchor_system=None
    if case.get('anchor_system'):
        from anchor_reference import System as AnchorSystem,coupled_forces as anchor_forces
        anchor_system=AnchorSystem(mesh,case['anchor_system'])
    joint_system=None
    if case.get('joint_system'):
        from joint_reference import System as JointSystem,coupled_forces as joint_forces
        joint_system=JointSystem(mesh,case['joint_system'],run['joint_integration_points'])
    general=run.get('geostatic') is not None;reference=np.asarray(run['reference_displacement_m']) if general else np.zeros_like(xyz)
    if general:
        from verify_geostatic import verify_initial
        verify_initial(mesh,ground,run,q,offset)
    seepage_pressures=None;seepage_reports=[]
    if run['schema']=='soil-fem-stage-result/12':
        from darcy_reference import stage_pressures
        seepage_pressures,seepage_reports=stage_pressures(mesh,ground,run,q,offset)
    for index,frame in enumerate(run['frames']):
        state=frame['state'];active=np.array(state['active_elements']);assert len(state['material_states'])==offset[-1]
        if seepage_pressures is not None:
            raw=meshio.read(exports/f'integration-points-{index:03}.vtu');ids=raw.point_data['integration_point_id']
            expected_pore=seepage_pressures[index][ids]
            assert np.allclose(raw.point_data['pore_pressure_kpa'],expected_pore,rtol=1e-10,atol=1e-8)
            expected_total=np.asarray([state['material_states'][i]['effective_stress_compression_positive_kpa'] for i in ids])
            expected_total[:,:3]+=expected_pore[:,None]
            assert np.allclose(raw.point_data['total_stress_compression_positive_kpa'],expected_total,rtol=1e-10,atol=1e-8)
        if run['schema']=='soil-fem-stage-result/11':
            expected_level=ground['initialization']['water_level_m']
            if expected_level is None:expected_level=float(xyz[:,2].min())
            for stage in case['stages'][:index]:
                if stage['kind']=='hydrostatic_water':expected_level=stage['target_level_m']
            assert np.isclose(state['hydrostatic_water_level_m'],expected_level,rtol=0,atol=1e-12)
            raw=meshio.read(exports/f'integration-points-{index:03}.vtu');ids=raw.point_data['integration_point_id']
            independent_positions=np.concatenate([v['position'] for v in q])
            expected_pore=np.maximum(expected_level-independent_positions[ids,2],0)*ground['initialization']['water_unit_weight_kn_m3']
            assert np.allclose(raw.point_data['pore_pressure_kpa'],expected_pore,rtol=1e-12,atol=1e-11)
            expected_total=np.asarray([state['material_states'][i]['effective_stress_compression_positive_kpa'] for i in ids])
            expected_total[:,:3]+=expected_pore[:,None]
            assert np.allclose(raw.point_data['total_stress_compression_positive_kpa'],expected_total,rtol=1e-12,atol=1e-10)
        frame_ground=ground
        if 'hydrostatic_water_level_m' in state:
            frame_ground={**ground,'initialization':{**ground['initialization'],'water_level_m':state['hydrostatic_water_level_m']}}
        material_ids=state.get('element_material_ids');force,yield_value=internal(mesh,frame_ground,state['material_states'],active,q,offset,material_ids,state.get('pore_pressure_kpa'))
        f,geometry=physical(mesh,frame_ground,active,state['surface_pressure_kpa'],case.get('pressure_on_original_surface_only',False),material_ids,joint_system)
        if case.get('boundary') and state.get('boundary_factor',0)!=0:
            f+=state['boundary_factor']*ref.conditions_load(mesh,case['boundary'],active)
        beam=None
        if beam_system is not None:
            force,f,beam=coupled_forces(beam_system,state,force,f,frame['beam_results'])
        if anchor_system is not None:force,f=anchor_forces(anchor_system,state,force,f,frame['anchor_results'])
        if joint_system is not None:force,f=joint_forces(joint_system,state,force,f,frame['joint_results'])
        if joint_system is not None and frame.get('joint_support') and beam_system is None and anchor_system is None:
            independent_support=joint_system.rigid_support(case,state,frame['joint_support'])
            assert independent_support['unrestrained_modes']==0
        assert np.allclose(f.reshape(-1),state['applied_nodal_load_kn'],rtol=1e-10,atol=1e-8)
        assert np.linalg.norm(state['excavation_release_load_kn'])<1e-9
        assert np.linalg.norm(state.get('transition_release_load_kn',[]))<1e-9
        used,fixed=constraints(mesh,active,case['settings']['kinematics']);residual=f-force
        if case.get('boundary') is not None:
            basis,constraint_error=ref.staged_constraint_space(mesh,case,state)
            reduced=np.concatenate([b.T@residual[n] for n,b in basis.items()])
            assert constraint_error<1e-10 and abs(constraint_error-frame['maximum_constraint_error_m'])<1e-11
        else:reduced=reduced_residual(xyz,residual,fixed,used,case['settings']['kinematics'])
        relative=np.linalg.norm(reduced)/max(np.linalg.norm(f),1);assert relative<case['settings']['relative_force_tolerance']*1.01
        reaction=np.where(fixed,force-f,0);reaction[np.setdiff1d(np.arange(len(xyz)),used)]=0
        if case.get('boundary') is not None:
            reaction[:]=0
            for n,b in basis.items():
                value=force[n]-f[n];reaction[n]=value-b@(b.T@value)
        applied=f.sum(axis=0);support=reaction.sum(axis=0)
        if beam is not None:
            applied+=beam['external_force']-sum(beam['attached_load'].values(),np.zeros(3))
            support+=beam['support_force']
        assert np.allclose(support,frame['support_reaction_kn'],rtol=1e-9,atol=1e-7)
        assert np.allclose(applied,frame['applied_force_kn'],rtol=1e-11,atol=1e-8)
        assert np.isclose(yield_value,frame['maximum_normalized_yield'],rtol=1e-8,atol=1e-8)
        center=(xyz.min(axis=0)+xyz.max(axis=0))/2
        moment=np.cross(xyz-center,f+reaction).sum(axis=0)
        if anchor_system is not None or joint_system is not None:moment=np.cross(xyz,f+reaction).sum(axis=0)
        if beam is not None:
            attached_moment=sum((np.cross(xyz[n],v) for n,v in beam['attached_load'].items()),np.zeros(3))
            moment=np.cross(xyz,f+reaction).sum(axis=0)+beam['external_moment']-attached_moment+beam['support_moment']
            assert np.allclose(moment,frame['global_moment_balance_error_kn_m'],atol=1e-6,rtol=1e-8)
        elif anchor_system is not None or joint_system is not None:assert np.allclose(moment,frame['global_moment_balance_error_kn_m'],atol=1e-6,rtol=1e-8)
        moment_scale=max(np.linalg.norm(f,axis=1).sum()*np.linalg.norm(np.ptp(xyz,axis=0)),1.)
        relative_moment=float(np.linalg.norm(moment)/moment_scale)
        if case['settings']['kinematics']!='oedometer':assert relative_moment<1e-7,relative_moment
        strain=0.;u=np.array(state['displacement_m'])
        for ei in np.flatnonzero(active):strain=max(strain,float(np.max(np.abs(np.einsum('qai,i->qa',q[ei]['b'],u[q[ei]['nodes']].reshape(-1))))))
        if general:
            assert 'maximum_absolute_engineering_strain_from_k0' not in frame
            assert np.isclose(strain,frame['maximum_absolute_engineering_strain_from_seed'],rtol=1e-10,atol=1e-12)
            construction=max(float(np.max(np.abs(np.einsum('qai,i->qa',q[ei]['b'],(u-reference)[q[ei]['nodes']].reshape(-1))))) for ei in np.flatnonzero(active))
            assert np.isclose(construction,frame['maximum_absolute_engineering_strain_from_initial_equilibrium'],rtol=1e-10,atol=1e-12)
            from verify_geostatic import verify_monitors
            monitor_origin=run['frames'][max(index-1,0)]
            if index and case['stages'][index-1]['kind']=='birth':monitor_origin={'state':next(a for a in run['transition_audits'] if a['stage']==frame['name'])['prepared_state']}
            verify_monitors(mesh,case,frame,monitor_origin,np.asarray(state.get('display_reference_displacement_m') or reference))
        else:assert np.isclose(strain,frame['maximum_absolute_engineering_strain_from_k0'],rtol=1e-10,atol=1e-12)
        if index:
            previous=run['frames'][index-1]['state']
            history_before=previous
            if case['stages'][index-1]['kind']=='reset_display':
                from verify_construction import verify_display_reset
                verify_display_reset(run,index)
            if case['stages'][index-1]['kind'] in ('birth','material_set','reset_state','install_beams','remove_beams','install_anchors','remove_anchors'):
                from verify_construction import verify_transition
                history_before=verify_transition(mesh,ground,run,index,q,offset)
            for ei in np.flatnonzero(~active):
                assert state['material_states'][offset[ei]:offset[ei+1]]==history_before['material_states'][offset[ei]:offset[ei+1]]
            if case['stages'][index-1]['kind']=='excavate':
                removed=np.array(previous['active_elements'])&~active;removed_force,_=internal(mesh,frame_ground,previous['material_states'],removed,q,offset,previous.get('element_material_ids'),previous.get('pore_pressure_kpa'))
                expected_holding=np.array(previous['applied_nodal_load_kn']).reshape(-1,3)-f-removed_force
                if joint_system is not None:
                    retained={**previous,'active_elements':active.tolist()}
                    expected_holding-=joint_system.audit(previous)-joint_system.audit(retained)
                expected_holding[np.setdiff1d(np.arange(len(xyz)),used)]=0
                audit=next(a for a in run['excavation_audits'] if a['stage']==frame['name'])
                assert np.allclose(expected_holding.reshape(-1),audit['initial_holding_load_kn'],rtol=1e-9,atol=1e-7)
                assert audit['zero_release_relative_force_residual']<case['settings']['relative_force_tolerance']
        else:assert np.array_equal(u,reference)
        reports.append({'stage':frame['name'],'active_elements':int(active.sum()),'active_material_points':int(sum(len(q[i]['weight']) for i in np.flatnonzero(active))),
                        'independent_force_residual':float(relative),'force_balance_error_kn':(applied+support).tolist(),
                        'moment_balance_error_kn_m':moment.tolist(),'relative_moment_error':relative_moment,
                        'maximum_normalized_yield':yield_value,('maximum_absolute_engineering_strain_from_seed' if general else 'maximum_absolute_engineering_strain_from_k0'):strain,**geometry})
    assert all(a['committed_state_unchanged_on_rejection'] for a in run['attempts'] if not a['accepted'])
    export=json.loads((exports/'export-verification.json').read_text());assert export['mesh_sha256']==hashfile(mesh_path) and export['result_sha256']==hashfile(result_path)
    if mesh.get('tetrahedra10'):assert export['integration_point_roundtrip']=='passed'
    report={'schema':'soil-fem-stage-verification/2','frames':reports,'vtu_roundtrip':export['vtu_roundtrip'],'raw_integration_point_roundtrip':'passed'}
    if seepage_reports:report['independent_darcy_solutions']=seepage_reports
    if joint_system is not None:
        assert export['joint_point_roundtrip']=='passed'
        report['joint_geometry_force_energy_audit']='passed'
        if all(f.get('joint_support') for f in run['frames']) and beam_system is None and anchor_system is None:
            report['joint_rigid_support_svd_audit']='passed'
    (result_path.parent/'stage-verification.json').write_text(json.dumps(report,indent=2))
    root=Path(__file__).resolve().parents[1]
    assert run['provenance']['solver_binary_sha256']==hashfile((solver or root/'target/release/soil-fem-3d.exe')),'result produced by a different binary'
    package={'schema':'soil-fem-accepted-staged-run/1','inputs_and_results_sha256':{str(p.resolve()):hashfile(p) for p in [mesh_path,ground_path,result_path,result_path.parent/'stages.json']},
             'solver_binary_sha256':run['provenance']['solver_binary_sha256']}
    (result_path.parent/'accepted-stage-run.json').write_text(json.dumps(package,indent=2));print(json.dumps(report,indent=2))
