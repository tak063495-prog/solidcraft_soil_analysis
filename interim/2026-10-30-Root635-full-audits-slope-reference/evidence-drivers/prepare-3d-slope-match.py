"""Prepare conforming Tet10 extrusions of the verified pinned2D slope."""
from pathlib import Path
import collections,hashlib,json,math,sys
import numpy as np
p=Path(__file__).resolve().parent
root=p.parents[1]
out=p/'plane-strain-slope-comparison-v1'
reference=out/'reference-2d-v2'
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
    with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
    with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
if (p/'pause-requested').exists():raise SystemExit(2)
verified=read(reference/'terminal.json');assert verified['passed']
binding=read(reference/'source-binding.json')['files_sha256']
assert all(sha(Path(q))==h for q,h in binding.items())
fixture=read(reference/'fixture-2d.json')
raw=np.asarray(fixture['nodes'],float)
assert len(raw)==661 and len(fixture['elements'])==200
assert np.max(np.abs(raw[:,5:7]))==0
source_nodes={int(row[0]):row for row in raw}
edge_owners=collections.Counter()
for element in fixture['elements']:
    v=element[1:5]
    for i in range(4):edge_owners[tuple(sorted((v[i],v[(i+1)%4])))]+=1
assert set(edge_owners.values())=={1,2}
fixed_segments=[[],[]]
for (a,b),count in edge_owners.items():
    if count!=1:continue
    for axis in range(2):
        if source_nodes[a][3+axis] and source_nodes[b][3+axis]:
            fixed_segments[axis].append((source_nodes[a][1:3],source_nodes[b][1:3]))
def fixed_mask(xz,axis):
    found=np.zeros(len(xz),bool)
    for a,b in fixed_segments[axis]:
        a=np.asarray(a);v=np.asarray(b)-a
        t=(xz-a)@v/(v@v)
        distance=np.linalg.norm(xz-a-t[:,None]*v,axis=1)
        found|=(t>=-1e-10)&(t<=1+1e-10)&(distance<1e-10)
    return found
assert all(np.array_equal(fixed_mask(raw[:,1:3],axis),raw[:,3+axis].astype(bool)) for axis in range(2))
protocol=dict(prepared=True,reference_fixed_factor=verified['fixed_factor_for_future3d_match'],
    reference_2d_source_sha256=sha(reference/'source-binding.json'),
    source_2d_Q8_four_point_hourglass_factor=.05,target_3d_Tet10_full_positive27_point_quadrature=True,
    geometry_exact_straight_2d_polylines_extruded_one_metre=True,
    all_transverse_displacements_fixed_same_continuum_plane_strain_problem=True,
    transverse_shear_and_Y_variation_measured_as_discretization_error_not_silently_assumed_zero=True,
    initial_stress_free_then_full_gravity_turn_on_no_K0_substitution=True,
    c_and_tan_phi_reduced_once_by_archived_FSS_psi0_and_no_extra_Rankine=True,
    compare_geostatic_total_displacement_from_seed_not_zeroed_construction_display=True,
    compare_all661_2d_nodes_front_middle_back_and_all200_original_quad_prism_mean_stresses=True,
    predeclared_acceptance=dict(relative_physical_force=1e-8,relative_global_force=1e-8,
        maximum_normalized_3d_yield=1e-7,displacement_sample_relative_L2=.05,
        maximum_displacement_error_over_reference_max=.1,prism_mean_stress_relative_L2=.05,
        transverse_shear_over_stress_L2=.01,depth_displacement_variation_over_max=.01,
        maximum_engineering_strain_from_seed=.05),
    no_tolerance_relaxation_or_Fs_upper_bound_from_numerical_failure=True,
    general_MC_3d_IDH_all_interfaces_SRM_and_allP6_not_proven_by_these_cases=True,p6_complete=False)
dump(out/'protocol.json',protocol)
edges=[(0,1),(1,2),(2,0),(0,3),(1,3),(2,3)]
cases=[]
for subdivisions,ny in [(1,1),(1,2),(2,2)]:
    folder=out/f'xy{subdivisions}-ny{ny}'
    folder.mkdir()
    points=[];lookup={};triangles=[];parents=[]
    def point(value):
        key=tuple(round(float(a),12) for a in value)
        if key not in lookup:lookup[key]=len(points);points.append(np.asarray(value,float))
        index=lookup[key]
        assert np.max(np.abs(points[index]-value))<1e-10
        return index
    for element in fixture['elements']:
        corner=np.asarray([source_nodes[k][1:3] for k in element[1:5]])
        def mapped(i,j):
            r=i/subdivisions;s=j/subdivisions
            return point((1-r)*(1-s)*corner[0]+r*(1-s)*corner[1]+r*s*corner[2]+(1-r)*s*corner[3])
        for i in range(subdivisions):
            for j in range(subdivisions):
                a,b,c,d=[mapped(*v) for v in [(i,j),(i+1,j),(i+1,j+1),(i,j+1)]]
                triangles.extend([[a,b,c],[a,c,d]]);parents.extend([element[0],element[0]])
    n=len(points)
    nodes=[[float(x),layer/ny,float(z)] for layer in range(ny+1) for x,z in points]
    tetrahedra=[];original_quads=[]
    for layer in range(ny):
        for triangle,parent in zip(triangles,parents,strict=True):
            a,b,c=sorted(triangle)
            a0,b0,c0=[v+layer*n for v in [a,b,c]]
            a1,b1,c1=[v+(layer+1)*n for v in [a,b,c]]
            for tet in [[a0,b0,c0,c1],[a0,b0,b1,c1],[a0,a1,b1,c1]]:
                xyz=np.asarray([nodes[k] for k in tet])
                determinant=float(np.linalg.det((xyz[1:]-xyz[0]).T))
                if determinant<0:tet[1],tet[2]=tet[2],tet[1]
                assert abs(determinant)>1e-10
                tetrahedra.append(tet);original_quads.append(parent)
    midpoint={};quadratic=[]
    for tet in tetrahedra:
        mids=[]
        for a,b in edges:
            key=tuple(sorted((tet[a],tet[b])))
            if key not in midpoint:
                midpoint[key]=len(nodes);nodes.append(((np.asarray(nodes[key[0]])+nodes[key[1]])/2).tolist())
            mids.append(midpoint[key])
        quadratic.append(tet+mids)
    xyz=np.asarray(nodes)
    ownership=collections.Counter(tuple(sorted(tet[i] for i in face)) for tet in tetrahedra
        for face in [(0,1,2),(0,1,3),(0,2,3),(1,2,3)])
    assert set(ownership.values())=={1,2}
    volumes=np.asarray([np.linalg.det((xyz[tet[1:]]-xyz[tet[0]]).T)/6 for tet in np.asarray(tetrahedra)])
    assert np.all(volumes>0) and abs(float(volumes.sum())-220.)<1e-9
    group_volume=np.bincount(original_quads,weights=volumes,minlength=201)[1:]
    reference_area=np.asarray([sum(fixture['checks']['gp_coordinates'][4*k+j][4] for j in range(4)) for k in range(200)])
    assert np.allclose(group_volume,reference_area,rtol=0,atol=1e-12)
    for cell in quadratic:
        for k,(a,b) in enumerate(edges):assert np.array_equal(xyz[cell[4+k]],(xyz[cell[a]]+xyz[cell[b]])/2)
    xz=xyz[:,[0,2]]
    node_sets=dict(all=list(range(len(nodes))),fixed_x=np.flatnonzero(fixed_mask(xz,0)).tolist(),
        fixed_z=np.flatnonzero(fixed_mask(xz,1)).tolist(),crest=np.flatnonzero(abs(xyz[:,2]-10)<1e-10).tolist())
    coordinate_map={tuple(round(float(v),12) for v in row):k for k,row in enumerate(xyz)}
    assert len(coordinate_map)==len(nodes),'coincident distinct topological nodes'
    sample_maps={}
    for y in [0.,.5,1.]:
        indices=[]
        for row in raw:
            key=tuple(round(float(v),12) for v in [row[1],y,row[2]])
            index=coordinate_map[key]
            assert np.max(np.abs(xyz[index]-[row[1],y,row[2]]))<1e-10
            indices.append(index)
        sample_maps[str(y)]=indices
    mesh=dict(schema='soil-fem-mesh/3',units='m-kN-kPa',axis='z_up',crs='LOCAL',origin_m=[0.,0.,0.],
        nodes=nodes,tetrahedra=tetrahedra,tetrahedra10=quadratic,quadrature_order=3,
        material_ids=[1]*len(tetrahedra),region_ids=[1]*len(tetrahedra),node_sets=node_sets,
        regions=[dict(id=1,body='GriffithsSlope',material_id=1,bounds_m=[xyz.min(axis=0).tolist(),xyz.max(axis=0).tolist()],
            volume_m3=float(volumes.sum()),source_sha256=sha(reference/'fixture-2d.json'))],excavation_sets={})
    material=verified['material']
    ground=dict(schema='soil-fem-ground/1',units='m-kN-kPa',initialization=dict(method='equilibrate',
        surface_surcharge_kpa=0.,water_level_m=None,water_unit_weight_kn_m3=9.80665,relative_equilibrium_tolerance=1e-8,
        geostatic=dict(seed=dict(kind='body_states',states={'GriffithsSlope':dict(effective_stress_compression_positive_kpa=[0.]*6,
            plastic_strain_compression_positive=[0.]*6,void_ratio=.8,equivalent_plastic_strain=0.,plastic_work_kpa=0.,last_return='elastic')}),
            initial_boundary_factor=0.,maximum_initial_engineering_strain=.05,relative_volume_tolerance=1e-8)),
        materials=[dict(id=1,name='Pinned2D MC same fixed factor',unit_weight_kn_m3=material['unit_weight_kn_m3'],
            mohr_coulomb={k:v for k,v in material.items() if k!='unit_weight_kn_m3'})],bodies=[dict(body='GriffithsSlope',material_id=1)])
    case=read(root/'examples/stages-mc-column.json')
    case['settings'].update(kinematics='explicit',relative_force_tolerance=1e-8,relative_linear_tolerance=1e-10,
        initial_increment=1/32,minimum_increment=1/4096,maximum_newton_iterations=35,material_maximum_substep=1e-4)
    case['boundary']=dict(constraints=[dict(node_set=group,direction=np.eye(3)[axis].tolist(),value_m=0.)
        for group,axis in [('all',1),('fixed_x',0),('fixed_z',2)]])
    case['monitors']={};case['monitor_sets']={'crest':'crest'}
    case['verification_reject_once_stages']=[]
    case['stages']=[dict(kind='boundary',name='retain_full_gravity_state',target_factor=0.,initial_increment=1.)]
    dump(folder/'mesh.json',mesh);dump(folder/'ground.json',ground);dump(folder/'stages.json',case)
    dump(folder/'reference-map.json',dict(original_quad_per_tet=original_quads,source2d_node_ids=raw[:,0].astype(int).tolist(),
        matching_3d_nodes_front_middle_back=sample_maps,source_quad_area_m2=reference_area.tolist(),
        source_reference_result_sha256=sha(reference/'published-2d-result.json')))
    row=dict(case=folder.name,mesh=str(folder/'mesh.json'),elements=len(tetrahedra),nodes=len(nodes),integration_points=27*len(tetrahedra),
        volume_m3=float(volumes.sum()),minimum_determinant_m3=float(6*volumes.min()),
        every_straight_Tet10_midpoint_and_all_positive_affine_Jacobians_verified=True,
        all200_original_prism_volumes_equal_reference_areas_to1e_minus12=True,
        all661_source_nodes_mapped_front_middle_back=True,
        input_sha256={name:sha(folder/name) for name in ['mesh.json','ground.json','stages.json','reference-map.json']},
        native_validation_and_actual_FEM_solve_and_comparison_pending=True)
    dump(folder/'prepared.json',row);cases.append(row)
assert all(sha(Path(q))==h for q,h in binding.items())
dump(out/'input-preparation.json',dict(passed=True,reference2d_independently_verified=True,cases=cases,
    protocol_sha256=sha(out/'protocol.json'),preparer_sha256=sha(Path(__file__)),
    Root_original_native_and_solvers_not_modified=True,no_new_solver_job_created_or_solved=True,
    actual3D_comparison_not_yet_proven=True,p6_complete=False))
print(json.dumps(dict(prepared=True,cases=[dict(case=c['case'],elements=c['elements'],nodes=c['nodes']) for c in cases],FEM_pending=True)),flush=True)
