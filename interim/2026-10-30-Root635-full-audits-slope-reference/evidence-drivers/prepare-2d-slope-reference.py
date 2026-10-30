"""Verify pinned published 2D slope fields before a new independent3D match.

This reassembles the published2D output. It does not execute Excel, certify
an Fs upper bound, or claim that the new3D comparison has already passed.
"""
from pathlib import Path,PurePosixPath
import csv,datetime,hashlib,io,json,math,os,subprocess,sys,zipfile
import numpy as np
p=Path(__file__).resolve().parent
root=p.parents[1]
repo=root/'external/excel-vba-soil-fem'
d=repo/'docs/inconsistent-evidence24'
out=p/'plane-strain-slope-comparison-v1/reference-2d-v2'
out.mkdir(parents=True)
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
    with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
    with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
if (p/'pause-requested').exists():raise SystemExit(2)
dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,
    script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
revision='440d246d9e32f5892c79539d4a7b797ad29442ee'
git=lambda *args:subprocess.check_output(['git','-C',str(repo),*args],encoding='utf-8',
    env=dict(os.environ,GIT_OPTIONAL_LOCKS='0'),creationflags=subprocess.CREATE_NO_WINDOW).strip()
assert git('rev-parse','HEAD')==revision
files=['fixtures.json','srm_cases_INCO24.zip','archive_entries.json','manifest.json','material_force_audit.json']
binding={}
for name in files:
    q=d/name;raw=q.read_bytes()
    blob=hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
    assert blob==git('rev-parse',revision+':docs/inconsistent-evidence24/'+name)
    binding[str(q)]=sha(q)
binding[str(Path(__file__).resolve())]=sha(Path(__file__))
metadata=read(d/'manifest.json')
for name in files:
    if name=='manifest.json':continue
    assert sha(d/name)==metadata[name]['sha256'] and (d/name).stat().st_size==metadata[name]['bytes']
fixture=next(f for f in read(d/'fixtures.json')['fixtures'] if f['name']=='griffiths_1999_coarse')
archive=d/'srm_cases_INCO24.zip'
with zipfile.ZipFile(archive) as z:
    names=z.namelist();assert len(names)==len(set(names))==163
    assert z.testzip() is None
    assert all(not PurePosixPath(n).is_absolute() and '..' not in PurePosixPath(n).parts for n in names)
    entries=read(d/'archive_entries.json')['entries']
    assert set(names)==set(entries)|{'archive_entries.json'} and len(entries)==162
    assert z.read('archive_entries.json')==(d/'archive_entries.json').read_bytes()
    for name,item in entries.items():
        raw=z.read(name)
        assert len(raw)==item['bytes'] and hashlib.sha256(raw).hexdigest()==item['sha256'],name
    source='source/candidate24.xlsm'
    source_sha=hashlib.sha256(z.read(source)).hexdigest()
    assert source_sha=='83b1e4d527c0a9920603d052743df829abc552a582d319e29f20f8776a132856'
    result_name='cases/griffiths_1999_coarse/results/griffiths_1999_coarse_critical24.json'
    raw=z.read(result_name);result=json.loads(raw)
    (out/'published-2d-result.json').write_bytes(raw)
    strain_name='cases/griffiths_1999_coarse/analysis/griffiths_1999_coarse_critical24_strain.csv'
    raw=z.read(strain_name);strain_rows=list(csv.DictReader(io.StringIO(raw.decode('utf-8-sig'))))
    (out/'published-2d-strain.csv').write_bytes(raw)
assert result['status']=='PASS|' and result['build']=='20261010_INCO_24'
assert result['source_sha256'].lower()==source_sha
assert result['fixture_sha256'].lower()==sha(d/'fixtures.json')
assert result['flow_policy']=='INCONSISTENT' and result['hourglass_factor']==.05
metrics=dict(zip(result['metric_names'],result['metrics'],strict=True))
assert metrics['analysis_ok'] and metrics['nodes']==661 and metrics['elements']==200
factor=metrics['fss'];assert 1.34<factor<1.36
material=fixture['materials'][0]
assert material[:8]==[1,100000,.3,1,20,20,10,0]
young,nu,gamma=material[1],material[2],material[4]
cohesion=material[6]/factor
friction=math.atan(math.tan(math.radians(material[5]))/factor)
nodes=np.asarray(fixture['nodes'],float)
ids=nodes[:,0].astype(int);assert np.array_equal(ids,np.arange(1,662))
xy=nodes[:,1:3]
u=np.asarray(result['displacements_flat'],float).reshape(-1,3)
assert np.array_equal(u[:,0],ids)
displacement=u[:,1:]
stresses=np.asarray(result['stresses_flat'],float).reshape(-1,5)
states=np.asarray(result['gauss_state_flat'],float).reshape(-1,8)
assert len(stresses)==len(states)==len(strain_rows)==800
assert result['gauss_state_columns']==['element','gauss','sigma_z','yield_function','yielded','principal_max','principal_mid','principal_min']
assert np.array_equal(stresses[:,:2],states[:,:2])
assert np.isfinite(np.r_[displacement.ravel(),stresses.ravel(),states.ravel()]).all()
fixed=nodes[:,3:5].astype(bool)
assert np.max(np.abs(displacement[fixed]))==0

def shape(r,s):
    values=[.25*(1+a*r)*(1+b*s)*(a*r+b*s-1) for a,b in [(-1,-1),(1,-1),(1,1),(-1,1)]]
    return np.asarray(values+[.5*(1-r*r)*(1-s),.5*(1+r)*(1-s*s),.5*(1-r*r)*(1+s),.5*(1-r)*(1-s*s)])

g=1/math.sqrt(3)
quadrature=[(-g,-g),(g,-g),(g,g),(-g,g)]
internal=np.zeros((2,len(nodes),2))
external=np.zeros((len(nodes),2))
volume=0.;strain_error=0.;deformed_det_error=0.;point_error=0.
principal_error=0.;maximum_yield=-math.inf;maximum_strain=0.
projected_raw=[]
for k,element in enumerate(fixture['elements']):
    assert element[0]==k+1 and element[-1]==1
    index=np.asarray(element[1:9],int)-1
    assert len(set(index))==8
    xyz=xy[index]
    for local,(r,s) in enumerate(quadrature):
        row=4*k+local
        assert tuple(stresses[row,:2])==(element[0],local+1)
        n=shape(r,s)
        derivative=np.column_stack((np.imag(shape(r+1e-30j,s))/1e-30,np.imag(shape(r,s+1e-30j))/1e-30))
        jacobian=xyz.T@derivative
        weight=float(np.linalg.det(jacobian));assert weight>0
        grad=derivative@np.linalg.inv(jacobian)
        volume+=weight
        np.add.at(external,index,np.column_stack((np.zeros(8),-gamma*weight*n)))
        sx,sy,tau=stresses[row,2:]
        for option,sign in enumerate([1.,-1.]):
            traction=-weight*np.column_stack((grad[:,0]*sx+grad[:,1]*tau*sign,grad[:,1]*sy+grad[:,0]*tau*sign))
            np.add.at(internal[option],index,traction)
        gradient=displacement[index].T@grad
        engineering=np.array([gradient[0,0],gradient[1,1],gradient[0,1]+gradient[1,0]])
        target=strain_rows[row]
        assert (int(target['element']),int(target['gauss']))==(element[0],local+1)
        expected=np.array([float(target[q]) for q in ['ex','ey','engineering_shear']])
        strain_error=max(strain_error,float(np.max(np.abs(engineering-expected))))
        deformed_det_error=max(deformed_det_error,abs(float(np.linalg.det(np.eye(2)+gradient))-float(target['deformed_det'])))
        maximum_strain=max(maximum_strain,float(np.max(np.abs(engineering))))
        point=fixture['checks']['gp_coordinates'][row]
        assert point[:2]==[element[0],local+1]
        point_error=max(point_error,float(np.max(np.abs(np.r_[n@xyz,weight]-point[2:]))))
        tensor=np.array([[sx,0,tau],[0,states[row,2],0],[tau,0,sy]])
        eigen=np.linalg.eigvalsh(tensor)[::-1]
        principal_error=max(principal_error,float(np.max(np.abs(eigen-states[row,5:8]))))
        yield_function=eigen[0]-eigen[-1]-(eigen[0]+eigen[-1])*math.sin(friction)-2*cohesion*math.cos(friction)
        maximum_yield=max(maximum_yield,yield_function/max(1.,float(np.max(np.abs(eigen))),cohesion))
        projected_raw.append(dict(element=element[0],gauss=local+1,position_xz_m=(n@xyz).tolist(),
            weight_m2=weight,stress_x_z_shear_kpa=[float(sx),float(sy),float(tau)],
            out_of_plane_stress_kpa=float(states[row,2]),compression_positive_engineering_strain_x_z_xz=(-engineering).tolist()))
norm=float(np.linalg.norm(external[~fixed]));assert norm>0
residuals=[float(np.linalg.norm((force-external)[~fixed])/norm) for force in internal]
matches=[i for i,value in enumerate(residuals) if value<1e-8]
assert matches==[0],residuals
assert abs(volume-220.)<1e-9 and np.allclose(external.sum(axis=0),[0,-4400],rtol=0,atol=1e-9)
assert strain_error<1e-12 and deformed_det_error<1e-12 and point_error<1e-12
assert principal_error<1e-9 and maximum_yield<1e-8
reaction=np.where(fixed,internal[0]-external,0).sum(axis=0)
assert np.allclose(reaction,[metrics['reaction_x'],metrics['reaction_y']],rtol=0,atol=1e-7)
assert all(sha(Path(q))==h for q,h in binding.items())
dump(out/'fixture-2d.json',fixture)
dump(out/'reference-points.json',projected_raw)
dump(out/'source-binding.json',dict(pinned_revision=revision,files_sha256=binding,
    original_2d_workbook_sha256=source_sha,archive_members=163,all162_member_hashes_and_manifest_and_CRC_verified=True,
    published_actual_Excel_result_not_new_Excel_run=True,root_license='Apache-2.0'))
(out/'REFERENCE-LICENSE.txt').write_bytes((repo/'LICENSE').read_bytes())
dump(out/'terminal.json',dict(passed=True,pinned_revision=revision,elements2d=200,nodes2d=661,points2d=800,
    fixed_factor_for_future3d_match=factor,original_search_numeric_bracket_not_adopted_as_physical_collapse=True,
    material=dict(young_kpa=young,poisson=nu,unit_weight_kn_m3=gamma,cohesion_kpa=cohesion,
        friction_deg=math.degrees(friction),dilation_deg=0.,flow='non_associated',tension_cutoff_kpa=None),
    Q8_four_point_reassembled_material_relative_residual=residuals[0],
    wrong_shear_sign_counterexample_relative_residual=residuals[1],
    original_stress_and_raw_out_of_plane_and_all3D_principal_stresses_verified=True,
    maximum_full_3d_MC_normalized_yield=maximum_yield,principal_max_absolute_error_kpa=principal_error,
    archived_strain_and_deformation_determinant_all800_within_1e_minus12=True,
    strain_component_and_deformed_determinant_absolute_threshold=1e-12,
    maximum_strain_component_error=strain_error,maximum_deformed_determinant_error=deformed_det_error,
    all800_reference_positions_and_Jacobians_verified=True,maximum_reference_coordinate_or_det_error=point_error,
    area_m2=volume,body_force_per_unit_depth_kn=external.sum(axis=0).tolist(),reaction_per_unit_depth_kn=reaction.tolist(),
    maximum_absolute_engineering_strain=maximum_strain,
    original_reference_fields_and_sources_unchanged=True,
    future3d_Tet10_match_and_mesh_depth_sensitivity_and_realFs_search_and_generalP6_not_yet_proven=True,p6_complete=False))
print(json.dumps(dict(passed=True,points=800,full_force_relative_residual=residuals[0],
    wrong_shear_sign_residual=residuals[1],fixed_factor=factor,actual3d_comparison_pending=True)),flush=True)
