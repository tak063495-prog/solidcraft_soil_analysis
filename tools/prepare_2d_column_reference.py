"""Pinned published 2D workbook plus matching 2D/3D MC column fixtures.

The external sparse checkout is read only. Workbook/DLL bytes come from the
specified Git objects and are checked against their Git blob identities.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import zipfile
import numpy as np
import fem_reference as ref
from verify_incompressibility import mesh as cube


def write(path,value):
    with path.open('x',encoding='utf-8') as f:json.dump(value,f,indent=2)


def prepare(out):
    root=Path(__file__).resolve().parents[1];repo=root/'external/excel-vba-soil-fem'
    revision='440d246d9e32f5892c79539d4a7b797ad29442ee'
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip();assert head==revision
    out.mkdir(exist_ok=False);source=out/'source';source.mkdir();objects=[]
    for name in ('2DSoilFEM_ui_20261024.xlsm','soilfem_native_x64.dll','BUILD_RECORD.json'):
        path='workbook/ui_20261024/'+name;spec=revision+':'+path
        identity=subprocess.check_output(['git','rev-parse',spec],cwd=repo,text=True).strip()
        target=source/name
        with target.open('xb') as f:subprocess.run(['git','show',spec],cwd=repo,stdout=f,check=True)
        data=target.read_bytes();assert hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==identity
        objects.append(dict(path=path,git_blob=identity,sha256=hashlib.sha256(data).hexdigest(),bytes=len(data)))
    with zipfile.ZipFile(source/'2DSoilFEM_ui_20261024.xlsm') as z:
        assert 'xl/vbaProject.bin' in z.namelist() and z.testzip() is None
    assert (source/'soilfem_native_x64.dll').read_bytes().startswith(b'MZ')
    width=.2;height=.2;gamma=20.;young=1000.;nu=.25;cohesion=2.;phi=10.;load=12.
    for ny in (2,4):
        case_dir=out/f'ny{ny}';case_dir.mkdir();nodes=[];lookup={}
        for iy in range(2*ny+1):
            for ix in range(3):
                if ix%2 and iy%2:continue
                index=len(nodes)+1;lookup[(ix,iy)]=index
                nodes.append([index,width*ix/2,height*iy/(2*ny),1,int(iy==0),0.,0.,0.,0.])
        elements=[]
        for iy in range(ny):
            keys=[(0,2*iy),(2,2*iy),(2,2*iy+2),(0,2*iy+2),(1,2*iy),(2,2*iy+1),(1,2*iy+2),(0,2*iy+1)]
            elements.append([iy+1,*[lookup[k] for k in keys],1])
        loads=[[1,lookup[(ix,2*ny)],'Y','FORCE',-load*width*weight,None,None,None,None,'consistent top traction'] for ix,weight in ((0,1/6),(1,2/3),(2,1/6))]
        fixture=dict(materials=[[1,young,nu,1.,gamma,phi,cohesion,0.,0,1,'SOIL',0.,0]],nodes=nodes,elements=elements,
            loads=loads,stages=[[1,'GRAVITY',None,None,1,'elastic gravity reference'],[2,'LOAD',1,None,1,'plastic preload'],[3,'UNLOAD',1,None,1,'remove top load; retain gravity']],
            settings=dict(FLOW_POLICY='INCONSISTENT',OUTPUT_STAGE_MODE='ALL',OUTPUT_INACTIVE_ELEMENTS='SKIP',EXPORT_LOAD='',EXPORT_MODE='STAGE',VIEW_DETAIL_MODE='SIMPLE'))
        write(case_dir/'fixture-2d.json',fixture)
    m=cube(2,True);xyz=np.asarray(m['nodes'])*[width,1.,height];m['nodes']=xyz.tolist()
    m['node_sets']=dict(all=list(range(len(xyz))),bottom=np.flatnonzero(xyz[:,2]==0.).tolist(),top=np.flatnonzero(xyz[:,2]==height).tolist())
    m['face_sets']=dict(top=[list(f) for f,owners in ref.topology(m).items() if len(owners)==1 and np.all(xyz[list(f),2]==height)])
    m['regions']=[dict(id=1,body='Clay',material_id=1,bounds_m=[[0.,0.,0.],[width,1.,height]],volume_m3=width*height,
        source_sha256=hashlib.sha256(b'analytic common 2D/3D plane-strain column').hexdigest())]
    write(out/'mesh-3d.json',m)
    ground=dict(schema='soil-fem-ground/1',units='m-kN-kPa',initialization=dict(method='horizontal_k0',surface_surcharge_kpa=0.,water_level_m=None,
        water_unit_weight_kn_m3=9.80665,relative_equilibrium_tolerance=1e-9),
        materials=[dict(id=1,name='MC comparison',unit_weight_kn_m3=gamma,mohr_coulomb=dict(young_kpa=young,poisson=nu,cohesion_kpa=cohesion,
            friction_deg=phi,dilation_deg=0.,tension_cutoff_kpa=0.,flow='non_associated'))],
        layers=[dict(body='Clay',material_id=1,k0_x=nu/(1-nu),k0_y=nu/(1-nu),consolidation=dict(definition='not_applicable'),void_ratio=dict(mode='explicit',value=.8))])
    write(out/'ground-3d.json',ground)
    case=json.loads((root/'examples/stages-mc-column.json').read_text(encoding='utf-8'))
    case['settings'].update(kinematics='explicit',initial_increment=1/16,minimum_increment=1/1024,maximum_newton_iterations=35,material_maximum_substep=1e-4)
    case['boundary']=dict(constraints=[dict(node_set=group,direction=np.eye(3)[a].tolist(),value_m=0.) for group,a in (('all',0),('all',1),('bottom',2))],
                          surface_loads=[dict(kind='pressure',face_set='top',pressure_kpa=load)])
    case['monitors']={};case['monitor_sets']={'top':'top'};case['verification_reject_once_stages']=[]
    case['stages']=[dict(kind='boundary',name='preload',target_factor=1.),dict(kind='boundary',name='unload',target_factor=0.)]
    write(out/'stages-3d.json',case)
    write(out/'provenance.json',dict(revision=revision,objects=objects,width_m=width,height_m=height,extrusion_depth_m=1.,
        young_kpa=young,poisson=nu,cohesion_kpa=cohesion,friction_deg=phi,dilation_deg=0.,unit_weight_kn_m3=gamma,top_pressure_kpa=load,
        reference_convention='2D gravity displacement is retained separately; 3D K0 has the same elastic stress and zero reference displacement',
        scope='Common MC plane-strain column only; no general 2D/3D equivalence claim',
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
    print('Prepared pinned workbook/DLL and MC gravity/preload/unload column inputs',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);a=p.parse_args();prepare(a.output.resolve())
