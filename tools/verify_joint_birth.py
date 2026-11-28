"""Generate small contact-BIRTH cases and record independent acceptance evidence."""
import argparse
import copy
import json
import os
from pathlib import Path
import subprocess
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from joint_birth_reference import point_path,audit_elastic,audit_hydrostatic_reassembly

ROOT=Path(__file__).resolve().parents[1]
def read(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def dump(path,value):Path(path).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def without_clock(value):
    if isinstance(value,dict):return {k:without_clock(v) for k,v in value.items() if k!='elapsed_seconds'}
    if isinstance(value,list):return [without_clock(v) for v in value]
    return value

def fixture(mode='transfer_force',quadratic=True,water=None):
    mesh=read(ROOT/'examples'/('joint-column-tet10-mesh.json' if quadratic else 'joint-column-mesh.json'))
    ground=read(ROOT/'examples/ground-joint-mc.json');case=read(ROOT/'examples/stages-joint-normal-column.json')
    ground['initialization']['geostatic']['base_loads']={}
    ground['initialization']['surface_surcharge_kpa']=10.
    for s in ground['initialization']['geostatic']['seed']['states'].values():s['effective_stress_compression_positive_kpa']=[0.]*6
    for m in ground['materials']:m['unit_weight_kn_m3']=1e-20
    case['verification_reject_once_stages']=[]
    case['joint_system']['joints'][0]['birth_policy']=mode
    case['stages']=[dict(kind='pressure',name='preload',target_pressure_kpa=50.),dict(kind='joint_birth',name='contact',joints=['Seam']),dict(kind='pressure',name='reload',target_pressure_kpa=80.),dict(kind='pressure',name='unload',target_pressure_kpa=20.)]
    if water:
        mesh=read(ROOT/'examples/joint-water-column-tet10-mesh.json');ground=read(ROOT/'examples/ground-joint-water-mc.json');ground['initialization']['water_level_m']=1.5
        original=read(ROOT/'examples'/f'stages-joint-water-{water if water!="hydrostatic" else "column"}.json')
        case=original;case['joint_system']['joints'][0]['birth_policy']=mode
        birth=dict(kind='joint_birth',name='contact',joints=['Seam'])
        first=dict(kind='pressure',name='preload',target_pressure_kpa=50.) if water=='hydrostatic' else copy.deepcopy(original['stages'][0])
        if first['kind']=='steady_seepage':first['exterior_water_level_m']=1.5
        case['stages']=[first,birth,dict(kind='pressure',name='reload',target_pressure_kpa=80.)]
        case['verification_reject_once_stages']=[]
    return mesh,ground,case

def shear_fixture():
    mesh,ground,case=fixture('transfer_force',True)
    ground['initialization']['surface_surcharge_kpa']=0.
    xyz=np.asarray(mesh['nodes']);base=[];increment=[]
    sigma=np.array([[-100.,0.,12.],[0.,-100.,0.],[12.,0.,-100.]])
    extra=np.array([[0.,0.,10.],[0.,0.,0.],[10.,0.,0.]])
    faces=copy.deepcopy(mesh['face_sets']['exterior'])
    for axis in range(3):
        for sign,level in [(-1,xyz[:,axis].min()),(1,xyz[:,axis].max())]:
            selected=[f for f in faces if np.all(np.abs(xyz[f,axis]-level)<1e-12)]
            name=f'outer-{axis}-{sign}';mesh['face_sets'][name]=selected
            normal=np.eye(3)[axis]*sign
            base.append(dict(kind='traction',face_set=name,traction_kpa=(sigma@normal).tolist()))
            increment.append(dict(kind='traction',face_set=name,traction_kpa=(extra@normal).tolist()))
    ground['initialization']['geostatic']['base_loads']={'surface_loads':base}
    pins=[([0.,0.,0.],[0,1,2]),([1.,0.,0.],[1,2]),([0.,1.,0.],[2])];constraints=[]
    for i,(position,axes) in enumerate(pins):
        name=f'pin-{i}';node=int(np.argmin(np.linalg.norm(xyz-position,axis=1)));assert np.linalg.norm(xyz[node]-position)<1e-12
        mesh['node_sets'][name]=[node]
        constraints.extend(dict(node_set=name,direction=np.eye(3)[a].tolist(),value_m=0.) for a in axes)
    case['settings']['kinematics']='explicit';case['boundary']=dict(constraints=constraints,surface_loads=increment)
    case['settings']['relative_force_tolerance']=1e-11;case['settings']['relative_linear_tolerance']=1e-12
    ground['initialization']['relative_equilibrium_tolerance']=1e-12
    case['stages']=[dict(kind='joint_birth',name='contact',joints=['Seam']),dict(kind='boundary',name='shear_load',target_factor=1.),dict(kind='boundary',name='shear_unload',target_factor=0.)]
    return mesh,ground,case

def run_case(binary,folder,mesh,ground,case):
    folder.mkdir(parents=True,exist_ok=False)
    for name,data in [('mesh',mesh),('ground',ground),('case',case)]:dump(folder/f'{name}.json',data)
    inputs=[str(folder/f'{n}.json') for n in ('mesh','ground','case')]
    def invoke(args,log):
        result=subprocess.run([str(binary),*args],text=True,capture_output=True,encoding='utf-8')
        (folder/log).write_text(result.stdout+result.stderr,encoding='utf-8')
        assert result.returncode==0,(folder,result.returncode,result.stderr)
    invoke(['solve-stages',*inputs,str(folder/'run.json')],'run.log');run=read(folder/'run.json');assert run['complete']
    return run,inputs,invoke

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);parser.add_argument('--solver',type=Path,default=ROOT/'target/release/soil-fem-3d');parser.add_argument('--point-driver',type=Path,default=ROOT/'target/release/examples/joint_birth_points');args=parser.parse_args()
    output=args.output.resolve();output.mkdir(parents=True,exist_ok=False);solver=args.solver.resolve()
    report=dict(passed=False,elastic=[],point_integrations=[],restarts=[],scope='small matching Tet4/Tet10 contact BIRTH',remaining=['wet soil BIRTH', 'declared crack-tip combination', 'Richards wet JOINT', 'large CAD and GUI/GPU', 'P6 low-preload NC removal and general wet sliding full-history convergence'])
    try:
        for mode in ('transfer_force','stress_free'):
            for quadratic in (False,True):
                mesh,ground,case=fixture(mode,quadratic);name=f'dry-{mode}-'+('tet10' if quadratic else 'tet4')
                run,_,_=run_case(solver,output/name,mesh,ground,case);report['elastic'].append(dict(name=name,**audit_elastic(mesh,ground,run)))
            for water in ('hydrostatic','flow','soil','prescribed'):
                mesh,ground,case=fixture(mode,water=water);name=f'wet-{mode}-{water}'
                run,_,_=run_case(solver,output/name,mesh,ground,case);report['elastic'].append(dict(name=name,**audit_elastic(mesh,ground,run)))
        mesh,ground,case=shear_fixture()
        run,_,_=run_case(solver,output/'dry-inherited-shear',mesh,ground,case)
        report['elastic'].append(dict(name='dry-inherited-shear',**audit_elastic(mesh,ground,run)))
        inherited=run['transition_audits'][0]['prepared_state']['joints']['points']
        assert max(np.linalg.norm(np.asarray(p['traction_kpa'])-[12.,0.,-100.]) for p in inherited)<1e-7
        data=dict(material=dict(normal_stiffness_kpa_m=1e6,tangential_stiffness_kpa_m=2e5,friction_coefficient=.3),normal=[0,0,1],reference=[.02,-.01,.03],traction=[12.,5.,-100.],maximum_increment_m=1e-5,
                  increments=[[1e-5,-2e-5,-1e-5],[-1e-5,2e-5,0],[.003,.001,-.0001],[-.003,-.002,-.0001],[.003,.001,.001],[.004,-.001,.001],[.004,-.001,-.0001]])
        result=subprocess.run([str(args.point_driver.resolve())],input=json.dumps(data),text=True,capture_output=True);assert result.returncode==0,result.stderr
        native=json.loads(result.stdout);expected=point_path(data);error=0.
        for actual,reference in zip(native['states'],expected):
            for key,value in reference.items():
                if key=='initial_energy':continue
                diff=np.max(np.abs(np.asarray(actual[key])-value));error=max(error,float(diff));assert diff<1e-7,(key,diff)
            balance=reference['initial_energy']+actual['endpoint_work_kpa_m']-np.linalg.norm(np.asarray(actual['traction_kpa'])[:2])**2/(4e5)-actual['traction_kpa'][2]**2/(2e6)-actual['plastic_dissipation_kpa_m']-actual['numerical_dissipation_kpa_m'];assert abs(balance)<1e-10
        dump(output/'point-input.json',data);dump(output/'point-native.json',native);report['point_integrations'].append(dict(passed=True,states=len(expected),maximum_absolute_error=error))
        # A supported upper body can open while the declared water face remains.
        mesh,ground,case=fixture('stress_free',water='hydrostatic');case['settings']['kinematics']='explicit';case['boundary']=read(ROOT/'examples/stages-joint-column.json')['boundary']
        case['stages']=[dict(kind='joint_birth',name='contact',joints=['Seam']),dict(kind='boundary',name='slide',target_factor=-.1),dict(kind='joint_death',name='death',joints=['Seam']),dict(kind='joint_birth',name='rebirth',joints=['Seam']),dict(kind='boundary',name='open',target_factor=1.)]
        case['verification_reject_once_stages']=['contact','slide','death','rebirth']
        run,inputs,invoke=run_case(solver,output/'wet-rebirth-open',mesh,ground,case)
        report['wet_reassembly']=audit_hydrostatic_reassembly(mesh,ground,run)
        assert all(p['status']=='open' for p in run['last_committed']['joints']['points'])
        assert any(abs(f)>0 for f in run['last_committed']['joint_water']['nodal_load_kn'])
        assert run['last_committed']['joints']['lifecycle'][0]['generation']==2
        old=run['frames'][2]['state']['joints'];new=run['transition_audits'][2]['prepared_state']['joints']
        assert any(p['accumulated_plastic_slip_m']>0 and p['plastic_dissipation_kpa_m']>0 for p in old['points'])
        assert all(p['plastic_slip_m']==[0.,0.,0.] and p['plastic_dissipation_kpa_m']==0 for p in new['points'])
        clean=copy.deepcopy(case);clean['verification_reject_once_stages']=[]
        for stage in clean['stages']:
            if stage['name'] in case['verification_reject_once_stages']:stage['initial_increment']=case['settings']['initial_increment']/2
        clean_run,_,_=run_case(solver,output/'wet-rebirth-clean',mesh,ground,clean)
        assert run['last_committed']==clean_run['last_committed'] and run['frames']==clean_run['frames']
        def accepted_path(r):return [(a['stage'],a['progress_from'],a['progress_to']) for a in r['attempts'] if a['accepted']]
        assert accepted_path(run)==accepted_path(clean_run)
        assert without_clock([a for a in run['attempts'] if a['accepted']])==without_clock([a for a in clean_run['attempts'] if a['accepted']])
        report['rejection_isolation']=dict(passed=True,same_accepted_path=True,exact_state_frames_history=True,exact_accepted_attempt_history_except_clock=True)
        init=sum(a['accepted'] for a in run['geostatic']['attempts']);count=init;locations=[init]
        seen=set();stage_ends={}
        for a in run['attempts']:
            if a['accepted']:
                count+=1;stage_ends[a['stage']]=count
                if a['stage'] not in seen:locations.append(count);seen.add(a['stage'])
        locations=sorted(set(locations+list(stage_ends.values())))
        for count in locations:
            checkpoint=output/'wet-rebirth-open'/f'checkpoint-{count}.json';resumed=checkpoint.with_name(f'resumed-{count}.json')
            invoke(['checkpoint-stages',*inputs,str(count),str(checkpoint)],f'checkpoint-{count}.log')
            invoke(['resume-stages',*inputs,str(checkpoint),str(resumed)],f'resume-{count}.log')
            assert without_clock(read(resumed))==without_clock(run)
            report['restarts'].append(dict(accepted_increment=count,exact_except_elapsed_seconds=True))
        report['passed']=True
    except Exception as error:
        report['failure']=repr(error)
        raise
    finally:
        dump(output/'acceptance.json',report)
    print(json.dumps({k:v for k,v in report.items() if k not in ('elastic',)},ensure_ascii=False,indent=2))
if __name__=='__main__':main()
