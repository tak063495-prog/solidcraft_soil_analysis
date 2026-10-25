"""Fresh-process CLI, export, discrete audit, failure and restart acceptance."""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import shutil
from export_richards import export
from verify_checkpoint import invoke
from verify_richards import audit


def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,data):
    with p.open('x',encoding='utf-8') as stream:json.dump(data,stream,indent=2,ensure_ascii=False)
def normalized(r):
    r=copy.deepcopy(r)
    for a in r['attempts']:
        for i in a['iterations']:
            if i['linear']:i['linear']['elapsed_seconds']=0.
    return r


def fixture(solver,mesh,case,output,export_result=True):
    output.mkdir(parents=True,exist_ok=False);write(output/'mesh.json',mesh);write(output/'case.json',case)
    inputs=[str(output/'mesh.json'),str(output/'case.json')]
    invoke([str(solver),'solve-richards',*inputs,str(output/'result.json')],output/'solve.log')
    result=read(output/'result.json');assert result['solver_binary_sha256']==sha(solver)
    checked=audit(mesh,result);write(output/'verification.json',checked)
    if export_result:export(output/'mesh.json',output/'result.json',solver,output/'results')
    return result


def verify(solver,output,cad_mesh=None):
    root=Path(__file__).resolve().parents[1];output.mkdir(parents=True,exist_ok=False);profiles={}
    for name,order,case_name in [('wet','tet4','richards-steady-wet'),('dry','tet4','richards-steady-dry'),('rain-tet4','tet4','richards-rain-pause'),('rain-tet10','tet10','richards-rain-pause')]:
        profiles[name]=fixture(solver,read(root/f'examples/seepage-column-{order}-mesh.json'),read(root/f'examples/{case_name}.json'),output/name)
    for name in ('wet','dry'):
        f=profiles[name]['frames'][-1];q=sum(a['outward_flux_m3_s'] for a in f['balance']['atmospheric_faces'])
        assert abs(q-(1e-6 if name=='wet' else 0.))<1e-13
        assert all(a['wet']==(name=='wet') for a in f['balance']['atmospheric_faces'])
        assert all(a['seeping']==(name=='wet') for a in f['balance']['atmospheric_faces'])
    unicode_mesh=read(root/'examples/seepage-column-tet4-mesh.json');unicode_case=read(root/'examples/richards-steady-wet.json')
    unicode_name='\u5730\u8868\u9762';unicode_mesh['face_sets'][unicode_name]=unicode_mesh['face_sets']['top']
    unicode_case['analysis']['boundaries'][1]['face_set']=unicode_name
    fixture(solver,unicode_mesh,unicode_case,output/'utf8')
    restarts=[]
    for name in ('rain-tet4','rain-tet10'):
        folder=output/name;result=profiles[name];frames=result['frames'];counts=[1]+[max(i for i,f in enumerate(frames) if f['name']==period) for period in ('rain','pause')]+[len(frames)-2]
        for number,count in enumerate(counts):
            cp=folder/f'checkpoint-{number}.json';resumed=folder/f'resumed-{number}.json';inputs=[str(folder/'mesh.json'),str(folder/'case.json')]
            invoke([str(solver),'checkpoint-richards',*inputs,str(count),str(cp)],folder/f'checkpoint-{number}.log')
            saved=read(cp);assert saved['run']['runtime'] and not saved['run']['complete']
            invoke([str(solver),'resume-richards',*inputs,str(cp),str(resumed)],folder/f'resumed-{number}.log')
            assert normalized(read(resumed))==normalized(result)
            restarts.append(dict(profile=name,count=count,time_s=saved['run']['last_committed']['time_s'],exact_all_states_frames_and_attempts=True))
        bad=read(folder/'checkpoint-0.json');bad['run']['last_committed']['cell_storage_offset_m3'][0]=0.;write(folder/'bad-checkpoint.json',bad)
        invoke([str(solver),'resume-richards',str(folder/'mesh.json'),str(folder/'case.json'),str(folder/'bad-checkpoint.json'),str(folder/'bad-resumed.json')],folder/'bad-checkpoint.log',1)
        assert not (folder/'bad-resumed.json').exists()
        changed=read(folder/'case.json');changed['materials'][0]['retention']['alpha_per_m']*=2.;write(folder/'changed-case.json',changed)
        invoke([str(solver),'resume-richards',str(folder/'mesh.json'),str(folder/'changed-case.json'),str(folder/'checkpoint-0.json'),str(folder/'changed-resumed.json')],folder/'changed-case.log',1)
        assert not (folder/'changed-resumed.json').exists()
    # Genuine non-convergence, then exact replay with the accepted physical dt.
    case=read(root/'examples/richards-rain-pause.json');case['analysis']['periods']=case['analysis']['periods'][:1];case['analysis']['periods'][0]['end_time_s']=10000.
    case['settings']=dict(initial_time_step_s=10000.,maximum_time_step_s=10000.,maximum_iterations=2,store_integration_points=False)
    mesh=read(root/'examples/seepage-column-tet4-mesh.json');cut=fixture(solver,mesh,case,output/'cutback',False)
    assert any(not a['accepted'] for a in cut['attempts'])
    clean=copy.deepcopy(case);bc=clean['analysis']['periods'][0]['boundaries'];clean['analysis']['periods']=[dict(name=f'clean-{i}',end_time_s=f['state']['time_s'],boundaries=bc) for i,f in enumerate(cut['frames'][1:])];clean['settings']['maximum_iterations']=100
    clean_result=fixture(solver,mesh,clean,output/'clean',False)
    assert [f['state'] for f in cut['frames']]==[f['state'] for f in clean_result['frames']]
    # Reproduce the pre-increment/pre-tangent failure with the same small flow budget.
    tight=read(root/'examples/richards-rain-pause.json');tight['analysis']['periods']=tight['analysis']['periods'][:1];tight['analysis']['periods'][0]['end_time_s']=100.
    tight['settings']=dict(initial_time_step_s=50.,maximum_time_step_s=50.,absolute_flow_tolerance_m3_s=1e-22,store_integration_points=False)
    fixture(solver,read(root/'tests/data/richards-tight-flow-mesh.json'),tight,output/'tight-flow',False)
    failed=copy.deepcopy(case);failed['settings'].update(maximum_iterations=1,minimum_time_step_s=10000.)
    write(output/'failed-case.json',failed)
    invoke([str(solver),'solve-richards',str(output/'wet/mesh.json'),str(output/'failed-case.json'),str(output/'failed-result.json')],output/'failed.log',1)
    report=read(output/'failed-result.json');assert not report['complete'] and report['failure_reason'] and len(report['frames'])==1
    assert report['last_committed']==report['frames'][0]['state']
    invalid=read(root/'examples/richards-steady-wet.json');invalid['analysis']['boundaries'][0]['face_set']='missing';write(output/'invalid-case.json',invalid)
    invoke([str(solver),'solve-richards',str(output/'wet/mesh.json'),str(output/'invalid-case.json'),str(output/'invalid-result.json')],output/'invalid.log',1)
    assert not (output/'invalid-result.json').exists()
    zero=read(root/'examples/richards-steady-wet.json')
    for material in zero['materials']:material['retention']['saturated_specific_storage_per_m']=0.
    zero['analysis']=dict(kind='transient',periods=[dict(name='closed',end_time_s=1.,boundaries=[])])
    write(output/'zero-storage-case.json',zero)
    invoke([str(solver),'solve-richards',str(output/'wet/mesh.json'),str(output/'zero-storage-case.json'),str(output/'zero-storage-result.json')],output/'zero-storage.log',1)
    zero_result=read(output/'zero-storage-result.json');assert not zero_result['complete'] and 'zero-storage' in zero_result['failure_reason']
    assert zero_result['last_committed']==zero_result['frames'][0]['state']
    existing=sha(output/'wet/result.json')
    invoke([str(solver),'solve-richards',str(output/'wet/mesh.json'),str(output/'wet/case.json'),str(output/'wet/result.json')],output/'existing-output.log',1)
    assert sha(output/'wet/result.json')==existing
    if cad_mesh:
        mesh=read(cad_mesh);case=read(root/'examples/richards-steady-wet.json');case['initial_head']['head_m']=8.5
        case['analysis']['boundaries']=[dict(face_set='bottom',condition=dict(kind='head',value=dict(head_m=8.5)))]
        cad=fixture(solver,mesh,case,output/'cad-curved-hydrostatic')
        assert all(p['darcy_flux_m_s']==[0.,0.,0.] for p in cad['frames'][-1]['points'])
    report=dict(schema='soil-fem-richards-cli-verification/1',passed=True,solver_sha256=sha(solver),restarts=restarts,
        exact_cutback_accepted_history=True,genuine_exhaustion_preserves_initial_state=True,invalid_input_no_result=True,existing_output_preserved=True,
        cad_curved_hydrostatic_exercised=bool(cad_mesh),scope='Hydraulic core only. No suction mechanics, pond storage, routing dynamics or arbitrary 3D seepage/free-surface validation is inferred.')
    report['zero_storage_missing_pressure_reference_rejected']=True
    write(output/'verification.json',report);print('Richards CLI acceptance passed:',len(restarts),'fresh-process restarts')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--cad-mesh',type=Path);a=p.parse_args()
    verify(a.solver.resolve(),a.output_dir.resolve(),a.cad_mesh.resolve() if a.cad_mesh else None)
