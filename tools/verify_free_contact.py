"""Free 3D contact initialization, independent solves, rollback and restart."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys
import numpy as np
from make_free_contact import make
from run_joint_fixture import run as fixture
from verify_checkpoint import invoke
from verify_joint_coupling import reference


def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,v):
    with p.open('x',encoding='utf-8') as stream:json.dump(v,stream,indent=2)


def accepted(run):
    return ([(a['progress_to']) for a in run['geostatic']['attempts'] if a['accepted']],
            [(a['stage'],a['progress_to']) for a in run['attempts'] if a['accepted']])


def clean_replay(solver,folder):
    native=read(folder/'stage-results.json');case=read(folder/'stages.json');ground=read(folder/'ground/ground.json')
    output=folder/'clean';output.mkdir(exist_ok=False)
    ground['initialization']['geostatic']['verification_reject_once']=False
    case['verification_reject_once_stages']=[]
    case['settings']['initial_increment']=next(a['progress_to'] for a in native['geostatic']['attempts'] if a['accepted'])
    for stage in case['stages']:stage['initial_increment']=next(a['progress_to'] for a in native['attempts'] if a['accepted'] and a['stage']==stage['name'])
    write(output/'ground.json',ground);write(output/'stages.json',case)
    invoke([str(solver),'solve-stages',str(folder/'ground/mesh.json'),str(output/'ground.json'),str(output/'stages.json'),str(output/'result.json')],output/'solve.log')
    clean=read(output/'result.json');assert clean['complete']
    assert accepted(clean)==accepted(native),'clean accepted path differs'
    assert clean['last_committed']==native['last_committed'] and clean['frames']==native['frames'],'rejected history leaked'
    return {'exact_accepted_path_states_and_frames':True}


def independent(folder,steps=100,output=None):
    mesh=read(folder/'ground/mesh.json');ground=read(folder/'ground/ground.json');native=read(folder/'stage-results.json')
    case=json.loads(json.dumps(native['case']));case['stages']=case['stages'][:-1]
    assert all(s['kind']=='boundary' for s in case['stages'])
    assert all(p['equivalent_plastic_strain']==0 for f in native['frames'] for p in f['state']['material_states'])
    output=output or folder/'independent';output.mkdir(exist_ok=False);results=[];residuals=[];diagnostics=[]
    for count in (steps,2*steps):
        diag={};frames,res=reference(mesh,ground,case,native,count,True,diag)
        write(output/f'reference-{count}.json',{'frames':[{k:v.tolist() for k,v in f.items()} for f in frames],
            'maximum_relative_force_residual':res,'diagnostics':diag})
        results.append(frames);residuals.append(res);diagnostics.append(diag)
    rows=[]
    for frame,a,b in zip(native['frames'][:-1],*results):
        u=np.asarray(frame['state']['displacement_m']).reshape(-1);scale=max(np.linalg.norm(b['u']),1e-12)
        error=float(np.linalg.norm(u-b['u'])/scale);refinement=float(np.linalg.norm(a['u']-b['u'])/scale)
        q=np.asarray([p['traction_kpa'] for p in frame['state']['joints']['points']]);traction=float(np.linalg.norm(q-b['traction'])/max(np.linalg.norm(b['traction']),1.))
        assert error<.005 and refinement<1e-4 and traction<.01,(folder.name,frame['name'],error,refinement,traction)
        rows.append({'frame':frame['name'],'relative_all_node_3d_displacement_error':error,'reference_refinement_error':refinement,'relative_traction_error':traction})
    assert max(residuals)<1e-8
    report={'schema':'soil-fem-independent-free-contact/1','frames':rows,'reference_subincrements':[steps,2*steps],
        'reference':'Input-only GL125 elastic MC bulk, independent friction projection, bilateral elastic starting iterate followed by dense SVD least-squares Newton with full physical force line search and independent rigid-support gate. Native accepted load fractions only; no native displacement, stress, tangent or history generates the reference.',
        'maximum_independent_force_residual':max(residuals),'diagnostics':diagnostics}
    write(output/'verification.json',report);return report


def summarize(folder):
    data=read(folder/'stage-results.json');assert data['complete']
    assert all(c['node_set']=='bottom' for c in data['case']['boundary']['constraints'])
    assert all(f['joint_support']['unrestrained_modes']==0 for f in data['frames'])
    attempts=data['geostatic']['attempts']+data['attempts'];rejected=[a for a in attempts if not a['accepted']]
    assert rejected and all(a['committed_state_unchanged_on_rejection'] for a in rejected)
    predictors=[p for a in attempts for p in a.get('contact_predictors',[])]
    assert predictors and any(p['activation_predictor'] for p in predictors)
    for p in predictors:
        assert p['temporary_linear_solve']['true_relative_residual']<=data['case']['settings']['relative_linear_tolerance']
        assert np.isfinite(p['original_relative_linear_residual'])
        if p['activation_predictor']:assert p['original_relative_linear_residual']>data['case']['settings']['relative_linear_tolerance']
        if p['accepted_alpha'] is not None:assert p['full_relative_force_after']<p['full_relative_force_before']
    for a in attempts:
        assert all(s['true_relative_residual']<=data['case']['settings']['relative_linear_tolerance'] for s in a['linear_solves'])
    assert data['frames'][-2]['state']['joints']==data['frames'][-1]['state']['joints']
    assert all(not p['active'] for p in data['frames'][-1]['joint_results']['points'])
    assert data['excavation_audits'][0]['zero_release_relative_force_residual']<1e-8
    return {'frames':len(data['frames']),'nodes':len(read(folder/'ground/mesh.json')['nodes']),
        'geostatic_accepted':sum(a['accepted'] for a in data['geostatic']['attempts']),
        'geostatic_natural_rejections':sum(not a['accepted'] and a['reason']!='verification_injected_rejection_after_convergence' for a in data['geostatic']['attempts']),
        'activation_predictors':sum(p['activation_predictor'] for p in predictors),
        'maximum_original_predictor_linear_residual':max(p['original_relative_linear_residual'] for p in predictors),
        'plastic_bulk_rejected_trial_changed':any(a['rejected_trial_plastic_changed'] for a in rejected)}


def dense_compare(solver,folder,output):
    source=read(folder/'stage-results.json');case=read(folder/'stages.json')
    case['settings']['linear_solver']={'method':'dense_lu'}
    case['settings']['maximum_free_dofs']=min(2000,case['settings']['maximum_free_dofs'])
    case_path=output.parent/(output.name+'-stages.json');write(case_path,case)
    fixture(solver,folder/'ground/ground.json',case_path,output,folder/'ground/mesh.json')
    dense=read(output/'stage-results.json');assert accepted(source)==accepted(dense)
    comparisons=[]
    for a,b in zip(source['frames'],dense['frames']):
        fields=[('u',a['state']['displacement_m'],b['state']['displacement_m']),
            ('stress',[p['effective_stress_compression_positive_kpa'] for p in a['state']['material_states']],
                      [p['effective_stress_compression_positive_kpa'] for p in b['state']['material_states']]),
            ('joint_traction',[p['traction_kpa'] for p in a['state']['joints']['points']],
                              [p['traction_kpa'] for p in b['state']['joints']['points']])]
        row={'frame':a['name']}
        for name,x,y in fields:
            x=np.asarray(x);y=np.asarray(y);error=float(np.linalg.norm(x-y)/max(np.linalg.norm(y),1e-12))
            assert error<1e-7,(name,error);row[name+'_relative_l2_error']=error
        comparisons.append(row)
    report={'same_accepted_path':True,'frames':comparisons}
    write(output/'sparse-dense-verification.json',report);return report


def verify(solver,output,cad=None):
    root=Path(__file__).resolve().parents[1];output.mkdir(parents=True,exist_ok=False)
    inputs=output/'inputs';inputs.mkdir();cases={};comparisons={};restarts=0
    specifications=[('horizontal',dict(force=(1.,.5,0.))),('inclined',dict(angle=10.,force=(0.,.5,0.))),
        ('near-cap',dict(angle=10.,force=(1.,.5,0.))),('stiff',dict(angle=5.,friction=.1,force=(0.,.2,0.),tangential_stiffness=1e7)),
        ('mcc-oc',dict(mcc=True,force=(1.,.5,0.))),('mcc-nc',dict(mcc=True,force=(1.,.5,0.),pc=100.,preload_kpa=100.))]
    folders=[]
    for name,options in specifications:
        prepared=inputs/name;prepared.mkdir();mesh,ground,case=make(**options)
        for label,value in (('mesh',mesh),('ground',ground),('stages',case)):write(prepared/(label+'.json'),value)
        folder=output/name;fixture(solver,prepared/'ground.json',prepared/'stages.json',folder,prepared/'mesh.json');folders.append((name,folder))
    if cad is not None:folders.append(('cad',cad))
    for name,folder in folders:
        cases[name]=summarize(folder);cases[name].update(clean_replay(solver,folder))
        if name=='mcc-nc':assert cases[name]['plastic_bulk_rejected_trial_changed'],'NC case must change bulk plastic history before rejection'
        if not name.startswith('mcc'):comparisons[name]=independent(folder)
        # All four phases for the new horizontal, inclined, near-cap and NC
        # problems; OC and stiff variants retain independent/clean checks.
        if name in ('horizontal','inclined','near-cap','mcc-nc','cad'):
            cp=output/(name+'-checkpoint')
            invoke([sys.executable,str(root/'tools/verify_construction_checkpoint.py'),'--solver',str(solver),'--run',str(folder),'--output-dir',str(cp)],output/(name+'-checkpoint.log'))
            record=read(cp/'verification.json');assert len(record['restarts'])==4
            restarts+=len(record['restarts'])
        print('Free contact verified:',name,flush=True)
    sparse_dense=dense_compare(solver,output/'horizontal',output/'dense')
    failures={}
    for name,mu,force in [('frictionless',0.,(0.,0.,0.)),('overcapacity',.3,(5.,0.,0.)),('critical',.3,(3.,0.,0.)),('opening',.3,(0.,0.,15.))]:
        folder=output/('failure-'+name);folder.mkdir();mesh,ground,case=make(friction=mu,force=force)
        for label,value in (('mesh',mesh),('ground',ground),('stages',case)):write(folder/(label+'.json'),value)
        invoke([str(solver),'solve-stages',*[str(folder/(n+'.json')) for n in ('mesh','ground','stages')],str(folder/'result.json')],folder/'solve.log',expected=1)
        data=read(folder/'result.json');assert not data['complete'] and not data['frames']
        assert data['last_committed']==data['geostatic']['seed']
        assert all(not a['accepted'] and a['committed_state_unchanged_on_rejection'] for a in data['geostatic']['attempts'])
        normal=10-force[2];tangent=float(np.linalg.norm(force[:2]));margin=mu*normal-tangent
        assert name=='frictionless' or normal<=0 or margin<=1e-12
        failures[name]={'incomplete':True,'seed_history_unchanged':True,'upper_normal_force_kn':normal,'upper_tangential_force_kn':tangent,'coulomb_capacity_margin_kn':margin,
            'certificate':'Global upper-body force balance: nonpositive tensionless normal capacity or T >= mu N gives no supported unique equilibrium; frictionless unloaded tangent gives three rigid modes. This does not interpret numerical nonconvergence alone as collapse.'}
    report={'schema':'soil-fem-free-contact-verification/1','passed':True,'solver_sha256':sha(solver),'cases':cases,'comparisons':comparisons,
        'separate_process_exact_restarts':restarts,'failures':failures,'sparse_dense':sparse_dense,
        'reference_scope':'Independent global 3D displacement covers elastic MC bulk with nonlinear JOINT. MCC plastic bulk has force, history, clean replay and fresh-process restart checks, not an independent plastic 3D displacement reference.',
        'nc_loading_scope':'NC seed p=pc=100 kPa, then 100 kPa preload before unloading and excavation. The 5 kPa-preload NC full-excavation prototype remains unresolved and is not certified by this suite.'}
    write(output/'verification.json',report)
    print(json.dumps({k:v for k,v in report.items() if k not in ('comparisons',)},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--cad-run',type=Path)
    a=p.parse_args();verify(a.solver.resolve(),a.output_dir.resolve(),a.cad_run.resolve() if a.cad_run else None)
