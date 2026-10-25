"""Fresh CLI wet interfaces, independent flow/force, history and restart gates."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import argparse
import copy
from pathlib import Path
from verify_one_way_cases import read,write,normalized
from verify_checkpoint import invoke
from verify_acceptance_packet import verify as packet,sha
from export_stages import export
from verify_joint_water import verify_package


def verify(solver,output,cad=False):
    root=Path(__file__).resolve().parents[1];output.mkdir(parents=True,exist_ok=False);checks={};restarts=[]
    def command(args,log,expected=0):invoke([str(solver)]+[str(a) for a in args],log,expected)
    def profile(name,mesh,ground,case):
        folder=output/name;folder.mkdir()
        for n,v in (('mesh',mesh),('ground',ground),('stages',case)):write(folder/(n+'.json'),v)
        command(['solve-stages',folder/'mesh.json',folder/'ground.json',folder/'stages.json',folder/'stage-results.json'],folder/'solve.log')
        native=read(folder/'stage-results.json');export(folder/'mesh.json',folder/'stage-results.json',folder/'results');checks[name]=verify_package(folder/'mesh.json',folder/'ground.json',folder/'stage-results.json',folder/'results',solver)
        return native,folder
    def restart(run,folder,counts,mesh=None,ground=None):
        mesh=mesh or folder/'mesh.json';ground=ground or folder/'ground.json'
        for count in counts:
            cp=folder/f'checkpoint-{count}.json';result=folder/f'resumed-{count}.json'
            command(['checkpoint-stages',mesh,ground,folder/'stages.json',count,cp],folder/f'checkpoint-{count}.log');command(['resume-stages',mesh,ground,folder/'stages.json',cp,result],folder/f'resume-{count}.log');assert normalized(read(result))==normalized(run)
            restarts.append(dict(profile=folder.name,count=count,all_states_frames_attempts_flow_maps_and_water_histories_exact=True))
    def selected_counts(run,first_and_last=True):
        geo=sum(a['accepted'] for a in run['geostatic']['attempts']);first=run['case']['stages'][0]['name'];before_last=sum(a['accepted'] for a in run['attempts'] if a['stage']!=run['case']['stages'][-1]['name'])
        return [geo+1,geo+before_last+1] if first_and_last else [1,geo+1,geo+before_last+1]
    for material in ('mc','mcc'):
        for kind in ('mesh','tet10-mesh'):
            run,folder=profile(f'{material}-hydro-{kind}',read(root/f'examples/joint-water-column-{kind}.json'),read(root/f'examples/ground-joint-water-{material}.json'),read(root/'examples/stages-joint-water-column.json'));restart(run,folder,selected_counts(run))
    for policy,name in (('continuous','flow'),('soil','soil'),('prescribed','prescribed')):
        for kind in ('mesh','tet10-mesh'):
            case=read(root/f'examples/stages-joint-water-{name}.json')
            if policy=='continuous' and kind=='tet10-mesh':case['verification_reject_once_stages']=['flow']
            run,folder=profile(f'{policy}-{kind}',read(root/f'examples/joint-water-column-{kind}.json'),read(root/'examples/ground-joint-water-mc.json'),case);restart(run,folder,selected_counts(run,False))
            if policy=='continuous' and kind=='tet10-mesh':
                assert any(not a['accepted'] and a['rejected_trial_water_changed'] and a['committed_state_unchanged_on_rejection'] for a in run['attempts']);clean=copy.deepcopy(case);clean['verification_reject_once_stages']=[];clean['stages'][0]['initial_increment']=.125
                direct,_=profile('continuous-clean',read(folder/'mesh.json'),read(folder/'ground.json'),clean);assert normalized(direct['frames'])==normalized(run['frames']) and direct['last_committed']==run['last_committed']
    mesh=read(root/'examples/joint-water-column-mesh.json');ground=read(root/'examples/ground-joint-water-mc.json');case=read(root/'examples/stages-joint-water-cycle.json');run,folder=profile('cycle',mesh,ground,case);restart(run,folder,selected_counts(run));assert any(not a['accepted'] and a['rejected_trial_joint_slip_changed'] for a in run['attempts'])
    opened=next(f for f in run['frames'] if f['name']=='open');assert all(p['opening_m']>1e-6 and p['normal_pressure_kpa']==0. and p['hydraulic']['physical_pore_pressure_kpa']==2.5 for p in opened['joint_results']['points'])
    mesh=read(root/'examples/joint-water-column-tet10-mesh.json');case=read(root/'examples/stages-joint-water-flow.json');case['stages']=case['stages'][:1]+[dict(kind='excavate',name='cut',bodies=['UpperClay'],remaining_surface_pressure_kpa=0.)];case['verification_reject_once_stages']=['cut'];run,folder=profile('death',mesh,ground,case);restart(run,folder,selected_counts(run));assert run['frames'][1]['state']['joints']==run['frames'][2]['state']['joints'] and all(v==0. for v in run['frames'][2]['state']['joint_water']['nodal_load_kn'])
    run,folder=profile('mcc-flow',mesh,read(root/'examples/ground-joint-water-mcc.json'),read(root/'examples/stages-joint-water-flow.json'));restart(run,folder,selected_counts(run))
    # Invalid declarations create no accepted result. A non-admissible saturated
    # flow solution preserves the actual preceding mechanical/water state.
    bad=read(output/'continuous-tet10-mesh/stages.json');bad['joint_system']['joints'][0].pop('water');write(output/'invalid-policy.json',bad)
    command(['solve-stages',output/'continuous-tet10-mesh/mesh.json',output/'continuous-tet10-mesh/ground.json',output/'invalid-policy.json',output/'invalid-result.json'],output/'invalid.log',1);assert not(output/'invalid-result.json').exists()
    bad=read(output/'continuous-tet10-mesh/stages.json');bad['verification_reject_once_stages']=[];bad['stages'][0]['flow']['heads'][1]['head_m']=-1.;write(output/'negative-target.json',bad)
    command(['solve-stages',output/'continuous-tet10-mesh/mesh.json',output/'continuous-tet10-mesh/ground.json',output/'negative-target.json',output/'negative-result.json'],output/'negative.log',1);failed=read(output/'negative-result.json');assert not failed['complete'] and failed['last_committed']==failed['frames'][0]['state'] and failed['flow_audits'] and failed['last_committed'].get('pore_pressure_kpa') is None
    # Wet fields are absent from the old dry inputs/results. Compare against
    # the sealed dry contact specimen rather than assuming compatibility.
    old=root/'runs/joint-test-20261024-045026-efdb1ae7';accepted=packet(old);legacy=[]
    for material in ('mc','mcc'):
        source=old/material;folder=output/f'dry-{material}';folder.mkdir();baseline=read(source/'stage-results.json')
        command(['solve-stages',source/'ground/mesh.json',source/'ground/ground.json',source/'stages.json',folder/'stage-results.json'],folder/'solve.log');current=read(folder/'stage-results.json');assert current['last_committed']==baseline['last_committed'] and current['frames']==baseline['frames'];legacy.append(dict(material=material,exact_frames_and_committed_state=True,baseline_result_sha256=sha(source/'stage-results.json'),new_result_sha256=sha(folder/'stage-results.json')))
    if cad:
        import subprocess
        folder=output/'cad-local';args=['pwsh','-NoProfile','-File',str(root/'scripts/Invoke-Joints.ps1'),'-Design',str(root/'examples/joint-water-local-cylinder.json'),'-GroundFile',str(root/'examples/ground-joint-water-local-mc.json'),'-StageFile',str(root/'examples/stages-joint-water-local.json'),'-SplitFile',str(root/'examples/joint-column-split.json'),'-ElementOrder','2','-MeshSizeM','4','-SolverPath',str(solver),'-OutputDirectory',str(folder)]
        invoke(args,output/'cad.log');run=read(folder/'analysis/stage-results.json');checks['cad-local']=read(folder/'analysis/stage-verification.json');assert run['provenance']['solver_binary_sha256']==sha(solver);restart(run,folder/'analysis',selected_counts(run),folder/'analysis/ground/mesh.json',folder/'analysis/ground/ground.json')
        points=run['joint_integration_points'];oldfield=run['frames'][-2]['state']['joint_water'];field=run['last_committed']['joint_water'];removed=sum(a>0. and b==0. for a,b in zip(oldfield['physical_pressure_kpa'],field['physical_pressure_kpa']));retained=sum(p>0. for p in field['physical_pressure_kpa']);assert removed>0 and retained>0
        mesh=read(folder/'split/mesh.json');write(folder/'contact-domain-verification.json',dict(passed=True,nodes=len(mesh['nodes']),tetrahedra=len(mesh['tetrahedra']),integration_points=len(points),removed_wet_contact_points=removed,retained_wet_contact_points=retained,fresh_solidcraft_capture=True,source_design_sha256=sha(root/'examples/joint-water-local-cylinder.json')))
    result=dict(schema='soil-fem-wet-joint-cli-verification/1',passed=True,solver_sha256=sha(solver),cad_exercised=bool(cad),profiles=checks,restarts=restarts,dry_joint_legacy_cases=legacy,dry_joint_parent_integrity=accepted,flow_cutback_clean_path_exact=True,wet_contact_open_pressure_and_rejected_slip_history_verified=True,negative_saturated_flow_target_preserves_committed_state=True,missing_pressure_policy_rejected_before_output=True)
    write(output/'verification.json',result);print('Wet JOINT CLI:',len(checks),'profiles,',len(restarts),'fresh-process restarts, 2 exact dry contact cases')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',required=True,type=Path);p.add_argument('--output-dir',required=True,type=Path);p.add_argument('--cad',action='store_true');a=p.parse_args();verify(a.solver.resolve(),a.output_dir.resolve(),a.cad)
