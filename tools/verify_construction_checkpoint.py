"""Fresh-process save/resume through placement, material changes and replacement."""
import argparse
import hashlib
import json
from pathlib import Path
from verify_checkpoint import invoke
from verify_geostatic_checkpoint import stable_report


def verify(solver,run,output):
    output.mkdir(parents=True,exist_ok=False);baseline=json.loads((run/'stage-results.json').read_text())
    assert baseline['complete'] and baseline['schema'] in ('soil-fem-stage-result/5','soil-fem-stage-result/6','soil-fem-stage-result/7','soil-fem-stage-result/8','soil-fem-stage-result/9','soil-fem-stage-result/10','soil-fem-stage-result/11','soil-fem-stage-result/12','soil-fem-stage-result/14','soil-fem-mixed-mcc-stage-result/2','soil-fem-projected-mcc-stage-result/2')
    digest=hashlib.sha256(solver.read_bytes()).hexdigest();assert baseline['provenance']['solver_binary_sha256']==digest
    inputs=[str(run/'ground/mesh.json'),str(run/'ground/ground.json'),str(run/'stages.json')]
    geo=sum(a['accepted'] for a in baseline['geostatic']['attempts']) if baseline.get('geostatic') else 0
    count=geo;stops=[('initialization',1)] if baseline.get('geostatic') else []
    for stage in baseline['case']['stages']:
        accepted=sum(a['accepted'] and a['stage']==stage['name'] for a in baseline['attempts'])
        assert accepted>0 and (stage['kind']=='reset_display' or accepted>1);stops.append((stage['name'],count+1));count+=accepted
    records=[]
    for name,stop in stops:
        cp=output/f'{name}-checkpoint.json';rp=output/f'{name}-resumed.json'
        invoke([str(solver),'checkpoint-stages',*inputs,str(stop),str(cp)],output/f'{name}-checkpoint.log')
        saved=json.loads(cp.read_text());assert not saved['run']['complete']
        state=saved['run']['last_committed'];release=state.get('transition_release_load_kn',[])
        invoke([str(solver),'resume-stages',*inputs,str(cp),str(rp)],output/f'{name}-resumed.log')
        resumed=json.loads(rp.read_text());assert resumed['complete']
        for key in ('last_committed','frames','excavation_audits'):
            assert resumed[key]==baseline[key],(name,key)
        assert resumed.get('reference_displacement_m')==baseline.get('reference_displacement_m')
        for key in ('transition_audits','joint_integration_points'):
            assert resumed.get(key,[])==baseline.get(key,[]),(name,key)
        def stable_flow(audits):
            data=json.loads(json.dumps(audits))
            for audit in data:
                if audit['solution']['linear'] is not None:audit['solution']['linear']['elapsed_seconds']=0.
                for solve in audit['solution'].get('refinement_solves',[]):solve['elapsed_seconds']=0.
                for balance in audit['solution']['balances']:
                    for solve in balance['linear_solves']:solve['elapsed_seconds']=0.
            return data
        assert stable_flow(resumed.get('flow_audits',[]))==stable_flow(baseline.get('flow_audits',[]))
        if baseline.get('geostatic'):assert stable_report(resumed['geostatic'])==stable_report(baseline['geostatic'])
        else:assert resumed.get('geostatic') is None
        assert stable_report({'attempts':resumed['attempts']})==stable_report({'attempts':baseline['attempts']})
        records.append({'stage':name,'stop_after_accepted_increments':stop,'saved_progress':state['stage_progress'],
            'saved_hydrostatic_water_level_m':state.get('hydrostatic_water_level_m'),
            'saved_nonzero_transition_holding':any(abs(x)>1e-9 for x in release),
            'saved_nonzero_excavation_holding':any(abs(x)>1e-9 for x in state.get('excavation_release_load_kn',[])),
            'saved_anchor_force_control':any(a['force_control_kn'] is not None for a in state.get('anchors',{}).get('anchors',[])),
            'exact_states_frames_material_ids_and_audits':True})
    transitions={s['name'] for s in baseline['case']['stages'] if s['kind'] in ('birth','material_set','reset_state','install_beams','remove_beams')}
    if transitions:assert any(r['stage'] in transitions and r['saved_nonzero_transition_holding'] for r in records)
    jacking={s['name'] for s in baseline['case']['stages'] if s['kind']=='jack_anchors'}
    if jacking:assert any(r['stage'] in jacking and r['saved_anchor_force_control'] for r in records)
    report={'schema':'soil-fem-construction-checkpoint-verification/1','solver_sha256':digest,'restarts':records,
        'sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in output.glob('*.json')}}
    with (output/'verification.json').open('x') as f:json.dump(report,f,indent=2)
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--run',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();verify(a.solver.resolve(),a.run.resolve(),a.output_dir.resolve())
