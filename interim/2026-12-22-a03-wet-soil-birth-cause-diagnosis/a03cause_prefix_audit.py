"""Post-solve read-only audit: dormant pressure samples are not new birth data."""
import argparse
import hashlib
import json
from pathlib import Path
from a03cause_metrics import metric
from a03cause_driver import write


def audit(run):
    rows=[]
    for policy in ('soil','continuous','prescribed'):
        folder=run/'refusals'/('darcy-before-birth-'+policy)
        cp=json.loads((folder/'flow-complete-checkpoint.json').read_text(encoding='utf-8'))
        mesh=json.loads((folder/'mesh.json').read_text(encoding='utf-8'))
        ground=json.loads((folder/'ground.json').read_text(encoding='utf-8'))
        case=json.loads((folder/'case.json').read_text(encoding='utf-8'))
        state=cp['run']['last_committed'];points=cp['run']['integration_points']
        levels=ground['initialization'];h=levels['water_level_m'];gamma=levels['water_unit_weight_kn_m3']
        body={r['id']:r['body'] for r in mesh['regions']};structures=set(case['joint_system'].get('structure_bodies',[]))
        samples={'dormant':[],'active_soil':[]};expected={k:[] for k in samples}
        for point,p in zip(points,state['pore_pressure_kpa']):
            e=point['element_index'];kind='dormant' if not state['active_elements'][e] else 'active_soil'
            if kind=='active_soil' and body[mesh['region_ids'][e]] in structures:continue
            samples[kind].append(p);expected[kind].append(gamma*max(h-point['position_m'][2],0.))
        dormant=metric(samples['dormant'],expected['dormant']);active=metric(samples['active_soil'],expected['active_soil'])
        assert dormant['maximum_absolute']==0. and active['maximum_absolute']>0.
        rows.append(dict(policy=policy,checkpoint_sha256=hashlib.sha256((folder/'flow-complete-checkpoint.json').read_bytes()).hexdigest(),
                         dormant_GP_count=len(samples['dormant']),dormant_range_kpa=[min(samples['dormant']),max(samples['dormant'])],
                         dormant_vs_original_hydrostatic=dormant,active_soil_vs_original_hydrostatic=active,
                         active_soil_GP_count=len(samples['active_soil']),birth_executed=False,
                         interpretation='stage_pore retains start on inactive elements; array presence is not a Darcy BIRTH initialization contract'))
    return dict(diagnostic_only=True,passed=False,a03_complete=False,kind='saved checkpoint observation, no new solve',rows=rows,
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                production_source='src/stages.rs::stage_pore inactive branch / prepare_transition after-seepage guard')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();write(a.output,audit(a.run))
