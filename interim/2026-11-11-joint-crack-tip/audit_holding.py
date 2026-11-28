"""Independent zero-release and holding-force audit for full one-side removal."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
import fem_reference as ref
from joint_reference import System
from verify_high_order import prepare,internal,physical


def audit(folder):
    mesh=json.loads((folder/'mesh.json').read_text());ground=json.loads((folder/'ground.json').read_text());run=json.loads((folder/'result.json').read_text())
    assert run['complete'];case=run['case'];model=System(mesh,case['joint_system']);q,offset=prepare(mesh)
    rows=[]
    for index,stage in enumerate(case['stages']):
        if stage['kind']!='excavate':continue
        before=run['frames'][index]['state'];after=run['frames'][index+1]['state'];active=np.asarray(after['active_elements']);removed=np.asarray(before['active_elements']) & ~active
        assert removed.any() and all(not(active[p['owner_a']] and active[p['owner_b']]) for p in model.points)
        level=before.get('hydrostatic_water_level_m',ground['initialization']['water_level_m'])
        current={**ground,'initialization':{**ground['initialization'],'water_level_m':level}}
        dead,_=internal(mesh,current,before['material_states'],removed,q,offset,before.get('element_material_ids'),before.get('pore_pressure_kpa'))
        retained,_=internal(mesh,current,before['material_states'],active,q,offset,before.get('element_material_ids'),before.get('pore_pressure_kpa'))
        contact_before=model.audit(before);kept={**before,'active_elements':active.tolist()};contact_after=model.audit(kept)
        load,_=physical(mesh,current,active,stage['remaining_surface_pressure_kpa'],True,after.get('element_material_ids'),model)
        load+=after['boundary_factor']*ref.conditions_load(mesh,case['boundary'],active)
        # Every contact pair is inactive: the retained paired water field is
        # exactly zero. Newly exposed hydrostatic faces are in physical().
        if 'joint_water' in after:assert not any(after['joint_water']['nodal_load_kn'])
        assert np.max(np.abs(load.ravel()-after['applied_nodal_load_kn']))<1e-10
        holding=np.asarray(before['applied_nodal_load_kn']).reshape(-1,3)-load-dead-(contact_before-contact_after)
        used=np.unique(ref.cells(mesh)[active]);unused=np.setdiff1d(np.arange(len(mesh['nodes'])),used);holding[unused]=0.
        actual=next(a for a in run['excavation_audits'] if a['stage']==stage['name'])
        error=float(np.max(np.abs(holding.ravel()-actual['initial_holding_load_kn'])))
        assert error<1e-10,(stage['name'],error)
        basis,_=ref.staged_constraint_space(mesh,case,kept)
        residual=np.concatenate([b.T@(retained+contact_after-load-holding)[n] for n,b in basis.items()])
        relative=float(np.linalg.norm(residual)/max(np.linalg.norm(load+holding),1.));assert relative<1e-10
        injected=[a for a in run['attempts'] if a['stage']==stage['name'] and a['reason']=='verification_injected_rejection_after_convergence']
        assert injected and all(a['rejected_trial_active_set_changed'] and a['rejected_trial_load_changed'] and a['committed_state_unchanged_on_rejection'] for a in injected)
        assert before['joints']==after['joints'] and not any(after['excavation_release_load_kn'])
        rows.append(dict(stage=stage['name'],holding_force_max_abs_difference_kn=error,zero_release_reassembled_residual=relative,removed_contact_force_norm_kn=float(np.linalg.norm(contact_before-contact_after)),inactive_contact_history_exact=True,completed_release_zero=True,intended_rejection_changes_active_set_and_load=True))
    assert rows
    return dict(schema='soil-fem-joint-crack-tip-holding/1',passed=True,rows=rows)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();result=audit(a.run.resolve());a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
