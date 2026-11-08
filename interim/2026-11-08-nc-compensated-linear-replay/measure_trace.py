"""Audit this bounded trial from stored stresses; not an independent plastic solve."""
import argparse
from collections import Counter
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
spec=importlib.util.spec_from_file_location('old_nc_diagnostic_measure',ROOT/'interim/2026-11-08-nc-midpoint-stop-diagnosis/measure.py')
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)


def measure(evidence, trial, output):
    outcome=json.loads((trial/'diagnostic.json').read_text())
    assert outcome['saved_checkpoint_unchanged'] and not outcome['outcome']['converged'], 'bounded failure changed; retain actual result without certifying it'
    folder=evidence/'prior/nc-resumed'
    mesh,ground,case=[json.loads((folder/'input'/f'{n}.json').read_text()) for n in ('mesh','ground','stages')]
    assert (mesh,ground,case)==old.inputs('midpoint')
    run=json.loads((folder/'job/recovery-checkpoint-0.json').read_text())['run']
    state=run['last_committed'];origin=run['frames'][2]['state'];c=run['continuation']
    q,offset=old.prepare(mesh,run['integration_points']);active=np.array(state['active_elements'])
    indices=np.concatenate([np.arange(offset[i],offset[i+1]) for i in np.flatnonzero(active)])
    params={m['id']:m['mcc'] for m in ground['materials']}
    point_params=[dict(params[mesh['material_ids'][ei]],integration='midpoint') for ei in range(len(q)) for _ in q[ei]['weight']]
    events=[json.loads(l) for l in (trial/'trace.jsonl').read_text().splitlines()]
    counts=Counter(e['kind'] for e in events);rows=[];material=[];spectra=[]
    for e in events:
        if e['kind'] in ('iterate','line_search'):
            assert 'error' not in e
            assert e['joint_state']==origin['joints']
            assert all(not p['active'] for p in e['joint_report']['points'])
            for ei in np.flatnonzero(~active):
                assert e['states'][offset[ei]:offset[ei+1]]==origin['material_states'][offset[ei]:offset[ei+1]]
        if e['kind']=='iterate':
            t=old.transform(e['dof_rows']);bulk,_=old.internal(mesh,ground,e['states'],active,q,offset)
            joint=np.array(e['joint_internal']);assert np.count_nonzero(joint)==0
            full=bulk.ravel()+joint;target=np.array(e['target']);end=np.array(c['end_load']);holding=(1-e['progress'])*np.array(c['holding_load'])
            assert np.array_equal(target,end+holding)
            delta=np.linalg.norm(full-np.array(e['internal']));assert delta<1e-10
            reduced=t.T@(target-full);assert np.linalg.norm(reduced-np.array(e['reduced_residual']))<1e-10
            points=old.material_rows(e['states'],state['material_states'],indices,point_params)
            material.append(dict(iteration=e['iteration'],points=points))
            rows.append(dict(iteration=e['iteration'],relative_force_residual=e['relative_force_residual'],
                independent_stored_stress_force_residual=float(np.linalg.norm(reduced)/max(np.linalg.norm(target),1.)),
                independent_internal_difference_kn=float(delta),
                partitions_min=int(min(e['partitions'][i] for i in indices)),partitions_max=int(max(e['partitions'][i] for i in indices)),
                yield_max=max(p['normalized_yield'] for p in points),softening_points=sum(p['delta_pc_kpa']<0 for p in points)))
        elif e['kind']=='tangent':
            a=np.array(e['matrix_columns']).T;s=np.linalg.svd(a,compute_uv=False)
            spectra.append(dict(iteration=e['iteration'],singular_values_kn_m=s.tolist(),condition_number_2=float(s[0]/s[-1])))
    linear=[e['event']['linear_solve'] for e in events if e['kind']=='progress' and e['event']['kind']=='linear_solved']
    assert len(linear)==60 and all(p['true_relative_residual']<=1e-10 for p in linear)
    result=dict(passed=False,history_converged=False,production_adopted=False,
                iterates=counts['iterate'],line_search_candidates=counts['line_search'],linear_solves=len(linear),
                accepted_new_increments=0,source_progress=state['stage_progress'],
                target_progress=state['stage_progress']+c['increment'],minimum_increment=case['settings']['minimum_increment'],
                source_increment_is_minimum=c['increment']==case['settings']['minimum_increment'],
                repartition_iterations=[e['iteration'] for e in events if e['kind']=='repartition'],
                maximum_original_linear_residual=max(p['true_relative_residual'] for p in linear),
                total_refinement_corrections=sum(p['iterations']-1 for p in linear),
                final_relative_force_residual=rows[-1]['relative_force_residual'],
                minimum_relative_force_residual=min(r['relative_force_residual'] for r in rows),
                maximum_contact_active_points=0,contact_history_unchanged=True,removed_material_history_unchanged=True,
                origin_and_checkpoint_lambda_zero='unchanged old sealed state; no new physical acceptance',
                independent_plastic_3d_solution=False,full_path_history=False,full_removal_reached=False,
                force_rows=rows,spectra=spectra)
    output.mkdir(parents=True,exist_ok=False)
    (output/'trial-measurements.json').write_text(json.dumps(result,indent=2)+'\n')
    (output/'material-history.json').write_text(json.dumps(material,indent=2)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--evidence',type=Path,required=True);p.add_argument('--trial',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=measure(a.evidence,a.trial,a.output)
    print(json.dumps({k:v for k,v in r.items() if k not in ('force_rows','spectra')},indent=2))
