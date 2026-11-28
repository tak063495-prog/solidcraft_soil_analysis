"""Read completed diagnostics; fixed-order observations are not new solves."""
import argparse
import copy
import math
from pathlib import Path
import numpy as np
from a03cause_binding import study, reference, fields
from a03cause_driver import write
from a03cause_metrics import compare, metric, norm, common


def load(folder, side):return list(fields.load_fields(folder/(side+'-fields.jsonl')))


def vector(value):
    a=np.asarray(value,float)
    return dict(l2=norm(a),maximum_absolute=float(np.max(np.abs(a),initial=0.)),
                resultant=[math.fsum(map(float,a[:,k])) for k in range(3)])


def order_assessment(x, y, z):
    a=np.asarray(x,float);b=np.asarray(y,float);c=np.asarray(z,float)
    d1=b-a;d2=c-b;n1=norm(d1);n2=norm(d2)
    return dict(coarse_to_middle_absolute_l2=n1,middle_to_fine_absolute_l2=n2,
                observed_order=math.log2(n1/n2) if n1>0 and n2>0 else None,
                difference_cosine=math.fsum(float(i)*float(j) for i,j in zip(d1.ravel(order='C'),d2.ravel(order='C')))/(n1*n2) if n1>0 and n2>0 else None,
                interpretation='roundoff-level/finite-grid diagnostic; no asymptotic certification or extrapolated acceptance')


def assess(run, destination):
    report=study.read(run/'driver-result.json')
    result=dict(passed=False,a03_complete=False,diagnostic_only=True,all_field_numerical_acceptance=False,
                fixed_order_formula='C-order time/GP-or-node/component; math.fsum; ||a-b||2/max(||b||2,1e-12)',cases={},reference_assessment={},initialization={},force_paths={},refusals=report['refusals'])
    names=sorted(report['cases'])
    for name in names:
        folder=run/name;native=load(folder,'native');independent=load(folder,'independent')
        compared=compare(native,independent,fields.CATEGORIES)
        result['cases'][name]=dict(fixed_order_comparison=compared,legacy_comparison=report['cases'][name]['comparison'],
                                  native_audits=report['cases'][name]['native'],independent=report['cases'][name]['independent'])
        if name.endswith('-n4'):
            base=study.read(folder/'native/run.json');prep=next(a for a in base['transition_audits'] if a['stage']==base['case']['stages'][0]['name'])
            model=reference.Model(*[study.read(p) for p in study.paths(folder)])
            selected=set(prep['selected_elements']);points=base['integration_points'];state=prep['prepared_state'];fresh=base['case']['stages'][0]['states']['UpperClay']
            selected_states=[s for p,s in zip(points,state['material_states']) if p['element_index'] in selected]
            assert all(study.exact(s)==study.exact(fresh) for s in selected_states)
            keyorder=compare(independent,load(run/(name+'-key-order'),'independent'),fields.CATEGORIES)
            assert all(v.get('maximum_absolute',v.get('mismatches',0))==0 for v in keyorder['fields'].values())
            p=np.asarray([model.gamma*(model.level-point['position_m'][2]) for point in points])
            ulps={}
            for direction in (-math.inf,math.inf):
                level=math.nextafter(model.level,direction)
                perturbed=[model.gamma*(level-point['position_m'][2]) for point in points]
                ulps[str(direction)]=dict(input_level_m=level,input_change_m=level-model.level,pore_pressure_change_kpa=metric(perturbed,p))
            result['initialization'][name]=dict(fresh_selected_points=len(selected_states),fresh_state_equals_input_exact=True,
                repeated_native_run_exact_except_elapsed=report['cases'][name]['native']['whole_exact_except_elapsed'],
                reordered_keys_native_full_run_exact_except_elapsed=report['cases'][name]['key_order_full_run_exact_except_elapsed'],
                reordered_keys_preparation_and_final_state_exact=True,independent_reordered_keys=keyorder,
                pressure_ulp_formula_only_not_equilibrium=ulps,pressure_definition='gamma_w*(explicit matched level-z); fully submerged only',
                seeded_effective_stress=fresh['effective_stress_compression_positive_kpa'],
                reversed_independent_volume_assembly=compare(independent,load(run/(name+'-reverse-volume'),'independent'),fields.CATEGORIES))
            # Measure force contributions on the actual saved full equilibrium path.
            full_ground=copy.deepcopy(model.ground);full_ground['initialization']['geostatic']['hydrostatic_face_sets']=['exterior']
            full_model=reference.Model(model.mesh,full_ground,model.case)
            rows=[]
            for row,f in native:
                active=f['active_elements'];buoyancy=np.zeros_like(f['u'])
                for e in np.flatnonzero(active):
                    q=model.ops[e][0];buoyancy[q['nodes'],2]+=model.gamma*(q['weight']@q['shape'])
                full_ext,_=full_model.external(active,0.)
                missing=full_ext['exterior_water_load']-f['exterior_water_load']
                balance=f['exterior_water_load']+f['joint_water_load']-f['pore_internal']-buoyancy
                physical=sum(f[k] for k in ('gravity_load','base_load','surcharge_load','exterior_water_load','joint_water_load'))
                rows.append(dict(row=row,components={k:vector(f[k]) for k in ('effective_internal','pore_internal','contact_internal','gravity_load','base_load','surcharge_load','exterior_water_load','joint_water_load','holding_load')},
                    physical_plus_holding_minus_applied=vector(physical+f['holding_load']-f['applied_load']),
                    raw_total_force_residual_includes_supports=vector(f['internal']-f['applied_load']),
                    buoyancy=vector(buoyancy),buoyancy_identity_defect=vector(balance),
                    omitted_exterior_force_observation_only=vector(missing),
                    defect_plus_omitted_exterior=vector(balance+missing),
                    active_elements=int(np.count_nonzero(active)),exterior_faces=int(np.count_nonzero(f['exterior_faces'])),wet_exterior_faces=int(np.count_nonzero(f['wet_exterior_faces']))))
            prepared=list(fields.load_fields(folder/'native-prepared.jsonl'))[0][1]
            physical=sum(prepared[k] for k in ('gravity_load','base_load','surcharge_load','exterior_water_load','joint_water_load'))
            result['force_paths'][name]=dict(rows=rows,birth_preparation=dict(holding=vector(prepared['holding_load']),
                holding_identity_defect=vector(prepared['internal']-physical-prepared['holding_load']),
                components={k:vector(prepared[k]) for k in ('effective_internal','pore_internal','contact_internal','gravity_load','base_load','surcharge_load','exterior_water_load','joint_water_load')}),
                maximum_buoyancy_identity_defect=max(r['buoyancy_identity_defect']['maximum_absolute'] for r in rows),
                maximum_defect_plus_omitted_exterior=max(r['defect_plus_omitted_exterior']['maximum_absolute'] for r in rows),
                omitted_exterior_is_force_accounting_not_a_new_solved_case=True)
            result['reference_assessment'][name+'-gl5-gl7']=compare(independent,load(run/(name+'-gl7'),'independent'),fields.CATEGORIES)
    for label in ('h1_5','h2_0'):
        for side in ('native','independent'):
            available=[int(n.split('-n')[1]) for n in names if n.startswith(label+'-n')]
            if not all(n in available for n in (2,4,8)):continue
            rows=[common(load(run/f'{label}-n{n}',side),2) for n in (2,4,8)]
            # compare checks exact coordinates and IDs inherited from each layout.
            values={'2_to_4':compare(rows[0],rows[1],fields.CATEGORIES),'4_to_8':compare(rows[1],rows[2],fields.CATEGORIES)}
            values['orders']={k:order_assessment(*[np.stack([f[k] for _,f in r]) for r in rows]) for k in sorted(rows[0][0][1]) if k not in fields.CATEGORIES}
            result['reference_assessment'][label+'-'+side+'-common-grid']=values
    mcc_folder=run/'refusals/fresh-MCC'
    if mcc_folder.exists():
        g=study.read(mcc_folder/'ground.json');s=g['initialization']['geostatic']['seed']['states']['UpperClay']
        sig=np.asarray(s['effective_stress_compression_positive_kpa'],float);p=math.fsum(map(float,sig[:3]))/3
        mat=next(m['mcc'] for m in g['materials'] if m['id']==2)
        result['mcc_seed_diagnosis']=dict(state=s,parameters=mat,initial_mean_effective_pressure_kpa=p,
            isotropic_pc_over_p=s['preconsolidation_pressure_kpa']/p,
            initial_bulk_modulus_kpa=(1+s['void_ratio'])*p/mat['kappa'],
            positive_pressure_seed_control_completed=True,guard_is_MC_only_before_state_validation=True,
            cannot_infer_MCC_birth_equilibrium_or_formation_history=True)
    result['additional_refusals']=report.get('additional_refusals',[])
    result['water_accounting']=dict(flux=None,storage=None,birth_water_transfer=None,balance_acceptance=None,
        reason='prescribed pressure has no water balance PDE or birth-transfer contract; null is not zero')
    result['fresh_completed_paths']=len(report['fresh_paths']);result['fresh_accepted_states']=sum(r['states'] for r in report['fresh_paths'])
    result['completed_flow_prefixes']=sum('completed_flow_prefix' in r for r in report['refusals'])
    result['valid_MCC_no_birth_controls']=sum('valid_MCC_seed_control' in r for r in report['refusals'])
    result['prefix_and_seed_control_paths_not_in_formal_matrix']=result['completed_flow_prefixes']+result['valid_MCC_no_birth_controls']
    result['prefix_and_seed_control_accepted_states']=result['completed_flow_prefixes']*8+result['valid_MCC_no_birth_controls']*4
    result['restored_checkpoints']=sum(len(c.get('restart',[])) for c in report['cases'].values())
    result['restored_suffix_states']=sum(r['suffix_states_exact'] for c in report['cases'].values() for r in c.get('restart',[]))
    result['rehash_tamper_refusals']=sum(len(c.get('tampering',[])) for c in report['cases'].values())
    write(destination,result)
    print(destination)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--run',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();assess(args.run,args.output)
