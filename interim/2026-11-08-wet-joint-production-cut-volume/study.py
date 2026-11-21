"""Existing wet JOINT only. Production re-solves; saved independent references.

Measurements never certify all-history convergence. No NC or large runs.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import copy
import hashlib
import json
from pathlib import Path
import sys
import time
import zipfile
import numpy as np
import meshio
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import fem_reference as fem
from diagnose_cut_volume_history import inputs, FIELDS, summary
from diagnose_wet_joint_convergence import configuration, schedule
from diagnose_wet_joint_cut_volume import components, error
from diagnose_wet_joint_simultaneous import recorded, force_alignment, NATIVE
from export_stages import export
from hydrostatic_volume_audit import correction, enabled
from joint_water_reference import Water
from verify_checkpoint import invoke
from verify_high_order import prepare, internal, physical
from verify_joint_water import audit
from verify_wet_joint_cycle import restarts
from verify_wet_structures import without_time
from wet_joint_cycle_reference import native_fields, differences, serializable, contact_work

BASE = '840de37cfc8b5ec7bab1fdf3ffed691c32a209ad'
INPUTS = ('joint-water-column-tet10-mesh', 'ground-joint-water-mc', 'stages-joint-water-cycle')
ROUTES = (('original', None, None), ('g256-m5e-6', 256, 5e-6),
          ('g512-m5e-6', 512, 5e-6), ('g512-m2.5e-6', 512, 2.5e-6))

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def write(p, value):
    p.write_text(json.dumps(serializable(value), ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')

def solve(solver, folder):
    folder.mkdir(parents=True, exist_ok=False)
    raw = [json.loads((ROOT/'examples'/(n+'.json')).read_text()) for n in INPUTS]
    timings = []
    for key, steps, maximum in ROUTES:
        out = folder/key; out.mkdir()
        case = raw[2] if steps is None else configuration(raw[2], steps, maximum)
        for n, value in zip(('mesh', 'ground', 'stages'), (raw[0], raw[1], case)):
            write(out/(n+'.json'), value)
        start = time.monotonic()
        invoke([str(solver), 'solve-stages', *[str(out/(n+'.json')) for n in ('mesh','ground','stages')], str(out/'stage-results.json')], out/'solve.log')
        run = json.loads((out/'stage-results.json').read_text())
        assert enabled(run) and run['complete'] and len(run['frames']) == 7, run.get('failure_reason')
        assert run['case']['settings']['relative_force_tolerance'] == 1e-8
        assert run['case']['settings']['relative_linear_tolerance'] == 1e-10
        if steps is not None: schedule(run, steps)
        timings.append(dict(route=key, elapsed_seconds=time.monotonic()-start, complete=True))
        write(folder/'timings.json', timings)
        print('completed', timings[-1], flush=True)
    # Every phase, strictly interior, sealed by this unchanged executable.
    run = json.loads((folder/'g256-m5e-6/stage-results.json').read_text())
    rows = restarts(solver, folder/'g256-m5e-6', run)
    write(folder/'restarts.json', rows)
    write(folder/'solver.json', dict(sha256=sha(solver), path=str(solver)))

def fields_and_audit(mesh, ground, run, folder):
    """15 stored fields + kinematic strain, displayed GP pressure, force reaction.

    Nodal reaction is reconstructed from stored stress/contact and exact pressure,
    not a separately saved native primary solution. Retain its free residual.
    """
    assert enabled(run)
    q, offset = prepare(mesh, run['integration_points'])
    water = Water(mesh, ground, run['case'], run['joint_integration_points'])
    result = native_fields(run)
    rows = []
    for index, (frame, target) in enumerate(zip(run['frames'], result['frames'])):
        state = frame['state']; active = np.asarray(state['active_elements']); u = np.asarray(state['displacement_m'])
        assert active.all() and state.get('pore_pressure_kpa') is None
        level = state.get('hydrostatic_water_level_m', ground['initialization']['water_level_m'])
        assert level == .75
        pore = np.maximum(level-np.asarray([p['position_m'][2] for p in run['integration_points']]),0)*10
        target['pore'] = pore
        target['strain'] = np.concatenate([np.einsum('qai,i->qa', v['b'], u[v['nodes']].ravel()) for v in q])
        force, _ = internal(mesh, ground, state['material_states'], active, q, offset, hydrostatic_cut=True)
        interface = water.hydrostatic(active, level)
        force += water.audit(state, frame['joint_results'], interface)
        load, _ = physical(mesh, ground, active, state['surface_pressure_kpa'], run['case']['pressure_on_original_surface_only'],joint_system=water.contact)
        load += interface['load']+state['boundary_factor']*fem.conditions_load(mesh, run['case']['boundary'], active)
        basis, constraint_error = fem.staged_constraint_space(mesh, run['case'], state)
        residual = force-load
        reaction = residual.copy()
        for n, b in basis.items(): reaction[n] -= b@(b.T@residual[n])
        # Match the reference definition: internal minus applied nodal load,
        # including the small free residual. Check projected supports separately.
        target['reaction'] = residual
        reduced = np.concatenate([b.T@residual[n] for n,b in basis.items()])
        relative = float(np.linalg.norm(reduced)/max(np.linalg.norm(load),1))
        assert relative < 1e-8 and constraint_error < 1e-10
        np.testing.assert_allclose(reaction.sum(axis=0), frame['support_reaction_kn'], atol=1e-7, rtol=1e-9)
        display_error = None
        if folder is not None:
            gp = meshio.read(folder/f'integration-points-{index:03}.vtu'); ids = gp.point_data['integration_point_id']
            np.testing.assert_allclose(gp.point_data['pore_pressure_kpa'], pore[ids], atol=1e-11, rtol=1e-12)
            display_error = float(np.max(abs(gp.point_data['pore_pressure_kpa']-pore[ids])))
        assert set(target)-{'name'} == set(FIELDS)
        rows.append(dict(frame=frame['name'], independent_cut_force_relative_residual=relative,
                         constraint_error_m=constraint_error, displayed_gp_pressure_error_kpa=display_error,
                         physical_interface_pressure_error_kpa=float(np.max(abs(np.asarray(target['physical_pressure'])-2.5)))))
    return result, rows

def holding_audit(mesh, ground, run):
    q, offset = prepare(mesh, run['integration_points']); g = run['geostatic']; seed = g['seed']
    water = Water(mesh, ground, run['case'], run['joint_integration_points']); active=np.asarray(seed['active_elements'])
    interface=water.hydrostatic(active,.75)
    # Seed contact force is zero. The independent effective seed + cut pressure
    # is assembled before looking at the stored holding load.
    force,_=internal(mesh,ground,seed['material_states'],active,q,offset,hydrostatic_cut=True)
    base,_=physical(mesh,ground,active,seed['surface_pressure_kpa'],run['case']['pressure_on_original_surface_only'],joint_system=water.contact)
    base+=interface['load']
    expected=(force-base).ravel(); saved=np.asarray(g['initial_holding_load_kn'])
    delta=error(saved,expected)
    assert delta['maximum_absolute'] < 1e-10
    np.testing.assert_allclose(seed['applied_nodal_load_kn'], force.ravel(),atol=1e-10,rtol=0)
    # At each restarted geostatic point verify the ORIGINAL saved holding load
    # times (1-progress), rather than inventing loads from the current stress.
    return dict(initial_holding_difference=delta, seed_internal_difference=error(seed['applied_nodal_load_kn'],force.ravel()),
                excavation_stages_present=False, final_holding_zero=not np.asarray(g['last_committed']['geostatic_release_load_kn']).any())

def native_gates(run):
    attempts=run['geostatic']['attempts']+run['attempts']
    accepted=[a for a in attempts if a['accepted']]
    linear=[s for a in accepted for s in a['linear_solves']]
    assert all(a['relative_force_residual'] < 1e-8 for a in accepted)
    assert all(s['true_relative_residual'] <= 1e-10 for s in linear)
    assert all(a['committed_state_unchanged_on_rejection'] for a in attempts if not a['accepted'])
    return dict(accepted=len(accepted),rejected=len(attempts)-len(accepted),linear_solves=len(linear),
                maximum_relative_force_residual=max(a['relative_force_residual'] for a in accepted),
                maximum_true_relative_linear_residual=max(s['true_relative_residual'] for s in linear),
                maximum_frame_relative_force_residual=max(f['relative_force_residual'] for f in run['frames']),
                all_rejections_preserve_committed_state=True)

def restart_state_audit(folder,mesh,ground,run,rows):
    q,offset=prepare(mesh,run['integration_points']);water=Water(mesh,ground,run['case'],run['joint_integration_points'])
    checks=[]
    for row in rows:
        p=folder/f"checkpoint-{row['stage']}.json";assert sha(p)==row['checkpoint_sha256']
        saved=json.loads(p.read_text());r=saved['run'];state=r['last_committed'];active=np.asarray(state['active_elements'])
        assert enabled(r) and state.get('pore_pressure_kpa') is None
        assert state.get('hydrostatic_water_level_m',ground['initialization']['water_level_m'])==.75
        pressure_error=float(np.max(abs(np.asarray(state['joint_water']['physical_pressure_kpa'])-2.5)))
        assert pressure_error < 1e-10
        force,_=internal(mesh,ground,state['material_states'],active,q,offset,hydrostatic_cut=True)
        traction=np.asarray([v['traction_kpa'] for v in state['joints']['points']])
        force+=np.einsum('qag,qa,q->g',water.contact.c,traction,water.contact.weight).reshape(-1,3)
        load=np.asarray(state['applied_nodal_load_kn']).reshape(-1,3)
        basis,ce=fem.staged_constraint_space(mesh,run['case'],state)
        residual=np.concatenate([b.T@(force-load)[n] for n,b in basis.items()])
        relative=float(np.linalg.norm(residual)/max(np.linalg.norm(load),1))
        assert relative < 1e-8 and ce < 1e-10
        check=dict(stage=row['stage'],relative_force_residual=relative,interface_pressure_error_kpa=pressure_error)
        if row['stage']=='geostatic':
            holding=np.asarray(run['geostatic']['initial_holding_load_kn'])
            interface=water.hydrostatic(active,.75)
            base,_=physical(mesh,ground,active,state['surface_pressure_kpa'],run['case']['pressure_on_original_surface_only'],joint_system=water.contact)
            base+=interface['load']
            expected=base.ravel()+(1-row['progress'])*holding
            check['saved_holding_path_difference']=error(load.ravel(),expected)
            assert check['saved_holding_path_difference']['maximum_absolute'] < 1e-10
        checks.append(check)
    return checks

def saved_references():
    prior=ROOT/'interim/2026-11-08-wet-joint-simultaneous-force'
    values, provenance=recorded(prior,['paths/g1024-l32.json','paths/g2048-l32.json'])
    results={Path(n).stem:v for n,v in values.items()}
    prior=ROOT/'interim/2026-11-08-wet-joint-finer-full-history'; manifest=json.loads((prior/'archive-sha256.json').read_text())
    wanted={'paths/g2048-l64/result.json':'g2048-l64','paths/g4096-l64/result.json':'g4096-l64'}
    used={}; archives={}
    for row in manifest['archives']:
        p=prior/row['file']; assert sha(p)==row['sha256'],('archive changed',p)
        archives[row['file']]=row['sha256']
        with zipfile.ZipFile(p) as z:
            for name in set(z.namelist())&set(wanted):
                raw=z.read(name); digest=hashlib.sha256(raw).hexdigest();assert digest==manifest['entries'][name]
                results[wanted[name]]=json.loads(raw);used[name]=digest
    assert len(results)==4
    return results, dict(simultaneous=provenance,finer=dict(archives=archives,entries=used,rerun=False))

def dissipation_definition(mesh,ground,case,refs):
    system=Water(mesh,ground,case).contact
    old=-.0001*system.normal
    q=-100*system.normal; target=old.copy();target[:,0]+=.001
    rows=[]
    for n in (8,16,32,64):
        pd,nd,work=contact_work(system,old,q,target,n)
        rows.append(dict(partitions=n,maximum_numerical_kpa_m=float(nd.max()),sum_numerical_kpa_m=float(nd.sum()),
                         sum_plastic_kpa_m=float(pd.sum()),sum_endpoint_work_kpa_m=float(work.sum())))
    pairs={}
    for a,b in [('g2048-l32','g2048-l64'),('g2048-l64','g4096-l64')]:
        r=differences(refs[a],refs[b]);pairs[a+'__'+b]=[dict(frame=v['frame'],**v['numerical_dissipation']) for v in r]
    return dict(definition='dw = q_end dot (closed_fraction * delta_jump); nd = dw - delta_elastic_energy - plastic_dissipation; accumulate max(nd,0) for every contact subsegment',
                source='tools/wet_joint_cycle_reference.py:contact_work',source_sha256=sha(ROOT/'tools/wet_joint_cycle_reference.py'),
                is_continuum_dissipation=False,partition_dependent_by_definition=True,
                fixed_endpoint_partition_ablation=rows,recorded_reference_differences=pairs,
                interpretation='Approximately 100% relative differences may occur as this partition error halves toward zero. This does not certify any other contact history.')

def measure(folder, do_export=False):
    mesh,ground,case=inputs();native={};checks={};schedules={};holds={}
    for key,steps,maximum in ROUTES:
        out=folder/key;run=json.loads((out/'stage-results.json').read_text())
        assert run['ground']==ground
        assert fem.normalize_stage_case(json.loads((out/'stages.json').read_text()))==run['case']
        if do_export:export(out/'mesh.json',out/'stage-results.json',out/'exports')
        checks[key]=dict(stored_field_audit=audit(mesh,ground,run,out/'exports'),
                        force_and_gp=None,accepted_gates=native_gates(run))
        native[key],checks[key]['force_and_gp']=fields_and_audit(mesh,ground,run,out/'exports')
        write(out/'fields-18.json',native[key]);holds[key]=holding_audit(mesh,ground,run)
        if steps:schedules[key]=schedule(run,steps)
    refs,provenance=saved_references();comparisons={}
    reference_sensitivity={}
    for a,b in [('g1024-l32','g2048-l32'),('g2048-l32','g2048-l64'),('g2048-l64','g4096-l64')]:
        rows=differences(refs[a],refs[b]);reference_sensitivity[a+'__'+b]=dict(fields=rows,maxima=summary(rows))
    for key,fields in native.items():
        comparisons[key]={}
        for ref,value in refs.items():
            rows=differences(fields,value);comparisons[key][ref]=dict(fields=rows,maxima=summary(rows))
    old,old_provenance=recorded(NATIVE,['separated/g256-m5e-06/stage-results.json'])
    original_alignment=force_alignment(mesh,ground,case,next(iter(old.values())))
    comp={k:components(mesh,ground,rule,7) for k,rule in [('sampled','gj'),('cut','cut')]}
    sensitivities={}
    for a,b in [('g256-m5e-6','g512-m5e-6'),('g512-m5e-6','g512-m2.5e-6')]:
        rows=differences(native[a],native[b]);sensitivities[a+'__'+b]=dict(fields=rows,maxima=summary(rows))
    rest=json.loads((folder/'restarts.json').read_text());assert len(rest)==7 and all(v['whole_run_exact_except_elapsed_seconds'] for v in rest)
    return dict(schema='wet-joint-production-cut-volume/1',starting_commit=BASE,passed=False,history_converged=False,
                production_code_corrected=True,limited_case_completed=True,full_history_acceptance=False,
                solver=json.loads((folder/'solver.json').read_text()),timings=json.loads((folder/'timings.json').read_text()),
                input_sha256={n:sha(ROOT/'examples'/(n+'.json')) for n in INPUTS},
                root_cause='The uncut material GP rule samples max(H-z,0) across partially wet affine cells. Analytic water pressure force requires a different equilibrium solution and a holding load assembled with that same rule.',
                original_saved_force_alignment=original_alignment,original_native_provenance=old_provenance,
                pressure_force_difference=error(comp['sampled']['pore_force'],comp['cut']['pore_force']),
                native_checks=checks,native_schedules=schedules,holding_audits=holds,restarts=rest,
                restart_state_audits=restart_state_audit(folder/'g256-m5e-6',mesh,ground,json.loads((folder/'g256-m5e-6/stage-results.json').read_text()),rest),
                native_vs_independent=comparisons,native_step_sensitivity=sensitivities,independent_step_sensitivity=reference_sensitivity,independent_provenance=provenance,
                field_sources=dict(stored=15,strain='derived from saved displacement at original GP',pore='original GP hydrostatic formula and VTU readback',reaction='internal minus applied nodal load, reconstructed from saved effective stress/contact + exact hydrostatic force, including free residual as in the reference; projected support sum checked against saved native reaction'),
                water_transition_scope='This unchanged cycle has fixed H=.75 and CONTINUOUS pressure=2.5 throughout all boundary stages and restarts; it contains no changing-water/flow/excavation stage. Existing regression tests are separate evidence, not an independent 3D certificate for those paths.',
                numerical_dissipation=dissipation_definition(mesh,ground,case,refs),
                missing_acceptance=['No existing all-history error criterion','No independent general plastic/structure/Darcy/Richards/partial-crack or BIRTH 3D evidence','Excavation holding under this cut-water cycle is not exercised because the original input has no excavation'],
                nc_rerun=False,large_runs_executed=False,windows_ci='unverified',p3b_p4_complete=False,structure_water_joint_complete=False,p6_complete=False)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--solver',type=Path);p.add_argument('--output',type=Path,required=True);p.add_argument('--remeasure',type=Path);a=p.parse_args()
    if a.remeasure:
        result=measure(a.remeasure.resolve());write(a.output,result)
    else:
        solve(a.solver.resolve(),a.output.resolve());result=measure(a.output.resolve(),True);write(a.output/'measurements.json',result)
    print(json.dumps({k:result[k] for k in ('passed','history_converged','limited_case_completed')}))
