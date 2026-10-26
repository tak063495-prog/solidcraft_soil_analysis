"""Original-input dry K0/MCC + plate/frame independent construction history."""
import mcc_structural_reference as coupled  # Limits BLAS threads first.
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import fem_reference as ref
from beam_reference import System as Members
from mcc_projected_cycle_diagnostic import integrals
from mcc_rate_fem import seed, physical


def solve(mesh_file, ground_file, case_file, output, steps):
    if output.exists() or steps<2:raise ValueError('new output and at least two steps required')
    mesh,ground,case=[json.loads(p.read_text(encoding='utf-8-sig')) for p in (mesh_file,ground_file,case_file)]
    ground=ref.normalize_ground(ground);case=ref.normalize_stage_case(case)
    projected=case['schema']=='soil-fem-projected-mcc-stages/1'
    if case['schema'] not in ('soil-fem-projected-mcc-stages/1','soil-fem-mixed-mcc-stages/1'):
        raise ValueError('explicit mixed MCC formulation required')
    if ground['initialization']['method']!='horizontal_k0' or ground['initialization']['water_level_m'] is not None:
        raise ValueError('independent member construction currently requires dry horizontal K0')
    if any(case.get(k) for k in ('joint_system','anchor_system','initially_inactive_bodies')):
        raise ValueError('unsupported independent contact/water/birth coupling')
    paths=[mesh_file,ground_file,case_file]+[Path(__file__).with_name(name) for name in (
        'run_mcc_structural_reference.py','mcc_structural_reference.py','mcc_projected_sparse_reference.py',
        'mcc_rate_fem.py','mcc_rate_midpoint.py','fem_reference.py','beam_reference.py','plate_reference.py',
        'mcc_projected_cycle_diagnostic.py','verify_high_order.py')]
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    bindings={str(p.resolve()):sha(p) for p in paths}
    initial=seed(mesh,ground);state=initial.copy();u=np.zeros(np.asarray(mesh['nodes']).size)
    active=np.ones(len(mesh['tetrahedra']),bool);members=Members(mesh,case['beam_system']);history=members.initial_state()
    factor=0.;pressure=ground['initialization']['surface_surcharge_kpa']
    system=coupled.StructuralSystem(mesh,ground,active,case,projected=projected)
    system.configure_structures(members,history,0.)
    offset=np.zeros(len(system.keys));frames=[];releases=[]

    def external():
        out=physical(mesh,ground,active,pressure,case.get('pressure_on_original_surface_only',False))[0]
        out+=factor*ref.conditions_load(mesh,case.get('boundary') or {},active,skip_removed_surfaces=True)
        return out.reshape(-1)

    def response():return members.solve_interface(history,u.reshape(-1,3),factor)

    def total_internal(local,member):
        force=system.force(local)
        for node,value in member['interface_residual'].items():force[3*node:3*node+3]+=value
        return force

    def frame(name,diagnostics):
        nonlocal history
        member=response();history=member['state'];local=state[active].reshape(-1,14)
        balance=external()-total_internal(local,member)
        residual=float(np.linalg.norm(system.reduce(balance))/max(np.linalg.norm(external()),1.))
        defect,volume_scale,theta=integrals(system,state,initial,u)
        compatibility=float(np.linalg.norm(defect-offset)/max(np.linalg.norm(volume_scale),system.weight.sum()*1e-12))
        if residual>=1e-8 or compatibility>=1e-8:
            raise ValueError(f'independent nonlinear frame residual failed: {name}, force={residual}, volume={compatibility}')
        primary=system.pressure(local)[0]
        p=local[:,:3].mean(axis=1).reshape(system.weight.shape)
        weak_mean=float(np.linalg.norm(system.mass@primary-system.pressure(local)[1])/max(np.linalg.norm(system.pressure(local)[1]),1.))
        gp_gap=float(abs(np.einsum('eqa,ea->eq',system.n,primary[system.pi])-p).max())
        if not projected and gp_gap>1e-7:raise ValueError('independent GP mean stress left the P1 pressure space')
        frames.append(dict(name=name,displacement_m=u.reshape(-1,3).tolist(),material_state=state.tolist(),
            active_elements=active.tolist(),beams=history,
            primary_pressure=[dict(node_index=k[0],material_id=k[1],compression_positive_kpa=float(v)) for k,v in zip(system.keys,primary,strict=True)],
            relative_force_residual=residual,relative_frozen_volume_residual=compatibility,
            relative_weak_mean_pressure_residual=weak_mean,maximum_gp_pressure_difference_kpa=gp_gap,
            plastic_points=int(np.count_nonzero(np.linalg.norm(local[:,8:],axis=1)>1e-10)),
            diagnostics=diagnostics,linear_diagnostics=system.solver_diagnostics()))
        print(json.dumps(dict(stage=name,force_residual=residual,volume_residual=compatibility,
            plastic_points=frames[-1]['plastic_points'])),flush=True)

    frame('initial_k0',{})
    for stage in case['stages']:
        previous_state=state.copy();factor0=factor
        if stage['kind']=='boundary':factor=stage['target_factor']
        elif stage['kind']=='pressure':pressure=stage['target_pressure_kpa']
        elif stage['kind'] in ('install_beams','remove_beams'):
            history=members.change(history,u.reshape(-1,3),stage)
        elif stage['kind']=='excavate':
            remove=set()
            for body in stage['bodies']:
                selected=mesh.get('excavation_sets',{}).get(body)
                if selected is None:
                    regions={r['id'] for r in mesh['regions'] if r['body']==body}
                    selected=[i for i,r in enumerate(mesh['region_ids']) if r in regions]
                if not selected or not active[selected].all():raise ValueError('empty/overlapping reference excavation')
                remove.update(selected)
            active[sorted(remove)]=False
            pressure=stage['remaining_surface_pressure_kpa']
            system=coupled.StructuralSystem(mesh,ground,active,case,projected=projected)
            offset=integrals(system,state,initial,u)[0]
            used=set(np.unique(ref.cells(mesh)[active]))
            if any(members.attach[n] is not None and members.attach[n] not in used for n in np.flatnonzero(members.active_nodes(history))):
                raise ValueError('active reference structure lost its soil attachment')
        else:raise ValueError('unsupported independent structural construction stage')
        # The starting material state and member references are our own prior
        # solution. Vary physical soil/member loads linearly over the interval.
        member_at_start=members.solve_interface(history,u.reshape(-1,3),factor0)
        local=state[active].reshape(-1,14)
        before=total_internal(local,member_at_start)
        after=external()
        if stage['kind'] in ('excavate','install_beams','remove_beams'):
            holding=before-after
            unused=np.ones(len(mesh['nodes']),bool);unused[np.unique(ref.cells(mesh)[active])]=False
            holding.reshape(-1,3)[unused]=0.
            releases.append(dict(stage=stage['name'],holding_nodal_load_kn=holding.tolist()))
        system.configure_structures(members,history,factor-factor0)
        local,u,diagnostics=system.advance(local,u,after-before,steps,
            progress=lambda t,n,e:print(f'{stage["name"]}: {t:.6f}, steps {n}, events {e}',flush=True))
        state[active]=local.reshape(state[active].shape)
        if not np.array_equal(state[~active],previous_state[~active]):raise ValueError('removed reference material history changed')
        frame(stage['name'],diagnostics)
    if bindings!={str(p.resolve()):sha(p) for p in paths}:raise ValueError('independent reference sources/inputs changed')
    result=dict(schema='soil-fem-independent-mcc-structures/1',projected_volume=projected,
        reference_steps=steps,state_origin='original-input analytical K0/OCR/MCC reference void ratio; independent interval histories and installation references',
        structural_unknowns='all free translations/rotations/end releases retained in mixed sparse KKT; independent polynomial DKT and section-flexibility matrices',
        frames=frames,release_audits=releases,input_and_source_sha256=bindings,
        spatial_convergence_accepted=False,goal_complete=False,p6_complete=False)
    with output.open('x',encoding='utf-8') as stream:json.dump(result,stream,allow_nan=False)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('mesh','ground','case','output'):p.add_argument(name,type=Path)
    p.add_argument('--steps',type=int,default=100);a=p.parse_args()
    solve(*(getattr(a,n).resolve() for n in ('mesh','ground','case','output')),a.steps)
