"""Audit the already saved original proportional NC fold, without path solving.

Native matrices use saved real constitutive tangents. The independent replay
uses its own saved history and unchanged strain-driven material map. Only six
accepted points per local path are replayed. No prehistory is recomputed.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np
from fold_stability import point,judge,low_stress,audit_load,spectrum

ROOT=Path(__file__).resolve().parents[1]
PRIOR=ROOT/'interim/2026-11-09-nc-pseudo-arclength'
PATHS=('native-coarse','native-fine','independent-map-arc-coarse','independent-map-arc-fine','independent-map-arc-material8')


def read(path):return json.loads(path.read_text(encoding='utf-8'))
def write(path,value):path.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')


def modules():
    sys.path.insert(0,str(PRIOR))
    from native_arc import Native
    from independent_map_arc import Model
    return Native,Model


def linear_audit(systems):
    results=[]
    for eq in systems:
        a=np.asarray(eq['matrix']);b=np.asarray(eq['rhs']);x=np.asarray(eq['solution'])
        residual=a@x-b;n=len(b)-1
        augmented=float(np.linalg.norm(residual)/max(np.linalg.norm(b),1e-30))
        # Corrector force block has its own RHS. A tangent system has zero
        # force RHS, so report its homogeneous backward error separately.
        force=float(np.linalg.norm(residual[:n])/max(np.linalg.norm(b[:n]),1e-30)) if np.linalg.norm(b[:n]) else float(
            np.linalg.norm((a@x)[:n])/max(np.linalg.norm(a[:n])*np.linalg.norm(x),1e-30))
        results.append(dict(augmented=augmented,force_block=force,
                            kind=eq.get('kind','corrector'),actual_augmented=spectrum(a)))
    return results


def material_points(mesh,ground,q,offsets,active,states,native,b_matrices,v):
    params={m['id']:m['mcc'] for m in ground['materials']};records=[];active_index=0
    data=None if native else np.asarray(states).reshape(-1,14)
    for ei,quadrature in enumerate(q):
        for qi,position in enumerate(quadrature['position']):
            pi=int(offsets[ei])+qi
            if native:
                state=states[pi];stress=state['effective_stress_compression_positive_kpa'];void=state['void_ratio']
            else:
                state=dict(stress=data[pi,:6].tolist(),plastic_strain=data[pi,8:].tolist(),pc=float(np.exp(data[pi,6])))
                stress=data[pi,:6].tolist();void=float(np.expm1(data[pi,7]))
            row=dict(point=pi,element=ei,position_m=np.asarray(position).tolist(),active=bool(active[ei]),
                     stress_kpa=stress,void_ratio=void,kappa=params[mesh['material_ids'][ei]]['kappa'],history=state)
            if active[ei]:
                row['weak_mode_strain_norm']=float(np.linalg.norm(b_matrices[active_index]@v));active_index+=1
            records.append(row)
    return records


def run(evidence,output):
    output.mkdir(parents=True,exist_ok=False);Native,Model=modules()
    from verify_high_order import prepare,physical
    from low_preload_nc_reference import ClosedStickSystem
    native=Native(evidence/'prior',Path('read-only-no-binary-called'),output/'native-geometry-only')
    independent=Model(evidence/'independent-input-rk4-fine',evidence/'independent-arc-history400',4)
    q,offsets=prepare(native.mesh)
    # Preserve the original constrained-DOF bookkeeping as well as free force.
    from nc_fold_path import native_loads
    active=np.asarray(native.active)
    origin=native.checkpoint['run']['continuation']['origin']
    final,holding=native_loads(native.mesh,native.ground,native.case,origin,active)
    load_comparison=float(np.linalg.norm(holding-native.holding)/max(np.linalg.norm(holding),1.))
    for lam in (0.,native.lam0,1.):audit_load(native.load(lam),final,holding,lam)
    paths={};stress_results={}
    for name in PATHS:
        native_map=name.startswith('native');m=native if native_map else independent
        samples=read(evidence/name/'samples.json');events=read(evidence/name/'events.json');systems=read(evidence/name/'linear-systems.json')
        peak=max(range(len(samples)),key=lambda i:samples[i]['progress']);indices=range(max(0,peak-2),min(len(samples),peak+4))
        tangent_systems=[s for s in systems if s.get('kind')=='tangent']
        if not native_map and len(tangent_systems)!=len(samples):raise ValueError('missing accepted independent tangent systems')
        audited=[]
        for index in indices:
            row=samples[index];previous=samples[index-1] if index else None
            step=row['step'];history=previous['states' if native_map else 'material_state'] if previous else (
                native.cp['material_states'] if native_map else independent.oldfull)
            selected_events=[e for e in events if e['step']==step];last=selected_events[-1]
            if native_map:
                response=read(evidence/name/f'call-{last["call"]:05d}-result.json')
                if not all(p['ad_primal_exact'] and p['substep_replay_exact'] for p in response['points']):raise ValueError('native material replay equality failed')
                ds=np.array([p['implicit_tangent'] for p in response['points']]);b=native.rb;weights=native.w
                if any(p['updated_state']!=row['states'][pi] for p,pi in zip(response['points'],native.pi)):raise ValueError('saved native GP history mismatch')
                plane=np.array(previous['tangent']) if previous else np.r_[np.linalg.solve(native.k0,native.L),1.]
                plane/=np.linalg.norm(plane)
                actual=[eq for eq in systems if eq['step']==step]
                h=native.R.T@holding;ffinal=final;hfull=holding
            else:
                independent.count=8 if name.endswith('material8') else 4
                integrations=[];original=independent.material.integrate
                def capture(*a,**kw):
                    value=original(*a,**kw);integrations.append(value.copy());return value
                independent.material.integrate=capture
                try:
                    f,replayed_k,new=independent.evaluate(np.asarray(last['reduced_increment_m']),np.asarray(history))
                finally:independent.material.integrate=original
                if not np.array_equal(new,np.asarray(row['material_state'])):raise ValueError('all independent material histories do not replay exactly')
                ds=np.stack([(integrations[1+2*j][:,:6]-integrations[2+2*j][:,:6])/(2e-8) for j in range(6)],axis=2)
                b=independent.B;weights=independent.w
                plane=np.asarray(tangent_systems[index]['matrix'])[-1]/1000
                # The saved tangent system is exactly the converged real matrix.
                actual=[tangent_systems[index]]
                # Associate preceding correctors with this accepted tangent row.
                ti=[i for i,e in enumerate(systems) if e.get('kind')=='tangent'][index]
                start=0 if index==0 else [i for i,e in enumerate(systems) if e.get('kind')=='tangent'][index-1]+1
                actual=systems[start:ti+1]
                h=-independent.L;hfull=-independent.delta;ffinal=independent.before+independent.delta
            matrices=np.einsum('qai,qab,qbj,q->qij',b,ds,b,weights,optimize=True)
            # Native uses exactly the saved driver's global contraction; the
            # independent side uses its actual saved accepted tangent matrix.
            # The separately assembled per-GP sum must reproduce those matrices.
            replayed_global=np.einsum('qai,qab,qbj,q->ij',b,ds,b,weights,optimize=True)
            k=replayed_global if native_map else np.asarray(tangent_systems[index]['matrix'])[:-1,:-1].copy()
            if np.linalg.norm(replayed_global-k)/np.linalg.norm(k)>1e-10:raise ValueError('accepted independent tangent replay differs')
            if np.linalg.norm(matrices.sum(0)-k)/np.linalg.norm(k)>1e-10:raise ValueError('per-GP tangent decomposition differs from global matrix')
            pieces=[dict(kind='material_gp',point=int(pi),matrix=matrix) for pi,matrix in zip(native.pi if native_map else np.concatenate([np.arange(offsets[i],offsets[i+1]) for i in np.flatnonzero(independent.active)]),matrices)]
            present={p['point'] for p in pieces}
            pieces.extend(dict(kind='material_gp',point=i,active=False,matrix=np.zeros_like(k)) for i in range(int(offsets[-1])) if i not in present)
            # All 32 current contacts are inactive. Their nonzero stored traction
            # is frozen history, never a current force or stiffness contribution.
            pieces.extend(dict(kind='contact',point=i,active=False,matrix=np.zeros_like(k)) for i in range(32))
            d=point(k,h,plane,contributions=pieces)
            element_sums={}
            for contribution in d['contributions']:
                if contribution['kind']=='material_gp':
                    element=int(np.searchsorted(offsets,contribution['point'],side='right')-1)
                    element_sums[element]=element_sums.get(element,0.)+contribution['w_k_v']
            d['element_w_k_v']=element_sums
            d['material_gp_w_k_v_sum']=sum(c['w_k_v'] for c in d['contributions'] if c['kind']=='material_gp')
            d['contact_w_k_v_sum']=sum(c['w_k_v'] for c in d['contributions'] if c['kind']=='contact')
            audits=linear_audit(actual)
            # Re-solve the small accepted native tangent system, whose original
            # raw vector was not saved before normalization. Keep this distinct.
            tangent_a=np.vstack((np.column_stack((k,h)),plane));rhs=np.r_[np.zeros(len(k)),1.]
            solution=np.linalg.solve(tangent_a,rhs)
            replay_audit=linear_audit([dict(kind='tangent',matrix=tangent_a.tolist(),rhs=rhs.tolist(),solution=solution.tolist())])[0]
            augmented=max([a['augmented'] for a in audits]+[replay_audit['augmented']])
            force_block=max([a['force_block'] for a in audits]+[replay_audit['force_block']])
            states=row['states' if native_map else 'material_state']
            material=np.array([states[i]['effective_stress_compression_positive_kpa'] for i in native.pi]) if native_map else np.asarray(states)[independent.active].reshape(-1,14)[:,:6]
            internal=-np.einsum('qai,qa,q->i',b,material,weights)
            load=m.load(row['progress']);r=np.linalg.norm(m.R.T@load-internal)/max(np.linalg.norm(load),1.)
            loadaudit=audit_load(load,ffinal,hfull,row['progress'])
            v=np.asarray(d['right_weak_mode'])
            records=material_points(m.mesh,m.ground,q,offsets,active if native_map else independent.active,states,native_map,b,v)
            if index==peak:
                stress_results[name]=low_stress(records)
                stress_results[name]['lambda']=row['progress']
                stress_results[name]['contacts']=dict(count=32,active_count=0,current_stiffness_contribution=0,
                    frozen_history=native.cp['joints'] if native_map else independent.seed['contacts'])
            audited.append(dict(step=step,**{'lambda':row['progress']},accepted=True,diagnostic=d,
                                relative_force_residual=float(r),augmented_linear_residual=augmented,
                                force_block_linear_residual=force_block,full_component_load_audit=loadaudit,
                                all_gp_histories_replay_exact=True,saved_system_audits=audits,
                                accepted_tangent_small_matrix_replay=replay_audit,
                                actual_accepted_matrix_columns=k.T.tolist(),reduced_holding_load_kn=h.tolist(),
                                saved_dlambda_ds=row['tangent'][-1],
                                oriented_tangent_difference=float(np.linalg.norm(np.asarray(d['scaled_tangent'])-row['tangent']))))
        paths[name]=dict(tangent_source=('existing native MCC material-map implicit tangent used by archived arc driver' if native_map else 'existing independent RK4 material-map central derivative used by archived arc driver'),production_cli_newton_matrix=False,peak_lambda=samples[peak]['progress'],peak_sample_index=peak,
                         sampled_points=len(samples),local_arc_step_m=read(evidence/name/'controls.json')['step_size_m'],
                         diagnostic_points=audited,verdict=judge(audited))
    result=dict(schema='soil-fem-nc-fold-stability/1',conclusion='undetermined',conclusion_ja='未確定',
                passed=False,nc01_complete=False,physical_stability='unconfirmed',
                scope='fixed original proportional load; six saved accepted points around each local turning point',
                reason='local fold evidence is separated from a missing original proportional prehistory 1/128 vs 1/256 comparison',
                native_holding_reconstruction_relative_difference=load_comparison,paths=paths,
                comparison=dict(native_local_step_peak_difference=abs(paths['native-coarse']['peak_lambda']-paths['native-fine']['peak_lambda']),
                    independent_local_step_peak_difference=abs(paths['independent-map-arc-coarse']['peak_lambda']-paths['independent-map-arc-fine']['peak_lambda']),
                    independent_material_4_8_peak_difference=abs(paths['independent-map-arc-coarse']['peak_lambda']-paths['independent-map-arc-material8']['peak_lambda']),
                    native_independent_peak_difference=abs(paths['native-fine']['peak_lambda']-paths['independent-map-arc-fine']['peak_lambda'])),
                evidence_gaps=['original proportional initial-step 1/128 versus 1/256 with fixed prehistory not established by nonproportional series',
                               'dynamic/thermodynamic stability is not implied by a nonsymmetric static tangent',
                               'a local maximum does not prove absence of disconnected equilibria at lambda=1'],
                prehistory_recomputed=False,production_solver_changed=False,legacy_evidence_changed=False,
                original_archive_sha256=read(evidence.parent/'extraction.json')['archive_sha256'],
                base_tip='5f564621c7396c2a605edbfe5451cd800cf318d3',
                fixed_inputs_sha256={name:hashlib.sha256((evidence/'prior/prior/fresh-head-midpoint/input'/(name+'.json')).read_bytes()).hexdigest() for name in ('mesh','ground','stages')},
                fixed_prehistory_sha256={name:hashlib.sha256((evidence/name).read_bytes()).hexdigest() for name in ('prior/prior/fresh-head-midpoint/job/recovery-checkpoint-0.json','independent-input-rk4-fine/reference.json','independent-input-rk4-fine/unload.json')},
                production_cli_newton_fold_confirmed=False,
                source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'tools/fold_stability.py',ROOT/'tools/nc_fold_path.py',PRIOR/'native_arc.py',PRIOR/'independent_map_arc.py']})
    write(output/'fold-results.json',result);write(output/'low-stress.json',dict(schema='soil-fem-nc-low-stress/1',diagnostic_only=True,paths=stress_results))
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--evidence',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    a=parser.parse_args();r=run(a.evidence.resolve(),a.output.resolve());print(json.dumps(dict(conclusion=r['conclusion_ja'],comparison=r['comparison']),ensure_ascii=False,indent=2))
