from common import HERE
from linear import diagnose
from load_split import diagnose_sources
"""Input-only GL125 soil + polynomial Tri6 + rigid-arm/frame direct solver.

Never imports production, native states, checkpoints or native accepted paths.
Existing independent element libraries are reused without modification.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import copy
import json
from pathlib import Path
import numpy as np
from scipy.linalg import lu_factor, lu_solve
from cases import ROOT, read
import fem_reference as ref
from plate6_offset_reference import System
from verify_high_order import physical


def solve(folder, output):
    m, g, c = [read(folder / (s+'.json')) for s in ('mesh', 'ground', 'stages')]
    g = ref.normalize_ground(g)
    assert not c.get('joint_system') and not c.get('anchor_system')
    assert g['initialization']['water_level_m'] is None
    assert all(s['kind'] in ('install_structures', 'remove_structures', 'boundary') for s in c['stages'])
    nstep = round(1 / c['settings']['initial_increment'])
    assert nstep * c['settings']['initial_increment'] == 1
    xyz = np.asarray(m['nodes']); nd = xyz.size
    if nd > 6500:
        raise ValueError('small/medium independent direct soil limit 6500 DOFs')
    active = np.ones(len(m['tetrahedra']), bool)
    mat = {v['id']: v['mohr_coulomb'] for v in g['materials']}
    regions = {r['id']: r['body'] for r in m['regions']}
    seeds = g['initialization']['geostatic']['seed']['states']
    k = np.zeros((nd, nd)); seed = np.zeros(nd); ops = []
    for ei in range(len(active)):
        d = ref.elastic_tensor(mat[m['material_ids'][ei]]['young_kpa'], mat[m['material_ids'][ei]]['poisson'])
        qp = ref.volume(m, ei, independent=True)
        idx = (qp['nodes'][:, None]*3+np.arange(3)).ravel()
        s0 = np.asarray(seeds[regions[m['region_ids'][ei]]]['effective_stress_compression_positive_kpa'])
        assert not s0.any(), 'this small case has zero input seed stress'
        k[np.ix_(idx, idx)] += np.einsum('qai,ab,qbj,q->ij', qp['b'], d, qp['b'], qp['weight'], optimize=True)
        seed[idx] -= np.einsum('qai,a,q->i', qp['b'], s0, qp['weight'], optimize=True)
        # Output at every production quadrature coordinate, assembled from input
        # rules; solve stiffness remains independently integrated with GL125.
        qout = ref.volume(m, ei)
        ops.append((idx, qout, d, s0, seeds[regions[m['region_ids'][ei]]]['void_ratio']))
    bottom = np.asarray(m['node_sets']['bottom'])
    fixed = (bottom[:, None]*3 + np.arange(3)).ravel()
    free = np.setdiff1d(np.arange(nd), fixed)
    soil_lu = lu_factor(k[np.ix_(free, free)])
    fsoil, _ = physical(m, g, active, g['initialization']['surface_surcharge_kpa'], c.get('pressure_on_original_surface_only', False))
    model = System(m, c['beam_system']); state = model.initial_state()
    u = np.zeros_like(xyz); factor = 0.; count = 0
    def observation(si, progress, f, br, ds, holding):
        nonlocal count
        count += 1
        stress = []; strain = []; void = []
        for idx, qp, d, s0, e0 in ops:
            eps = np.einsum('qai,i->qa', qp['b'], u.ravel()[idx])
            strain.extend(eps.tolist()); stress.extend((s0-eps@d.T).tolist())
            void.extend(((1+e0)*np.exp(eps[:, :3].sum(axis=1))-1).tolist())
        soil_internal = k@u.ravel()+seed
        interface = np.zeros(nd)
        for n, v in br['interface_residual'].items():
            interface[3*n:3*n+3] += v+br['attached_load'][n]
        total = soil_internal+interface
        applied = f.ravel().copy()
        for n, v in br['attached_load'].items(): applied[3*n:3*n+3] += v
        reaction = np.zeros(nd); reaction[fixed] = (total-applied)[fixed]
        residual = total-applied
        fb = model.assemble(state, factor)[2]
        attached = sum(float(v@u[n]) for n, v in br['attached_load'].items())
        return dict(accepted_index=count, stage_index=si, progress=progress,
                    displacement_m=u.tolist(), beams=copy.deepcopy(state),
                    stress_kpa=stress, engineering_strain=strain, void_ratio=void,
                    soil_internal_kn=soil_internal.tolist(), total_internal_kn=total.tolist(),
                    applied_nodal_load_kn=applied.tolist(), transition_release_load_kn=holding.ravel().tolist(),
                    support_reaction_by_node_kn=reaction.reshape(-1, 3).tolist(),
                    free_residual_relative=float(np.linalg.norm(residual[free])/max(np.linalg.norm(applied), 1.)),
                    load_work_kn_m=float(applied@u.ravel()+fb@ds-attached),
                    structural_load_work_kn_m=float(fb@ds),
                    beam_results={k: v for k, v in br.items() if k in ('plates','members','external_force','external_moment','support_force','support_moment')})
    def encoder(v):
        if isinstance(v, np.ndarray): return v.tolist()
        if isinstance(v, np.generic): return v.item()
        raise TypeError(type(v).__name__)
    with output.open('x', encoding='utf-8') as out:
        for j in range(1, nstep+1):
            lam = j/nstep
            u.ravel()[free] = lu_solve(soil_lu, (lam*fsoil.ravel()-seed)[free])
            br = model.solve_interface(state, u, 0.)
            # Inactive structural nodes remain frozen in the existing state.
            ds = np.c_[state['node_displacement_m'], state['node_rotation_rad']].ravel()
            rec = observation(0, lam, lam*fsoil, br, ds, np.zeros(nd))
            out.write(json.dumps(rec, default=encoder, allow_nan=False)+'\n')
        for si, stage in enumerate(c['stages'], 1):
            old_factor = factor; target = stage.get('target_factor', factor)
            transition = stage['kind'] != 'boundary'
            if transition:
                state = model.change(state, u, stage)
                at_origin = model.solve_interface(state, u, factor)
                holding = k@u.ravel()+seed-fsoil.ravel()
                for n, v in at_origin['interface_residual'].items(): holding[3*n:3*n+3] += v
            else: holding = np.zeros(nd)
            # Full explicit soil/frame unknowns, unlike production condensation.
            # Translational soil binding is injected directly on the free indices.
            active_nodes = model.active_nodes(state)
            binding = np.zeros((model.size, len(free)))
            reverse = {int(v): i for i, v in enumerate(free)}
            for i in np.flatnonzero(active_nodes):
                n = model.attach[i]
                if n is not None:
                    for a in range(3):
                        if 3*n+a in reverse: binding[6*i+a, reverse[3*n+a]] = 1.
            offset, _, qb, normals = model.space(state, u)
            # Attached translations are unknowns in binding, not an input solution.
            for i in np.flatnonzero(active_nodes):
                if model.attach[i] is not None: offset[6*i:6*i+3] = 0.
            kb, rest, *_ = model.assemble(state, old_factor)
            kk = np.block([[k[np.ix_(free, free)]+binding.T@kb@binding, binding.T@kb@qb],
                           [qb.T@kb@binding, qb.T@kb@qb]])
            factored = lu_factor(kk)
            source_probes = diagnose_sources(model,state,binding,qb,factored,m,c,active,free) if si == 2 else None
            for j in range(1, nstep+1):
                lam = j/nstep; factor = old_factor+(target-old_factor)*lam
                kb, rest, fb, eff, extf, extm, info = model.assemble(state, factor)
                f = fsoil + (factor*ref.conditions_load(m, c['boundary'], active))
                h = (1-lam)*holding
                f = f+h.reshape(-1, 3)
                rhs_b = fb-rest-kb@offset
                rhs = np.r_[(f.ravel()-seed)[free]+binding.T@rhs_b, qb.T@rhs_b]
                answer = lu_solve(factored, rhs)
                u.ravel()[free] = answer[:len(free)]
                ds = offset+binding@answer[:len(free)]+qb@answer[len(free):]
                br = model.results(state, ds, kb, rest, fb, eff, extf, extm, info, normals)
                state = br['state']
                rec = observation(si, lam, f, br, ds, h)
                if source_probes is not None and lam == 1.: rec['source_probes'] = source_probes
                rec['linear_diagnostics'] = diagnose(kk, rhs, answer, factored)
                physical_residual = np.asarray(rec['total_internal_kn'])-np.asarray(rec['applied_nodal_load_kn'])
                direct_residual = kk@answer-rhs
                rec['linear_diagnostics']['physical_minus_direct_soil_relative'] = float(np.linalg.norm(physical_residual[free]-direct_residual[:len(free)])/max(np.linalg.norm(rhs),1.))
                rec['direct_linear_residual'] = float(np.linalg.norm(kk@answer-rhs)/max(np.linalg.norm(rhs), 1.))
                out.write(json.dumps(rec, default=encoder, allow_nan=False)+'\n')
                out.flush()


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('folder', type=Path); p.add_argument('output', type=Path)
    a = p.parse_args(); solve(a.folder, a.output)

