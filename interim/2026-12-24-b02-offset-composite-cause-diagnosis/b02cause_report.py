"""Read-only whole-field, force-path and work accounting; no acceptance rewrite."""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from b02cause_binding import legacy
from b02cause_driver import write
from b02cause_metrics import common, dot, metric, norm

System = legacy['reference'].System
read = legacy['cases'].read


def rows(folder, kind):
    source = folder/('native/accepted.jsonl' if kind == 'native' else kind+'.jsonl')
    return [json.loads(line) for line in source.read_text(encoding='utf-8').splitlines()]


def factor(case, row):
    value = 0.
    si = row['stage_index']
    for stage in case['stages'][:max(0, si-1)]:
        value = stage.get('target_factor', value)
    if si and case['stages'][si-1]['kind'] == 'boundary':
        value += row['progress']*(case['stages'][si-1]['target_factor']-value)
    return value


def sum_vectors(vectors, size):
    return np.array([math.fsum(float(v[i]) for v in vectors) for i in range(size)])


def state_fields(row):
    s = row.get('state', row)
    return dict(u=s['displacement_m'], rotation=s['beams']['node_rotation_rad'],
                plate_reference=s['beams']['plate_reference_global_kinematics'],
                beam_reference=s['beams']['reference_global_kinematics'],
                stress=([v['effective_stress_compression_positive_kpa'] for v in s['material_states']]
                        if 'state' in row else row['stress_kpa']))


def compare_paths(a, b):
    aa, bb = common(a), common(b)
    if aa.keys() != bb.keys():
        raise ValueError('path coordinate mismatch')
    fa = [state_fields(aa[k]) for k in sorted(aa)]
    fb = [state_fields(bb[k]) for k in sorted(bb)]
    return {name: metric([v[name] for v in fa], [v[name] for v in fb]) for name in fa[0]}


def force_work_energy(folder, kind):
    mesh, case = [read(folder/(name+'.json')) for name in ('mesh', 'stages')]
    model = System(mesh, case['beam_system'])
    nd = 3*len(mesh['nodes'])
    fixed = {3*n+a for n in mesh['node_sets']['bottom'] for a in range(3)}
    free = [i for i in range(nd) if i not in fixed]
    result = []
    for row in rows(folder, kind):
        s = row.get('state', row)
        bs = s['beams']
        d = np.c_[bs['node_displacement_m'], bs['node_rotation_rad']].ravel()
        u = np.asarray(s['displacement_m']).ravel()
        lam = factor(case, row)
        k, rest, load, eff, _, _, info = model.assemble(bs, lam)
        plate_terms, beam_terms, plate_loads, beam_loads = [], [], [], []
        energy, work_identity = [], []
        for i, (nodes, section, kernel) in enumerate(model.plate_data):
            if not bs['active_plates'][i]:
                continue
            indices = np.array([6*n+a for n in nodes for a in range(6)])
            kin = d[indices]-bs['plate_reference_global_kinematics'][i]
            f = np.zeros(model.size)
            f[indices] = kernel.stiffness@kin
            plate_terms.append(f)
            plate = model.plates[i]
            traction = np.asarray(plate.get('constant_traction_kpa', [0.]*3), float)+lam*np.asarray(plate.get('boundary_traction_kpa', [0.]*3))
            traction[2] -= section['unit_weight_kn_m3']*section['thickness_m']
            fp = kernel.neutral.load(traction)
            mapped = kernel.map.T@fp
            source = np.zeros(model.size)
            source[indices] = mapped
            plate_loads.append(source)
            work_identity.append(dot(mapped, kin)-dot(fp, kernel.map@kin))
            local_s = kernel.neutral.transform@kin
            local_p = kernel.neutral.transform@kernel.map@kin
            terms = [[] for _ in range(7)]
            strains = []
            for bary, weight in kernel.neutral.quadrature:
                bm, bb, bq, bd = kernel.neutral.operators(bary)
                es, ep, curvature = bm@local_s, bm@local_p, bb@local_p
                strains.append(float(np.max(np.abs(ep-(es-kernel.e*curvature)))))
                values = [.5*dot(es, kernel.neutral.membrane@es),
                          -kernel.e*dot(es, kernel.neutral.membrane@curvature),
                          .5*kernel.e**2*dot(curvature, kernel.neutral.membrane@curvature),
                          .5*dot(curvature, kernel.neutral.bending@curvature),
                          .5*kernel.neutral.shear*dot(bq@local_p, bq@local_p),
                          .5*kernel.neutral.drilling*float(bd@local_p)**2,
                          .5*dot(ep, kernel.neutral.membrane@ep)]
                for target, value in zip(terms, values, strict=True):
                    target.append(weight*value)
            sums = [math.fsum(values) for values in terms]
            energy.append(dict(id=plate['id'], attachment_membrane=sums[0], cross=sums[1], arm_squared=sums[2],
                               neutral_bending=sums[3], shear=sums[4], drilling=sums[5], neutral_membrane=sums[6],
                               strain_identity_max=max(strains),
                               membrane_identity_residual=math.fsum(sums[:3])-sums[6]))
        for mi in sorted(info):
            p, local_k, ref, pre, f, *_ = info[mi]
            beam_terms.append(p.T@(local_k@(p@d-ref)+pre))
            beam_loads.append(p.T@f)
        plate_internal = sum_vectors(plate_terms, model.size)
        beam_internal = sum_vectors(beam_terms, model.size)
        plate_load = sum_vectors(plate_loads, model.size)
        beam_load = sum_vectors(beam_loads, model.size)
        nodal_load = load-plate_load-beam_load
        structural_internal = sum_vectors([plate_internal, beam_internal], model.size)
        mapped_internal = np.zeros(nd)
        mapped_load = np.zeros(nd)
        # Structure node ID order is the declared immutable order, not hash order.
        for i, node in enumerate(model.attach):
            if node is not None:
                mapped_internal[3*node:3*node+3] += structural_internal[6*i:6*i+3]
                mapped_load[3*node:3*node+3] += eff[i]
        soil = np.asarray(row['soil_internal_kn'])
        applied = np.asarray(s['applied_nodal_load_kn'])
        holding = np.asarray(s['transition_release_load_kn'])
        if holding.size == 0:
            holding = np.zeros(nd)
        rhs_without_attachment = applied-mapped_load
        global_residual = soil+mapped_internal-applied
        # This saved-field reassembly is coupled to input-only full solves elsewhere.
        # It is not a new solution, and does not stand in for global equilibrium.
        structure_free = model.space(bs, np.asarray(s['displacement_m']))[2]
        rotation_residual = structure_free.T@(structural_internal-load)
        source_work = dict(soil_physical=dot(rhs_without_attachment-holding, u),
                           construction_holding=dot(holding, u), plate=dot(plate_load, d),
                           beam=dot(beam_load, d), nodal=dot(nodal_load, d))
        result.append(dict(stage_index=row['stage_index'], progress=row['progress'],
            force_vectors_kn=dict(soil=soil.tolist(), plate_structure=plate_internal.tolist(), beam_structure=beam_internal.tolist(),
                                  structure_load=load.tolist(), mapped_structure_internal=mapped_internal.tolist(), applied=applied.tolist()),
            structural_decomposition_max=float(np.max(np.abs(structural_internal-(k@d+rest)), initial=0.)),
            global_free_residual_l2=norm(global_residual[free]), original_normalization=max(norm(applied), 1.),
            global_free_residual_relative=norm(global_residual[free])/max(norm(applied), 1.),
            explicit_free_structure_residual_l2=norm(rotation_residual),
            work_sources_kn_m=source_work, instantaneous_work_kn_m=math.fsum(source_work.values()),
            plate_load_virtual_work_max=max(map(abs, work_identity), default=0.), energy_kn_m=energy,
            unsupported_contributions=dict(joint=None, water=None, plastic_soil=None, mcc=None, idh=None, anchor=None, soil_birth_death_matset=None)))
    return result


def plate_probes(folder, coarse_folder, kind):
    mesh, case = [read(folder/(n+'.json')) for n in ('mesh', 'stages')]
    cm, cc = [read(coarse_folder/(n+'.json')) for n in ('mesh', 'stages')]
    model, coarse = System(mesh, case['beam_system']), System(cm, cc['beam_system'])
    samples = legacy['cases'].ref  # Input geometry; no nearest-neighbour mapping.
    del samples
    from plate6_reference import SAMPLES
    queries = [np.asarray(l)@coarse.positions[nodes][:3] for nodes, _, _ in coarse.plate_data for l in SAMPLES]
    owners = []
    for query in queries:
        found = []
        for pi, (nodes, _, kernel) in enumerate(model.plate_data):
            points = model.positions[nodes][:3]
            r = np.linalg.solve((points[1:, :2]-points[0, :2]).T, query[:2]-points[0, :2])
            bary = np.r_[1-r.sum(), r]
            if bary.min() >= -1e-8:
                found.append((pi, nodes, kernel, bary))
        if len(found) != 1:
            raise ValueError('fixed probe owner is not unique')
        owners.append(found[0])
    out = {}
    for key, row in common(rows(folder, kind)).items():
        s = row.get('state', row)
        bs = s['beams']
        d = np.c_[bs['node_displacement_m'], bs['node_rotation_rad']]
        qvalues, raw, reconstruction = [], [], []
        if any(bs['active_plates']):
            for pi, nodes, kernel, bary in owners:
                kin = d[nodes].ravel()-bs['plate_reference_global_kinematics'][pi]
                local = kernel.neutral.transform@kernel.map@kin
                xy = bary@kernel.neutral.xy[:3]
                direct = kernel.neutral.shear*(kernel.neutral.operators(bary)[2]@local)
                shear = direct
                if kind == 'native':
                    plate = next(v for v in row['beam_results']['plates'] if v['id'] == model.plates[pi]['id'])
                    qxy = [np.asarray(v['barycentric'])@kernel.neutral.xy[:3] for v in plate['samples']]
                    basis = np.vstack([kernel.neutral.basis(*v) for v in qxy])
                    coeff = np.linalg.lstsq(basis, np.asarray([v['transverse_shear_force_kn_m'] for v in plate['samples']]).ravel(), rcond=None)[0]
                    shear = kernel.neutral.basis(*xy)@coeff
                qvalues.extend(kernel.axes.T@np.r_[shear, 0.])
                raw.extend(kernel.axes.T@np.r_[kernel.neutral.shear*(kernel.neutral.raw_shear(*xy)@local), 0.])
                reconstruction.append(float(np.max(np.abs(shear-direct))))
        out[key] = dict(q=qvalues, raw_gradient_q=raw,
                        recovered_vs_direct_max=max(reconstruction, default=0.))
    return out


def run(formal, output):
    plan = read(formal/'plan-results.json')
    cases, detailed, probes = {}, {}, {}
    for rec in plan['cases']:
        name = rec['case']
        folder = formal/name
        cases[name] = dict(legacy_existing_gates=rec['comparison']['existing_field_gates_passed'],
                           input_only_vs_native=compare_paths(rows(folder, 'native'), rows(folder, 'independent')))
        detailed[name] = {kind: force_work_energy(folder, kind) for kind in ('native', 'independent')}
        probes[name] = {kind: plate_probes(folder, formal/'tet12-n4', kind) for kind in ('native', 'independent')}
    reference_sensitivity = {}
    for size in (12, 96):
        folder = formal/f'tet{size}-n4'
        reference_sensitivity[str(size)] = {label: compare_paths(rows(folder, label), rows(folder, 'independent')) for label in ('gl7', 'reverse-q')}
    series = {}
    work = {rec['case']: {kind: legacy['work_audit'].calculate(formal/rec['case'], kind) for kind in ('native', 'independent')} for rec in plan['cases']}
    for size in (12, 96):
        series[str(size)] = {}
        for kind in ('native', 'independent'):
            pairs = []
            for a, b in ((4, 8), (8, 16)):
                na, nb = f'tet{size}-n{a}', f'tet{size}-n{b}'
                pa, pb = probes[na][kind], probes[nb][kind]
                wa, wb = common(work[na][kind]), common(work[nb][kind])
                pairs.append(dict(steps=[a, b], state=compare_paths(rows(formal/na, kind), rows(formal/nb, kind)),
                                  q=metric([pa[k]['q'] for k in sorted(pa)], [pb[k]['q'] for k in sorted(pb)]),
                                  instantaneous_work=metric([wa[k]['global_F_dot_d_kn_m'] for k in sorted(wa)], [wb[k]['global_F_dot_d_kn_m'] for k in sorted(wb)]),
                                  cumulative_work=metric([wa[k]['cumulative_total_work_kn_m'] for k in sorted(wa)], [wb[k]['cumulative_total_work_kn_m'] for k in sorted(wb)])))
            series[str(size)][kind] = pairs
    spatial = {}
    for kind in ('native', 'independent'):
        a, b = probes['tet12-n4'][kind], probes['tet96-n4'][kind]
        spatial[kind] = metric([a[k]['q'] for k in sorted(a)], [b[k]['q'] for k in sorted(b)])
    maximum = dict(plate_virtual_work=0., membrane_strain_identity=0., membrane_energy_identity=0., structural_decomposition=0.,
                   global_free_residual_native=0., global_free_residual_independent=0., explicit_structure_residual=0., recovered_vs_direct=0.)
    source_differences = {}
    for name, kinds in detailed.items():
        for kind, values in kinds.items():
            for row in values:
                maximum['plate_virtual_work'] = max(maximum['plate_virtual_work'], row['plate_load_virtual_work_max'])
                maximum['structural_decomposition'] = max(maximum['structural_decomposition'], row['structural_decomposition_max'])
                maximum['global_free_residual_'+kind] = max(maximum['global_free_residual_'+kind], row['global_free_residual_relative'])
                maximum['explicit_structure_residual'] = max(maximum['explicit_structure_residual'], row['explicit_free_structure_residual_l2'])
                for energy in row['energy_kn_m']:
                    maximum['membrane_strain_identity'] = max(maximum['membrane_strain_identity'], energy['strain_identity_max'])
                    maximum['membrane_energy_identity'] = max(maximum['membrane_energy_identity'], abs(energy['membrane_identity_residual']))
            maximum['recovered_vs_direct'] = max(maximum['recovered_vs_direct'], max(v['recovered_vs_direct_max'] for v in probes[name][kind].values()))
        source_differences[name] = {field: metric([r['work_sources_kn_m'][field] for r in kinds['native']], [r['work_sources_kn_m'][field] for r in kinds['independent']]) for field in kinds['native'][0]['work_sources_kn_m']}
    write(formal/'force-work-energy-raw.json', detailed)
    write(formal/'work-raw.json', work)
    result = dict(schema='b02cause-measurements/1', passed=False, b02_complete=False, diagnostic_only=True,
                  numerical_acceptance_all_fields=None, cases=cases, reference_sensitivity=reference_sensitivity,
                  load_step_sensitivity=series, spatial_12_to_96=spatial, work_source_differences=source_differences, maximum=maximum,
                  restarts=plan['restarts'], total_state_pairs=sum(r['accepted_states'] for r in plan['cases']),
                  scope=read(formal/'scope/refusals.json'), mutations=read(formal/'mutations/mutations.json'),
                  residual_768='quotation only: 5e3d191; no new 768 solve or precision/stopping diagnostic',
                  unsupported_combined_equilibrium=None, single_cause_identified=False,
                  order='stage, exact progress, declared ID, C-order component; math.fsum norms/dots',
                  limitations='two spatial levels cannot give asymptotic order; no SPR replacement, extrapolation or new acceptance definition')
    write(output, result)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    run(a.run, a.output)
