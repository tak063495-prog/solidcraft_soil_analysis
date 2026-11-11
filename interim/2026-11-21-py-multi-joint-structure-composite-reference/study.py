"""Small raw-history/order/load-step/GP studies. No extrapolation or acceptance."""
import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import platform
import sys
import time
import numpy as np
import scipy
from reference import CASES, dump, run
from compare import TOLERANCES, CITATION


def summarize(doc):
    frames = doc['frames']
    diag = [f['diagnostics'] for f in frames]
    final = frames[-1]
    return dict(nodes=len(doc['node_ids']), joints=doc['joint_ids'],
                contact_gp=len(final['contacts']), elements=len(final['elements']), frames=len(frames),
                maximum_node_residual_kn=max(d['max_node_residual_kn'] for d in diag),
                maximum_matrix_vs_stress_scatter_kn=max(d['matrix_vs_stress_scatter_kn'] for d in diag),
                maximum_resultant_kn=max(np.max(abs(np.array(d['resultant_kn']))) for d in diag),
                maximum_moment_kn_m=max(np.max(abs(np.array(d['external_moment_kn_m']))) for d in diag),
                maximum_contact_moment_kn_m=max(np.max(abs(np.array(d['contact_moment_kn_m']))) for d in diag),
                maximum_member_moment_kn_m=max(np.max(abs(np.array(d['member_moment_kn_m']))) for d in diag),
                maximum_work_defect_kn_m=max(abs(d['work_defect_kn_m']) for d in diag),
                final_work_defect_kn_m=final['diagnostics']['work_defect_kn_m'],
                external_work_kn_m=final['external_work_kn_m'],
                stored_energy_kn_m=sum(final[k] for k in ('bulk_energy_kn_m', 'contact_energy_kn_m', 'member_energy_kn_m')),
                friction_dissipation_kn_m=final['friction_dissipation_kn_m'],
                maximum_accumulated_slip_m=max(r['accumulated_slip_m'] for r in final['contacts']),
                final_member_forces_kn={r['id']: r['axial_force_kn'] for r in final['members']},
                states_observed=sorted({r['state'] for f in frames for r in f['contacts']}),
                # Quoted force/work checks are diagnostics, never A-01 expansion.
                quoted_force_diagnostic=max(d['max_node_residual_kn'] for d in diag) <= TOLERANCES['nodal_force_kn'],
                quoted_work_diagnostic=max(abs(d['work_defect_kn_m']) for d in diag) <= TOLERANCES['work_kn_m'],
                accepted=None)


def fields(doc, frames=None):
    output = defaultdict(list)
    states = []
    for f in doc['frames'] if frames is None else frames:
        for n in sorted(f['nodes'], key=lambda r: r['id']):
            for key in ('displacement_m', 'bulk_internal_kn', 'contact_internal_kn', 'structure_internal_kn', 'reaction_kn'):
                output[key].extend(n[key])
        for p in sorted(f['contacts'], key=lambda r: (r['joint_id'], r['gp_id'])):
            for key in ('traction_kpa', 'plastic_slip_m', 'jump_m'):
                output[key].extend(p[key])
            output['accumulated_slip_m'].append(p['accumulated_slip_m'])
            states.append((p['generation'], p['state'], p['open']))
        for m in sorted(f['members'], key=lambda r: r['id']):
            output['axial_force_kn'].append(m['axial_force_kn'])
        for key in ('external_work_kn_m', 'bulk_energy_kn_m', 'contact_energy_kn_m', 'member_energy_kn_m', 'friction_dissipation_kn_m'):
            output[key].append(f[key])
    return {k: np.array(v) for k, v in output.items()}, states


def differences(a, b):
    require_same = a.keys() == b.keys()
    if not require_same:
        raise ValueError('field identities differ')
    return {k: dict(maximum_absolute=float(np.max(abs(a[k]-b[k]), initial=0.)),
                    l2_absolute=float(np.linalg.norm(a[k]-b[k])),
                    l2_relative=float(np.linalg.norm(a[k]-b[k])/np.linalg.norm(b[k])) if np.linalg.norm(b[k]) else None,
                    tolerance=None, accepted=None) for k in a}


def step_difference(coarse, fine):
    ns = len(coarse['frames'])//3; nf = len(fine['frames'])//3
    if nf % ns:
        raise ValueError('nonmatching phase increments; no temporal interpolation')
    fine_frames = [f for f in fine['frames'] if round(f['progress']*ns) == f['progress']*ns]
    a, sa = fields(coarse); b, sb = fields(fine, fine_frames)
    return dict(common_actual_endpoints=len(fine_frames), differences=differences(a, b),
                open_state_differences=sum(x != y for x, y in zip(sa, sb)), accepted=None)


def gp_fields(doc):
    # Noncoincident GPs are NOT compared as pointwise values. Only nodal fields,
    # work, section forces and explicit JOINT resultants/integrals are aligned.
    out = defaultdict(list)
    for f in doc['frames']:
        for n in sorted(f['nodes'], key=lambda r: r['id']):
            for key in ('displacement_m', 'reaction_kn'):
                out[key].extend(n[key])
        for jid in sorted(doc['joint_ids']):
            points = [r for r in f['contacts'] if r['joint_id'] == jid]
            out['joint_resultant_kn'].extend(np.sum([r['force_kn'] for r in points], axis=0))
            out['joint_slip_integral_m3'].extend(np.sum([np.array(r['plastic_slip_m'])*r['weight_m2'] for r in points], axis=0))
        for m in sorted(f['members'], key=lambda r: r['id']):
            out['axial_force_kn'].append(m['axial_force_kn'])
        out['friction_dissipation_kn_m'].append(f['friction_dissipation_kn_m'])
    return {k: np.array(v) for k, v in out.items()}


def run_study(out, steps=(8, 16, 32), gp_orders=(1, 2, 3), cases=CASES):
    out = Path(out); out.mkdir(parents=True, exist_ok=False)
    report = dict(schema='dry-multi-joint-study/1', complete=False, accepted=None,
                  A_01='未達', NC_01='未達・対象外', P6='未達・対象外',
                  tolerance_source=CITATION, quoted_diagnostic_tolerances=TOLERANCES,
                  new_acceptance_tolerances='未定義', cases={}, raw_sha256={},
                  environment=dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__, platform=platform.platform()))
    started = time.monotonic()
    try:
        for name in cases:
            cache = {}
            def get(ns, ng, variant='identity'):
                key = ns, ng, variant
                if key not in cache:
                    doc = run(name, ns, ng, variant=variant)
                    filename = f'{name}-steps{ns}-gp{ng}-{variant}.json'
                    dump(out/filename, doc)
                    report['raw_sha256'][filename] = hashlib.sha256((out/filename).read_bytes()).hexdigest()
                    cache[key] = doc
                    print(json.dumps(dict(case=name, steps=ns, gp=ng, variant=variant, seconds=doc['elapsed_seconds'])), flush=True)
                return cache[key]
            basesteps = steps[-1]
            base = get(basesteps, 2)
            entry = dict(reference=summarize(base), order={}, load_steps={}, gp={}, step_pairs=[], gp_pairs=[])
            a, sa = fields(base)
            for variant in ('contacts', 'joints', 'nodes'):
                doc = get(basesteps, 2, variant)
                b, sb = fields(doc)
                entry['order'][variant] = dict(differences=differences(a, b),
                                               open_state_differences=sum(x != y for x, y in zip(sa, sb)),
                                               discrete_history_equal=sa == sb, accepted=None)
            for ns in steps:
                entry['load_steps'][str(ns)] = summarize(get(ns, 2))
            for ns, nf in zip(steps, steps[1:]):
                entry['step_pairs'].append(dict(coarse=ns, fine=nf, **step_difference(get(ns, 2), get(nf, 2))))
            for ng in gp_orders:
                entry['gp'][str(ng)] = summarize(get(basesteps, ng))
            for ng, nf in zip(gp_orders, gp_orders[1:]):
                entry['gp_pairs'].append(dict(coarse=ng, fine=nf,
                                               differences=differences(gp_fields(get(basesteps, ng)), gp_fields(get(basesteps, nf))),
                                               pointwise_convergence='未確認; GP locations differ', accepted=None))
            report['cases'][name] = entry
            dump(out/'study.json', report)
        report['complete'] = True
    finally:
        report['elapsed_seconds'] = time.monotonic()-started
        dump(out/'study.json', report)
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--steps', type=int, nargs='+', default=[8, 16, 32])
    p.add_argument('--gp-orders', type=int, nargs='+', default=[1, 2, 3])
    p.add_argument('--cases', choices=CASES, nargs='+', default=list(CASES))
    a = p.parse_args()
    if sorted(set(a.steps)) != a.steps or len(a.steps) < 3 or sorted(set(a.gp_orders)) != a.gp_orders or len(a.gp_orders) < 3:
        p.error('at least three increasing levels; no silently reduced series')
    run_study(a.output, tuple(a.steps), tuple(a.gp_orders), tuple(a.cases))


if __name__ == '__main__':
    main()
