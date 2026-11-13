"""Tensor/divergence, Coulomb equations and work audit of raw reference histories.

Uses geometry only from make_model; never calls run, assemble or contact_return.
The physical external-work defect stays uncorrected. The separate discrete
Coulomb identity is an algebra diagnostic, not a corrected reference value.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from reference import make_model, dump
from compare import validate, TOLERANCES, CITATION


def audit(doc):
    validate(doc)
    meta = doc['metadata']; dim = meta['dimension']
    model = make_model(meta['case'], meta['gp_order'], meta['mesh_per_m'])
    coords = {nid: model.x[i] for i, nid in enumerate(model.node_ids)}
    pmap = {(p['joint_id'], p['gp_id']): p for p in model.contacts}
    old = {key: dict(slip=np.zeros(dim), jump=np.zeros(dim), q=np.zeros(dim), energy=0., accumulated=0.) for key in pmap}
    maxima = dict(bulk_tensor_scatter_kn=0., contact_scatter_kn=0., section_force_kn=0.,
                  contact_law_kpa=0., flow_slip_m=0., friction_work_identity_kn_m=0.,
                  closed_discrete_work_identity_kn_m=0., member_work_identity_kn_m=0.,
                  global_work_accounting_kn_m=0., raw_physical_work_defect_kn_m=0.,
                  plastic_update_m=0., accumulated_slip_m=0.)
    applicable = 0; inapplicable = 0; state_mismatches = 0; frames = []
    young, nu = meta['E_kpa'], meta['nu']
    mu_el = young/(2*(1+nu)); lam = young*nu/((1+nu)*(1-2*nu))
    kn, kt, friction = meta['kn_kpa_per_m'], meta['kt_kpa_per_m'], meta['mu']
    for f in doc['frames']:
        nodes = {n['id']: n for n in f['nodes']}
        bulk = {n: np.zeros(dim) for n in nodes}; contact = {n: np.zeros(dim) for n in nodes}
        displacement = {n: np.array(row['displacement_m']) for n, row in nodes.items()}
        for cell in doc['bulk_connectivity']:
            x = np.array([coords[n] for n in cell]); u = np.array([displacement[n] for n in cell])
            edge = x[1:]-x[0]
            # grad_u = (du/dnatural) * (dx/dnatural)^-1, tensor form.
            gu = np.linalg.solve(edge, u[1:]-u[0]).T
            eps = (gu+gu.T)/2
            stress = 2*mu_el*eps+lam*np.trace(eps)*np.eye(dim)
            gradients = np.linalg.inv(np.column_stack((np.ones(dim+1), x)))[1:, :].T
            volume = abs(np.linalg.det(edge))/(2 if dim == 2 else 6)
            for n, grad in zip(cell, gradients):
                bulk[n] += volume*stress@grad
        friction_identity = 0.; discrete_identity = 0.; frame_defects = []
        for r in f['contacts']:
            key = (r['joint_id'], r['gp_id']); p = pmap[key]; prev = old[key]
            normal = np.array(r['normal']); jump = np.array(r['jump_m'])
            q = np.array(r['traction_kpa']); slip = np.array(r['plastic_slip_m'])
            gap = normal@jump; jt = jump-gap*normal
            if gap >= 0:
                expected = np.zeros(dim)
                expected_slip = prev['slip']; expected_state = 'open'
            else:
                pn = -kn*gap; trial = kt*(jt-prev['slip']); norm = np.linalg.norm(trial)
                expected_t = trial if norm <= friction*pn else friction*pn*trial/norm
                expected = -pn*normal+expected_t
                expected_slip = prev['slip'] if norm <= friction*pn else jt-expected_t/kt
                expected_state = 'stick' if norm <= friction*pn else 'slip'
                qt = q-(q@normal)*normal
                ds = slip-prev['slip']
                maxima['flow_slip_m'] = max(maxima['flow_slip_m'], float(np.linalg.norm(jt-slip-qt/kt)))
                friction_identity += r['weight_m2']*(qt@ds-friction*pn*np.linalg.norm(ds))
            maxima['contact_law_kpa'] = max(maxima['contact_law_kpa'], float(np.max(abs(q-expected))))
            maxima['plastic_update_m'] = max(maxima['plastic_update_m'], float(np.max(abs(slip-expected_slip))))
            maxima['accumulated_slip_m'] = max(maxima['accumulated_slip_m'],
                abs(r['accumulated_slip_m']-prev['accumulated']-float(np.linalg.norm(slip-prev['slip']))))
            state_mismatches += int(r['state'] != expected_state)
            # Reassemble emitted tractions, preserving the physical side sign.
            for side, sign in [('minus', -1), ('plus', 1)]:
                for i, sh in zip(p[side], p['shape']):
                    contact[model.node_ids[i]] += sign*r['weight_m2']*sh*q
            # For closed endpoints with their elastic spring relation, the
            # trapezoid defect equals -1/2 (q_t,new-q_t,old).delta_plastic_slip.
            # Retain both signed defects. Never add this term to exported work.
            old_gap = normal@prev['jump']
            if gap <= 0 and old_gap <= 0:
                qt = q-(q@normal)*normal
                qt_old = prev['q']-(prev['q']@normal)*normal
                algebra = r['work_defect_increment_kn_m'] + .5*r['weight_m2']*(qt-qt_old)@(slip-prev['slip'])
                discrete_identity += algebra; applicable += 1
            else:
                inapplicable += 1
            frame_defects.append(r['work_defect_increment_kn_m'])
            old[key] = dict(slip=slip, jump=jump, q=q, energy=r['elastic_energy_kn_m'], accumulated=r['accumulated_slip_m'])
        for nid, row in nodes.items():
            maxima['bulk_tensor_scatter_kn'] = max(maxima['bulk_tensor_scatter_kn'], float(np.max(abs(bulk[nid]-row['bulk_internal_kn']))))
            maxima['contact_scatter_kn'] = max(maxima['contact_scatter_kn'], float(np.max(abs(contact[nid]-row['contact_internal_kn']))))
        for r in f['members']:
            a, b = r['node_ids']; direction = coords[b]-coords[a]; length = np.linalg.norm(direction); direction /= length
            extension = direction@(displacement[b]-displacement[a])
            ea = 600. if r['kind'] == 'beam_axial' else 300.
            axial = ea/length*extension
            maxima['section_force_kn'] = max(maxima['section_force_kn'], abs(axial-r['axial_force_kn']))
            maxima['member_work_identity_kn_m'] = max(maxima['member_work_identity_kn_m'], abs(.5*axial*extension-r['elastic_energy_kn_m']))
        energy = f['bulk_energy_kn_m']+f['contact_energy_kn_m']+f['member_energy_kn_m']
        defect = f['external_work_kn_m']-energy-f['friction_dissipation_kn_m']
        maxima['global_work_accounting_kn_m'] = max(maxima['global_work_accounting_kn_m'], abs(defect-f['diagnostics']['work_defect_kn_m']))
        maxima['raw_physical_work_defect_kn_m'] = max(maxima['raw_physical_work_defect_kn_m'], abs(defect))
        maxima['friction_work_identity_kn_m'] = max(maxima['friction_work_identity_kn_m'], abs(float(friction_identity)))
        maxima['closed_discrete_work_identity_kn_m'] = max(maxima['closed_discrete_work_identity_kn_m'], abs(float(discrete_identity)))
        frames.append(dict(step=f['step'], friction_work_identity_kn_m=float(friction_identity),
                           closed_discrete_identity_kn_m=float(discrete_identity), raw_work_defect_kn_m=defect))
    return dict(schema='dry-multi-joint-independent-audit/1', complete=True, accepted=None,
                case=meta['case'], maxima=maxima, closed_identity_evaluations=applicable,
                open_endpoint_identity_not_applied=inapplicable, state_mismatches=state_mismatches, frames=frames,
                tolerance_source=CITATION, quoted_diagnostic_tolerances=TOLERANCES,
                slip_moment_order_tolerance='未定義', physical_work_corrected=False)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    dump(a.output, audit(json.loads(a.input.read_text(encoding='utf-8'))))


if __name__ == '__main__':
    main()
