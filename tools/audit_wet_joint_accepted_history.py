"""Replay existing wet JOINT histories; consistency does not certify convergence.

Only the existing small CONTINUOUS boundary-cycle case and its lossless saved
reference fields are supported. This is not a new contact specification, an
independent native solution, or a full-history acceptance threshold.
"""
import argparse
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
HISTORY = ROOT/'interim/2026-11-08-wet-joint-finer-full-history'
sys.path.insert(0, str(HISTORY))
spec = importlib.util.spec_from_file_location('wet_accepted_history_codec', HISTORY/'study.py')
history = importlib.util.module_from_spec(spec)
spec.loader.exec_module(history)
from joint_water_reference import Water
from verify_joint_coupling import contact
from wet_joint_cycle_reference import contact_work


def audit(folder, result, values):
    if values != history.inputs():
        raise ValueError('existing small wet JOINT cycle only')
    existing = history.check_history(folder, result, values)
    system = Water(*values).contact
    old_jump = np.zeros((32, 3))
    old_q = np.zeros_like(old_jump)
    plastic = np.zeros_like(old_jump)
    accumulated = np.zeros(32)
    pd = np.zeros(32)
    nd = np.zeros(32)
    work = np.zeros(32)
    maximum = {key: 0. for key in ('traction_kpa', 'plastic_slip_m', 'accumulated_m',
                                  'plastic_dissipation_kpa_m', 'numerical_dissipation_kpa_m', 'endpoint_work_kpa_m')}
    count = 0
    for row, fields in history.accepted_states(folder):
        jump = fields['jump']
        q, flow, length = contact(system, old_jump, old_q, jump, result['contact_steps'])
        pdi, ndi, wi = contact_work(system, old_jump, old_q, jump, result['contact_steps'])
        plastic += flow[:, 0]
        accumulated += length[:, 0]
        pd += pdi
        nd += ndi
        work += wi
        pairs = [('traction_kpa', fields['traction'], q[:, 0], 1e-7, 1e-9),
                 ('plastic_slip_m', fields['plastic_slip'], plastic, 1e-10, 0.),
                 ('accumulated_m', fields['accumulated'], accumulated, 1e-10, 0.),
                 ('plastic_dissipation_kpa_m', fields['plastic_dissipation'], pd, 1e-9, 1e-8),
                 ('numerical_dissipation_kpa_m', fields['numerical_dissipation'], nd, 1e-9, 1e-8),
                 ('endpoint_work_kpa_m', fields['endpoint_work'], work, 1e-9, 1e-8)]
        for key, observed, expected, atol, rtol in pairs:
            maximum[key] = max(maximum[key], float(np.max(np.abs(observed-expected))))
            np.testing.assert_allclose(observed, expected, atol=atol, rtol=rtol,
                                       err_msg=f'accepted history propagation {key} at {row["stage"]}:{row["progress"]}')
        old_jump = jump.copy()
        old_q = q[:, 0].copy()
        count += 1
    return dict(passed=False, history_consistent=True, accepted_states=count, contact_points=32,
                old_field_checks=existing, maximum_absolute_differences=maximum,
                independent_native_solution=False, full_history_convergence_accepted=False,
                p3b_p4_complete=False, structure_water_joint_complete=False, p6_complete=False,
                scope='replay of existing reference contact and work history; no new physical acceptance')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--history', type=Path, required=True)
    p.add_argument('--result', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if args.output.exists():
        p.error('new output required')
    report = audit(args.history, json.loads(args.result.read_text()), history.inputs())
    args.output.write_text(json.dumps(report, indent=2)+'\n')
