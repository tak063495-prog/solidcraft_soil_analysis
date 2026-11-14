"""Actual perturbed residual-domain lambda=0 Newton solves, holding unchanged."""
import argparse
import json
from pathlib import Path
import sys
import time
from types import SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [str(HERE), str(ROOT/'tools'), str(ROOT/'interim/2026-11-10-nc-mean-stress-release')]
from bounded_study import Model, solve
from diagnose import write


def run(evidence, observer, out):
    out.mkdir(parents=True, exist_ok=False)
    binary = observer/'target/release/examples'/('material_diagnostic.exe' if sys.platform == 'win32' else 'material_diagnostic')
    records = {}
    for kind in ('native', 'independent'):
        folder = out/kind
        folder.mkdir()
        args = SimpleNamespace(map=kind, evidence=evidence/'prepared/origin', binary=binary,
                               reference=evidence/'prepared/reference400', order='complement-first', material_steps=4)
        m = Model(args, folder)
        direction = np.sin(np.arange(m.R.shape[1])+1.)
        initial = 1e-4*direction/np.linalg.norm(direction)
        old, du, events, failure = solve(m, m.old, 0., folder, 'lambda-zero', time.monotonic()+60, initial)
        result = dict(passed=False, initial_guess_m=initial.tolist(), events=events, stop_reason=failure,
                      holding_load_unchanged=True, not_a_zero_displacement_identity_check=True)
        if failure is None:
            force = -np.einsum('qai,qa,q->i', m.B, m.stress(old), m.w)
            residual = float(np.linalg.norm(m.R.T@m.load(0.)-force)/max(np.linalg.norm(m.load(0.)), 1.))
            result.update(relative_force_residual=residual, states=m.serial(old),
                          displacement_m=(m.u+(m.R@du).reshape(m.u.shape)).tolist(),
                          contacts=m.contacts, removed_material_frozen=m.frozen(old))
            if residual >= 1e-8 or not m.frozen(old):
                raise ValueError('lambda-zero equilibrium/history audit failed')
        else:
            raise ValueError('lambda-zero solve incomplete: '+failure)
        write(folder/'result.json', result)
        records[kind] = {k: v for k, v in result.items() if k not in ('states', 'displacement_m', 'contacts', 'initial_guess_m')}
    write(out/'measurements.json', dict(passed=False, residual_domain_lambda_zero_verified=True, maps=records,
                                      full_removal_accepted=False))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('evidence', 'observer', 'output'):
        p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args()
    run(a.evidence.resolve(), a.observer.resolve(), a.output.resolve())
