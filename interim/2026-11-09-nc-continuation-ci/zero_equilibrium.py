"""Solve the original remaining-domain lambda=0 load from a nonzero perturbation.

Diagnostic only. Holding loads are still those defined by the pre-cut state.
This is neither full release nor an independent proof of the holding-load model.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent/'2026-11-09-nc-beyond-fold'))
from continue_path import Adapter, digest, read, write, append


def run(args):
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    m = Adapter(args.map, args.evidence, args.binary, args.reference, args.arc, out, 4)
    if args.map == 'native':
        source = m.m.checkpoint['run']['continuation']['origin']
        old = source['material_states']
        u = np.array(source['displacement_m'])
        contacts = source['joints']
    else:
        source = read(args.reference/'unload.json')
        old = np.array(source['material_state'])
        u = np.array(source['displacement_m'])
        contacts = source['contacts']
    if contacts != m.contacts:
        raise ValueError('pre-cut/removed contact history differs')
    target = m.m.load(0.)
    # Existing 1e-4 m root-search amplitude, no physical input/load adjustment.
    direction = m.t[:-1]/np.linalg.norm(m.t[:-1])
    du = 1e-4*direction
    write(out/'controls.json', dict(map=args.map, lambda_fixed=0., initial_amplitude_m=1e-4,
        initial_reduced_increment_m=du.tolist(), material_steps=4, force_tolerance=1e-8,
        true_linear_tolerance=1e-10, correction_tolerance_m=1e-11, maximum_newton_iterations=20,
        minimum_correctors=1, source_sha256=digest(source), model_sha256=digest(m.source),
        driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        source_uses_native_fields=args.map == 'native', passed=False,
        note='pre-cut-derived holding target retained; nonzero numerical initial guess; not full release'))
    write(out/'initial.json', dict(states=m.serial(old), displacement_m=u.tolist(), contacts=contacts))
    solved = False
    stop = None
    try:
        for iteration in range(21):
            f, K, new = m.evaluate(du, old, f'lambda-zero:{iteration}')
            rhs = m.R.T@target-f
            rel = float(np.linalg.norm(rhs)/max(np.linalg.norm(target), 1.))
            dx = np.linalg.solve(K, rhs)
            lr = float(np.linalg.norm(K@dx-rhs)/max(np.linalg.norm(rhs), 1e-30))
            append(out/'linear-systems.jsonl', dict(iteration=iteration, matrix=K.tolist(),
                rhs=rhs.tolist(), solution=dx.tolist(), true_relative_residual=lr))
            append(out/'events.jsonl', dict(iteration=iteration, reduced_increment_m=du.tolist(),
                relative_force_residual=rel, correction_norm_m=float(np.linalg.norm(dx))))
            if lr > 1e-10:
                raise ValueError('true linear residual gate; inconclusive')
            if iteration >= 1 and rel < 1e-8 and np.linalg.norm(dx) < 1e-11:
                write(out/'solution.json', dict(progress=0., states=m.serial(new), contacts=contacts,
                    displacement_m=(u+(m.R@du).reshape(u.shape)).tolist(),
                    reduced_increment_m=du.tolist(), relative_force_residual=rel))
                solved = True
                stop = 'lambda=0 perturbed equilibrium solved; full removal inconclusive'
                break
            if iteration == 20:
                raise ValueError('Newton limit; inconclusive')
            du += dx
    except (ValueError, np.linalg.LinAlgError, subprocess.SubprocessError) as exc:
        stop = str(exc)
    report = dict(passed=False, history_converged=False, static_limit_confirmed=False,
        diagnostic_equilibrium_solved=solved, lambda_fixed=0., stop_reason=stop,
        full_removal_reached=False, full_release_freeze_verified=False, production_adopted=False)
    write(out/'report.json', report)
    if args.map == 'native': write(out/'material-processes.json', m.m.processes)
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--map', choices=['native', 'independent'], required=True)
    for name in ['evidence', 'binary', 'reference', 'arc']:
        p.add_argument('--'+name, type=Path)
    p.add_argument('--output', type=Path, required=True)
    print(json.dumps(run(p.parse_args()), indent=2))
