"""Input-only attempt at the existing 5 kPa MCC/contact 3D path.

Combines the existing independent MCC continuum rate equations/midpoint event
integrator with the *closed sticking* branch of the existing JOINT spring law.
The branch is checked after every accepted reference step. Opening or sliding
ends this restricted reference; it is never silently replaced by a new law.
No native state, tangent, displacement, fraction or history generates a field.
This diagnostic deliberately does not certify acceptance.
"""
import argparse
import json
from pathlib import Path
import time

import numpy as np
import fem_reference as fem
import mcc_rate_fem as rate
from mcc_rate_midpoint import advance_cached
from mcc_rate_rk4 import advance_rk4
from joint_reference import System as Joint
from probe_low_preload_nc import expected_inputs, write
from verify_high_order import physical


class ClosedStickSystem(rate.System):
    def __init__(self, mesh, ground, active, case):
        super().__init__(mesh, ground, active, case)
        self.joint = Joint(mesh, case['joint_system'])
        self.live = np.array([active[p['owner_a']] and active[p['owner_b']] for p in self.joint.points])
        n = self.joint.normal
        nn = np.einsum('qa,qb->qab', n, n)
        self.spring = self.joint.kn[:, None, None]*nn + self.joint.kt[:, None, None]*(np.eye(3)-nn)
        self.kjoint = np.einsum('qag,qab,qbh,q->gh', self.joint.c, self.spring,
                               self.joint.c, self.joint.weight*self.live, optimize=True)
        transform = np.zeros((self.xyz.size, len(self.free)))
        for i, dof in enumerate(self.free):
            v = np.zeros(self.xyz.size); v[dof] = 1.
            transform[:, i] = self.lift(v)
        self.kr = transform.T@self.kjoint@transform
        self.maximum_linear_residual = 0.

    def linear_solve(self, local_stiffness, rhs):
        n = len(self.free)
        values = np.bincount(self.inverse, weights=local_stiffness.reshape(-1)[self.mask],
                             minlength=len(self.sparse_rows))
        matrix = self.kr.copy()
        matrix[self.sparse_rows, self.sparse_cols] += values
        solution = np.linalg.solve(matrix, rhs)
        error = float(np.linalg.norm(matrix@solution-rhs)/max(np.linalg.norm(rhs), 1e-30))
        self.maximum_linear_residual = max(self.maximum_linear_residual, error)
        if error > 1e-10:
            raise ValueError('independent contact/MCC linear residual exceeded existing tolerance: '+str(error))
        return solution

    def contacts(self, displacement):
        jump = np.einsum('qag,g->qa', self.joint.c, displacement)
        traction = np.einsum('qab,qb->qa', self.spring, jump)
        gap = np.einsum('qa,qa->q', jump, self.joint.normal)
        pressure = -np.einsum('qa,qa->q', traction, self.joint.normal)
        shear = traction + pressure[:, None]*self.joint.normal
        # Existing JOINT audit absolute tolerances, not new acceptance values.
        if np.any(gap[self.live] > 1e-10):
            raise ValueError('independent path leaves closed-contact branch')
        if np.any(np.linalg.norm(shear[self.live], axis=1) >
                  self.joint.mu[self.live]*pressure[self.live]+1e-7):
            raise ValueError('independent path leaves sticking Coulomb branch')
        energy = .5*np.einsum('qa,qa->q', jump, traction)
        return dict(jump_m=jump.tolist(), traction_kpa=traction.tolist(),
                    plastic_slip_m=np.zeros_like(jump).tolist(), elastic_energy_kpa_m=energy.tolist())

    def total_force(self, material, displacement):
        return self.force(material)+self.kjoint@displacement


def solve(output, steps, method='midpoint', budget_seconds=600.):
    if type(steps) is not int or steps < 2:
        raise ValueError('at least two reference steps required')
    output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    if method not in ('midpoint', 'rk4') or not np.isfinite(budget_seconds) or budget_seconds <= 0:
        raise ValueError('invalid reference method or budget')
    mesh, ground, case = expected_inputs()
    for name, value in zip(('mesh', 'ground', 'stages'), (mesh, ground, case)):
        write(output/(name+'.json'), value)
    active = np.ones(len(mesh['tetrahedra']), dtype=bool)
    state = rate.seed(mesh, ground)
    u = np.zeros(np.asarray(mesh['nodes']).size)
    frames = []; segments = []; last = None; contacts = None
    system = ClosedStickSystem(mesh, ground, active, case)
    factor = 0.
    stage_name = 'initial_equilibrated'
    maximum_residual = 0.

    def load():
        force = physical(mesh, ground, active, 0., True, joint_system=system.joint)[0].reshape(-1)
        return force+factor*fem.conditions_load(mesh, case['boundary'], active, skip_removed_surfaces=True).reshape(-1)

    def advance(before, after):
        nonlocal state, u, contacts, last, maximum_residual
        shape = state[active].shape
        def observe(progress, material, displacement):
            nonlocal contacts, last, maximum_residual
            if time.monotonic()-start > budget_seconds:
                raise ValueError('independent reference wall time budget exhausted')
            current = state.copy(); current[active] = material.reshape(shape)
            evaluated = system.contacts(displacement)
            if contacts is None: contacts = evaluated
            else:
                for key, value in evaluated.items():
                    previous = np.asarray(contacts[key]); previous[system.live] = np.asarray(value)[system.live]
                    contacts[key] = previous.tolist()
            external = before+progress*(after-before)
            residual = float(np.linalg.norm(system.reduce(system.total_force(material, displacement)-external)) /
                             max(np.linalg.norm(external), 1.))
            maximum_residual = max(maximum_residual, residual)
            last = dict(name=stage_name, progress=progress, active_elements=active.tolist(),
                        material_state=current.tolist(), displacement_m=displacement.reshape(-1, 3).tolist(),
                        contacts=json.loads(json.dumps(contacts)), relative_free_force_residual=residual)
        targets = (.875, .9375, 1.) if stage_name == 'excavate_upper' else (1.,)
        previous = 0.
        for target in targets:
            width = target-previous
            def segment_observe(progress, material, displacement):
                observe(previous+width*progress, material, displacement)
            integrator = advance_cached if method == 'midpoint' else advance_rk4
            material, u, diagnostics = integrator(system, state[active].reshape(-1, 14), u,
                                                      (after-before)*width, max(2, round(steps*width)),
                                                      accepted_callback=segment_observe)
            state[active] = material.reshape(shape)
            diagnostics['maximum_true_linear_residual'] = system.maximum_linear_residual
            segments.append(dict(name=stage_name, progress_from=previous, progress_to=target, **diagnostics))
            if target < 1.:
                write(output/(stage_name+'-progress-'+str(target)+'.json'), last)
            previous = target
        frames.append(last)
        write(output/(stage_name+'.json'), last)

    error = None
    try:
        advance(system.total_force(state.reshape(-1, 14), u), load())
        strain = np.max(np.abs(np.einsum('eqai,ei->eqa', system.b, u[system.dofs])))
        if strain > ground['initialization']['geostatic']['maximum_initial_engineering_strain']:
            raise ValueError('independent initial strain exceeded existing input limit')
        for stage in case['stages']:
            stage_name = stage['name']; before = load()
            if stage['kind'] == 'boundary':
                factor = stage['target_factor']
            else:
                regions = {r['id'] for r in mesh['regions'] if r['body'] in stage['bodies']}
                active[np.isin(mesh['region_ids'], list(regions))] = False
                system = ClosedStickSystem(mesh, ground, active, case)
                before = system.total_force(state[active].reshape(-1, 14), u)
            advance(before, load())
    except (ValueError, np.linalg.LinAlgError) as exc:
        error = str(exc)
    report = dict(schema='soil-fem-low-preload-independent-rate-attempt/1', passed=False,
                  complete=error is None, steps=steps, integration_method=method,
                  budget_seconds=budget_seconds, elapsed_seconds=time.monotonic()-start,
                  failure_reason=error, last_accepted=last, completed_frames=[f['name'] for f in frames],
                  segments=segments, maximum_relative_free_force_residual=maximum_residual,
                  last_system_maximum_true_linear_residual=system.maximum_linear_residual,
                  input_only=True, uses_native_fields_or_fractions=False,
                  method='existing compatible 3D MCC continuum rate + '+method+' events; independently assembled closed-stick JOINT spring branch',
                  limitations=['閉鎖固着枝のみ。開離・滑りは拒否', '独立参照の刻み感度・全場比較前は受入不可',
                               '接触の数値散逸・端点仕事は離散経路依存。本体の履歴監査と区別'])
    write(output/'reference.json', report)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--steps', type=int, required=True)
    parser.add_argument('--method', choices=('midpoint', 'rk4'), default='midpoint')
    parser.add_argument('--budget-seconds', type=float, default=600.)
    args = parser.parse_args()
    result = solve(args.output, args.steps, args.method, args.budget_seconds)
    print(json.dumps({k: v for k, v in result.items() if k != 'last_accepted'}, ensure_ascii=False, indent=2))
