"""Bounded research diagnostic: Tet10 u / material-split P1 volume projection.

This is NOT a production solver or an acceptance of the existing GP-volume
backend. It changes the strain space and the work-conjugate equilibrium force.
It reuses the independent analytic MCC rate equations and event integrator,
not Rust constitutive updates, histories or solved fields. No stabilization,
stress clipping, artificial stiffness or pressure floor is introduced.
"""
import os
for _name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[_name] = '1'
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import mcc_rate_fem as rate
import fem_reference as ref
import mcc_rate_midpoint
import verify_high_order
import high_order_io


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class ProjectedSystem(rate.System):
    """Eliminate M theta = D u and M p = integral(N mean(sigma)).

    Bbar = dev(B) + I/3 N M^-1 D. Equilibrium is Bbar.T sigma,
    equivalently dev(B).T sigma + D.T p; raw GP sigma need not balance
    B.T sigma. Pressure consistency is weak, explicitly not pointwise.
    """

    def __init__(self, mesh, ground, active, case=None):
        super().__init__(mesh, ground, active, case)
        if len(self.free) > 1500:
            raise ValueError('dense diagnostic is limited to 1500 free displacement DOFs')
        keys = sorted({(node, mesh['material_ids'][int(e)])
                       for e in self.element_ids for node in mesh['tetrahedra'][int(e)]})
        self.keys = keys
        lookup = {key: i for i, key in enumerate(keys)}
        self.pi = np.asarray([[lookup[(node, mesh['material_ids'][int(e)])]
                               for node in mesh['tetrahedra'][int(e)]] for e in self.element_ids])
        q = [ref.volume(mesh, int(e)) for e in self.element_ids]
        self.n = np.stack([np.column_stack((1-v['natural'].sum(axis=1), v['natural'])) for v in q])
        self.mass = np.zeros((len(keys), len(keys)))
        self.divergence = np.zeros((len(keys), self.xyz.size))
        mass_local = np.einsum('eqa,eqb,eq->eab', self.n, self.n, self.weight)
        div_local = np.einsum('eqa,eqi,eq->eai', self.n, self.b[:, :, :3].sum(axis=2), self.weight)
        for e in range(len(q)):
            self.mass[np.ix_(self.pi[e], self.pi[e])] += mass_local[e]
            self.divergence[np.ix_(self.pi[e], self.dofs[e])] += div_local[e]
        self.a_full = np.linalg.solve(self.mass, self.divergence)
        self.a = self.a_full
        if self.rotation is not None:
            self.a = np.einsum('pni,nij->pnj', self.a.reshape(len(keys), -1, 3), self.rotation).reshape(self.a.shape)
        self.a = self.a[:, self.free]
        self.bd = self.b.copy()
        self.bd[:, :, :3] -= self.b[:, :, :3].sum(axis=2)[:, :, None, :]/3
        self.rbd = self.rb.copy()
        self.rbd[:, :, :3] -= self.rb[:, :, :3].sum(axis=2)[:, :, None, :]/3
        reduced = np.full(self.xyz.size, -1)
        reduced[self.free] = np.arange(len(self.free))
        self.local_free = reduced[self.dofs]
        self.maximum_true_linear_residual = 0.
        self.checks = self.audit_operators()

    def strain(self, displacement):
        theta = self.a_full @ displacement
        volumetric = np.einsum('eqa,ea->eq', self.n, theta[self.pi])
        out = np.einsum('eqai,ei->eqa', self.bd, displacement[self.dofs])
        out[:, :, :3] += volumetric[:, :, None]/3
        return out

    def volume(self, displacement):
        return self.a_full @ displacement

    def mass_solve(self, rhs):
        return np.linalg.solve(self.mass, rhs)

    def project_volume(self, theta):
        rhs = np.zeros(len(self.keys))
        local = np.einsum('eqa,eq,eq->ea', self.n, theta, self.weight)
        np.add.at(rhs, self.pi.reshape(-1), local.reshape(-1))
        return self.mass_solve(rhs)

    def pressure(self, state):
        mean = state[:, :3].mean(axis=1).reshape(self.weight.shape)
        rhs = np.zeros(len(self.keys))
        local = np.einsum('eqa,eq,eq->ea', self.n, mean, self.weight)
        np.add.at(rhs, self.pi.reshape(-1), local.reshape(-1))
        nodal = np.linalg.solve(self.mass, rhs)
        return nodal, rhs

    def force(self, state):
        stress = state[:, :6].reshape(*self.weight.shape, 6)
        local = -np.einsum('eqai,eqa,eq->ei', self.bd, stress, self.weight, optimize=True)
        out = np.zeros(self.xyz.size)
        np.add.at(out, self.dofs.reshape(-1), local.reshape(-1))
        return out-self.a_full.T @ self.pressure(state)[1]

    def condensed_tangent(self, ep):
        ep = ep.reshape(*self.weight.shape, 6, 6)
        local = np.einsum('eqai,eqab,eqbj,eq->eij', self.rbd, ep, self.rbd, self.weight, optimize=True)
        nf, nt = len(self.free), len(self.keys)
        values = np.bincount(self.inverse, weights=local.reshape(-1)[self.mask], minlength=len(self.sparse_rows))
        dd = np.zeros((nf, nf))
        dd[self.sparse_rows, self.sparse_cols] = values
        cm = ep[:, :, :, :3].sum(axis=3)/3
        mc = ep[:, :, :3, :].sum(axis=2)/3
        mcm = cm[:, :, :3].sum(axis=2)/3
        dt_local = np.einsum('eqai,eqa,eqb,eq->eib', self.rbd, cm, self.n, self.weight, optimize=True)
        td_local = np.einsum('eqa,eqb,eqbi,eq->eai', self.n, mc, self.rbd, self.weight, optimize=True)
        tt_local = np.einsum('eqa,eq,eqb,eq->eab', self.n, mcm, self.n, self.weight, optimize=True)
        dt, td, tt = np.zeros((nf, nt)), np.zeros((nt, nf)), np.zeros((nt, nt))
        for e, indices in enumerate(self.local_free):
            used = indices >= 0
            dt[np.ix_(indices[used], self.pi[e])] += dt_local[e, used]
            td[np.ix_(self.pi[e], indices[used])] += td_local[e][:, used]
            tt[np.ix_(self.pi[e], self.pi[e])] += tt_local[e]
        return dd+dt@self.a+self.a.T@td+self.a.T@tt@self.a

    def rate(self, state, force_rate, plastic):
        self.calls += 1
        yf, d, n, dn, den = rate.terms(state, self.parameters)
        if np.any(den[plastic] <= 0):
            raise ValueError('nonpositive active MCC consistency denominator')
        ep = d.copy()
        ep[plastic] -= np.einsum('qa,qb,q->qab', dn[plastic], dn[plastic], 1/den[plastic])
        matrix = self.condensed_tangent(ep)
        rhs = self.reduce(force_rate)
        reduced = np.linalg.solve(matrix, rhs)
        error = float(np.linalg.norm(matrix@reduced-rhs)/max(np.linalg.norm(rhs), 1e-30))
        self.maximum_true_linear_residual = max(self.maximum_true_linear_residual, error)
        if error > 1e-10:
            raise ValueError('projected-volume diagnostic true linear residual exceeded tolerance')
        coordinates = np.zeros(self.xyz.size)
        coordinates[self.free] = reduced
        du = self.lift(coordinates)
        de = -self.strain(du).reshape(-1, 6)
        divisor = np.where(plastic, den, np.maximum(np.abs(den), 1.))
        potential = np.einsum('qa,qa->q', dn, de)/divisor
        flow = (potential*plastic)[:, None]*n
        dy = np.zeros_like(state)
        dy[:, :6] = np.einsum('qab,qb->qa', d, de-flow)
        dy[:, 6] = np.exp(state[:, 7])/(self.parameters[:, 1]-self.parameters[:, 2])*flow[:, :3].sum(axis=1)
        dy[:, 7] = -de[:, :3].sum(axis=1)
        dy[:, 8:] = flow
        return du, dy, potential*np.linalg.norm(n, axis=1), yf

    def audit_operators(self):
        gradient = np.array([[.1, .2, -.3], [.4, -.2, .5], [.2, -.3, .6]])
        affine = (self.xyz@gradient.T+np.array([.2, -.5, .3])).reshape(-1)
        error = float(abs(self.strain(affine)-np.einsum('eqai,ei->eqa', self.b, affine[self.dofs])).max())
        rng = np.random.default_rng(1251025)
        stress = rng.normal(size=(*self.weight.shape, 6))
        state = np.zeros((self.weight.size, 14));state[:, :6] = stress.reshape(-1, 6)
        virtual = rng.normal(size=self.xyz.size)
        work = float(np.einsum('eqa,eqa,eq->', self.strain(virtual), stress, self.weight))
        adjoint = abs(float(virtual@self.force(state))+work)/max(abs(work), 1.)
        c = np.eye(6)*2.;c[:3, :3] += 1.
        matrix = self.condensed_tangent(np.broadcast_to(c, (self.weight.size, 6, 6)))
        coordinates = np.zeros(self.xyz.size);coordinates[self.free] = rng.normal(size=len(self.free))
        direction = self.lift(coordinates)
        trial = state.copy();trial[:, :6] = np.einsum('ab,eqb->eqa', c, -self.strain(direction)).reshape(-1, 6)
        action = self.reduce(self.force(trial))
        tangent = float(np.linalg.norm(action-matrix@coordinates[self.free])/max(np.linalg.norm(action), 1.))
        assert error < 1e-11 and adjoint < 1e-11 and tangent < 1e-11
        return dict(affine_patch_maximum_strain_error=error, virtual_work_relative_error=adjoint,
                    independently_applied_tangent_relative_error=tangent)


def solve(mesh_file, ground_file, case_file, identity_file, output, fraction, steps,
          system_factory=ProjectedSystem, extra_sources=()):
    if not 0 < fraction <= 1 or steps < 2:
        raise ValueError('fraction in (0,1] and at least two reference steps required')
    mesh, ground, case, identity = map(read, (mesh_file, ground_file, case_file, identity_file))
    cfg = ground['initialization'].get('geostatic')
    if cfg is None or cfg['initial_boundary_factor'] != 0 or ground['initialization']['water_level_m'] is not None:
        raise ValueError('dry geostatic zero initial construction-boundary factor required')
    if any(case.get(k) for k in ('joint_system', 'beam_system', 'anchor_system', 'initially_inactive_bodies')):
        raise ValueError('diagnostic excludes contacts, members and initially inactive domains')
    sources = [mesh_file, ground_file, case_file, identity_file, Path(__file__), Path(rate.__file__),
               Path(ref.__file__), Path(mcc_rate_midpoint.__file__), Path(verify_high_order.__file__), Path(high_order_io.__file__)]
    sources += list(extra_sources)
    fingerprints = {str(p): sha(p) for p in sources}
    system = system_factory(mesh, ground, np.ones(len(mesh['tetrahedra']), dtype=bool), case)
    initial = rate.seed(mesh, ground, identity).reshape(-1, 14)
    start = system.force(initial)
    physical = rate.physical(mesh, ground, system.active, ground['initialization']['surface_surcharge_kpa'],
                             case.get('pressure_on_original_surface_only', False))[0].reshape(-1)
    change = fraction*(physical-start)
    final, u, diagnostics = system.advance(initial, np.zeros(system.xyz.size), change, steps,
        progress=lambda t, n, e: print(f'Projected-volume initialization {fraction*t:.9f}, steps {n}, events {e}', flush=True))
    residual = float(np.linalg.norm(system.reduce(start+change-system.force(final)))/max(np.linalg.norm(start+change), 1.))
    pressure, rhs = system.pressure(final)
    pg = np.einsum('eqa,ea->eq', system.n, pressure[system.pi])
    pm = final[:, :3].mean(axis=1).reshape(system.weight.shape)
    theta = system.volume(u)
    report = dict(schema='soil-fem-projected-volume-mcc-diagnostic/1', diagnostic_only=True,
        production_backend=False, goal_complete=False, p6_complete=False, full_initialization_accepted=False,
        load_fraction=fraction, reference_steps=steps, operator_checks=system.checks,
        volume_space='Material-split continuous P1; L2 projection of compatible volumetric strain',
        equilibrium_definition='Work-conjugate Bbar.T sigma, or dev(B).T sigma + D.T projected mean pressure',
        pressure_definition='P1 L2 projection; GP material pressure and P1 pressure are distinct fields',
        displacement_m=u.reshape(-1, 3).tolist(), material_state=final.reshape(*system.weight.shape, 14).tolist(),
        volume_dofs=[dict(node_index=k[0], material_id=k[1], tension_positive=float(v)) for k, v in zip(system.keys, theta)],
        pressure_dofs=[dict(node_index=k[0], material_id=k[1], compression_positive_kpa=float(v)) for k, v in zip(system.keys, pressure)],
        relative_work_conjugate_force_residual=residual,
        relative_compatible_stress_force_residual=float(np.linalg.norm(system.reduce(start+change-rate.System.force(system, final)))/max(np.linalg.norm(start+change), 1.)),
        maximum_gp_pressure_projection_difference_kpa=float(abs(pg-pm).max()),
        pressure_projection_relative_residual=float(np.linalg.norm(system.mass@pressure-rhs)/max(np.linalg.norm(rhs), 1.)),
        volume_projection_relative_residual=float(np.linalg.norm(system.mass@theta-system.divergence@u)/max(np.linalg.norm(system.divergence@u), 1.)),
        maximum_material_strain=float(abs(system.strain(u)).max()),
        maximum_kinematic_strain=float(abs(np.einsum('eqai,ei->eqa', system.b, u[system.dofs])).max()),
        plastic_points=int(np.count_nonzero(np.linalg.norm(final[:, 8:], axis=1)>1e-10)),
        minimum_mean_stress_kpa=float(pm.min()), maximum_normalized_yield=float(rate.terms(final, system.parameters)[0].max()),
        maximum_true_linear_residual=system.maximum_true_linear_residual, diagnostics=diagnostics,
        input_and_source_sha256=fingerprints,
        limitations='Changed diagnostic discretization, no production integration, no pointwise pressure-consistency acceptance, no excavation/restart/mesh-convergence acceptance')
    report['reference_backend'] = system_factory.__module__+'.'+system_factory.__name__
    if hasattr(system, 'solver_diagnostics'):
        report['sparse_solver_diagnostics'] = system.solver_diagnostics()
    assert residual < 1e-8 and report['maximum_material_strain'] < cfg['maximum_initial_engineering_strain']
    assert fingerprints == {str(p): sha(p) for p in sources}, 'source/input changed during diagnostic'
    with output.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, allow_nan=False)
    print(json.dumps({k: report[k] for k in ('load_fraction', 'reference_steps', 'plastic_points',
        'maximum_material_strain', 'minimum_mean_stress_kpa', 'relative_work_conjugate_force_residual',
        'relative_compatible_stress_force_residual', 'maximum_gp_pressure_projection_difference_kpa', 'operator_checks')}), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for name in ('mesh', 'ground', 'case', 'identity', 'output'):
        p.add_argument(name, type=Path)
    p.add_argument('--fraction', type=float, required=True)
    p.add_argument('--steps', type=int, default=50)
    a = p.parse_args()
    solve(*(getattr(a, k).resolve() for k in ('mesh', 'ground', 'case', 'identity', 'output')), a.fraction, a.steps)
