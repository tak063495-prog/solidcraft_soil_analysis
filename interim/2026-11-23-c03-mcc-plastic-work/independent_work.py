"""Accepted-substep work of the independent backward-Euler MCC reference."""
import math

import numpy as np

from c03diag_mcc_point import step


def components(stress, increment):
    """Voigt engineering strains; epsilon_s is the equivalent deviatoric increment."""
    sigma = np.asarray(stress)
    dep = np.asarray(increment)
    p = sigma[:, :3].mean(axis=1)
    ev = dep[:, :3].sum(axis=1)
    dev = sigma[:, :3] - p[:, None]
    shear_strain = dep[:, :3] - ev[:, None] / 3
    q = np.sqrt(1.5 * ((dev**2).sum(axis=1) + 2 * (sigma[:, 3:]**2).sum(axis=1)))
    es = np.sqrt((2/3) * ((shear_strain**2).sum(axis=1) + .5 * (dep[:, 3:]**2).sum(axis=1)))
    return np.column_stack(((sigma*dep).sum(axis=1), p*ev, q*es))


class AcceptedWork:
    def __init__(self, model):
        self.model = model
        self.states = model.seeds
        self.previous_u = np.zeros(model.nd)

    def observe(self, row, fields):
        u = fields["u"].ravel()
        increments = -self.model.strains(u - self.previous_u)
        sigma, plastic, weights = [], [], []
        states = []
        for i, (mat, state, de, count, weight) in enumerate(zip(
                self.model.mats, self.states, increments, row["material_partitions"],
                self.model.weights, strict=True)):
            for _ in range(count):
                next_state = step(mat, state, de/count)
                sigma.append(next_state["effective_stress_compression_positive_kpa"])
                plastic.append(np.subtract(next_state["plastic_strain_compression_positive"],
                                          state["plastic_strain_compression_positive"]))
                weights.append(weight)
                state = next_state
            for key, field in (("effective_stress_compression_positive_kpa", "stress"),
                               ("plastic_strain_compression_positive", "plastic_strain"),
                               ("preconsolidation_pressure_kpa", "pc"), ("void_ratio", "void")):
                if not np.array_equal(state[key], fields[field][i]):
                    raise ValueError(f"independent substeps differ from accepted {field} at GP {i}")
            states.append(state)
        self.previous_u = u.copy()
        self.states = states
        terms = components(sigma, plastic)
        local = terms[:, 0] - (terms[:, 1] + terms[:, 2])
        totals = [math.fsum(column) for column in (np.asarray(weights)[:, None]*terms).T]
        return dict(stage=row["stage"], progress=row["progress"], work_kn_m=totals,
                    substeps=len(terms), plastic_substeps=int(np.any(np.array(plastic) != 0, axis=1).sum()),
                    max_local_residual_kpa=float(np.abs(local).max()),
                    sum_abs_residual_kn_m=math.fsum(np.abs(local)*weights))
