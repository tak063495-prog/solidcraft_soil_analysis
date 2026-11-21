"""Incremental work on the dry, elastic-bulk, two-interface test case."""

import numpy as np


def interval(model, basis, parent, candidate, old_load, load):
    sys = model.system
    u0, t0 = parent["u"], parent["traction"]
    u, t = candidate["u"], candidate["traction"]
    du = u - u0
    bulk0 = model.k @ u0 + model.seed
    bulk = model.k @ u + model.seed

    def reaction(b, traction, external):
        total = b + np.einsum("qag,qa,q->g", sys.c, traction, sys.weight) - external
        return total - basis @ (basis.T @ total)

    r0 = reaction(bulk0, t0, old_load)
    r = reaction(bulk, t, load)
    external = 0.5 * (old_load + load) @ du
    supports = 0.5 * (r0 + r) @ du
    elastic = 0.5 * (bulk0 + bulk) @ du

    def contact_energy(traction):
        qn = np.einsum("qa,qa->q", traction, sys.normal)
        qt = traction - qn[:, None] * sys.normal
        return qn**2 / (2 * sys.kn) + np.sum(qt**2, axis=1) / (2 * sys.kt)

    ce = sys.weight @ (contact_energy(t) - contact_energy(t0))
    friction, numerical, endpoint = (
        sys.weight @ (candidate[k] - parent[k]) for k in ("pd", "nd", "work")
    )
    return dict(
        external_trapezoid=external,
        support_trapezoid=supports,
        bulk_elastic=elastic,
        contact_elastic=ce,
        bulk_plastic=0.0,
        friction=friction,
        contact_numerical=numerical,
        contact_endpoint=endpoint,
        artificial=0.0,
        contact_identity_residual=endpoint - ce - friction - numerical,
        global_residual=external + supports - elastic - ce - friction - numerical,
    )


def native_state(parent, order):
    pts = [parent["joints"]["points"][i] for i in order]
    result = {
        key: np.array([p[name] for p in pts])
        for key, name in [
            ("traction", "traction_kpa"),
            ("pd", "plastic_dissipation_kpa_m"),
            ("nd", "numerical_dissipation_kpa_m"),
            ("work", "endpoint_work_kpa_m"),
        ]
    }
    return dict(u=np.asarray(parent["displacement_m"]).ravel(), **result)
