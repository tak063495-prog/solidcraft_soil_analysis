"""Moving-tangent-plane continuation of a single stage in declared scaled coordinates."""

from dataclasses import dataclass
import hashlib
import json

import numpy as np
from scipy.linalg import ldl
from scipy.optimize import root


def serial(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def identity(value):
    return hashlib.sha256(
        json.dumps(
            value,
            default=serial,
            allow_nan=False,
            sort_keys=True,
            separators=(",", ":"),
        ).encode()
    ).hexdigest()


class TrialRejected(ArithmeticError):
    """A numerical/domain rejection from a supported model, distinct from an IO failure."""


@dataclass
class Evaluation:
    residual: np.ndarray
    k: np.ndarray
    q: np.ndarray
    candidate: dict
    denominator: float
    eligible: bool = True
    event: tuple = ()
    solver_k: np.ndarray | None = None
    details: dict | None = None


@dataclass(frozen=True)
class Controls:
    ds: float
    steps: int
    force_tol: float
    linear_tol: float
    iterations: int
    arc_tol: float = 1e-12
    event_bisections: int = 8


def roundoff(matrix):
    return np.finfo(float).eps * max(matrix.shape) * np.linalg.norm(matrix, 2)


def spectrum(k, q, residual):
    """Use the physical work-conjugate coordinates supplied by the caller."""
    k, q, residual = map(np.asarray, (k, q, residual))
    if not np.isfinite(k).all() or k.ndim != 2 or k.shape[0] != k.shape[1]:
        raise ValueError("finite square physical tangent required")
    u, singular, vt = np.linalg.svd(k)
    tol = roundoff(k)
    weak = singular <= tol
    # Report the weakest pair too, even when the matrix is regular.
    selected = weak.copy()
    selected[-1] = True
    left, right = u[:, selected], vt[selected].T
    err = np.linalg.norm(k - k.T, 2)
    result = dict(
        asymmetry_norm=err,
        rank_roundoff=tol,
        nullity=int(sum(weak)),
        singular_values=singular,
        left=left,
        right=right,
        left_weak_projector=u[:, weak] @ u[:, weak].T,
        right_weak_projector=vt[weak].T @ vt[weak],
        left_q=left.T @ q,
        left_residual=left.T @ residual,
        null_load_norm=np.linalg.norm(u[:, weak].T @ q),
        null_residual_norm=np.linalg.norm(u[:, weak].T @ residual),
        right_pair_residual=np.linalg.norm(
            k @ right - left * singular[selected], axis=0
        ),
        left_pair_residual=np.linalg.norm(
            k.T @ left - right * singular[selected], axis=0
        ),
    )
    if err <= tol:
        eigen, modes = np.linalg.eigh(k)
        _, block, perm = ldl(k, hermitian=True)
        signs = np.linalg.eigvalsh(block)
        # Congruence preserves exact inertia. Near a singularity the two
        # numerical estimates remain separately visible.
        result.update(
            eigenvalues=eigen,
            minimum_mode=modes[:, 0],
            negative_eigen=int(sum(eigen < -tol)),
            negative_ldlt=int(sum(signs < -roundoff(block))),
            permutation=perm,
        )
    return result


def direction(ev, previous):
    jac = np.column_stack((ev.k, ev.q))
    _, s, vt = np.linalg.svd(jac, full_matrices=True)
    if s[-1] <= roundoff(jac):
        raise ArithmeticError("non-unique branch tangent")
    tangent = vt[-1]
    return tangent if tangent @ previous >= 0 else -tangent


def trace(
    evaluate,
    checkpoint,
    z0,
    forward,
    controls,
    *,
    local_control=None,
    stop=None,
    restart=None,
):
    """evaluate(parent,z) must integrate all history from parent on every call."""
    z = np.array(z0, dtype=float)
    parent = checkpoint
    point = evaluate(parent, z)
    rows, critical = [], {}
    accepted = [dict(z=z.copy(), state=parent, event=point.event, s=0.0)]
    old = np.asarray(forward, dtype=float)
    if (
        not point.eligible
        or np.linalg.norm(point.residual) > controls.force_tol * point.denominator
    ):
        raise ValueError(
            "starting checkpoint is not an accepted equilibrium on this stage"
        )
    incoming = point
    if restart is not None:
        payload = restart["payload"]
        if (
            identity(payload) != restart["sha256"]
            or payload["parent"] != identity(parent)
            or not np.array_equal(payload["z"], z)
        ):
            raise ValueError("restart checkpoint or predictor digest differs")
        old = np.asarray(payload["previous_direction"])
        incoming = Evaluation(
            point.residual,
            np.asarray(payload["k"]),
            np.asarray(payload["q"]),
            parent,
            point.denominator,
            event=tuple(payload["event"]),
        )
    reason = "step budget"
    s = 0.0
    for _ in range(controls.steps):
        parent_id = identity(parent)
        try:
            tangent = direction(incoming, old)
        except ArithmeticError as e:
            reason = str(e)
            critical["stop"] = dict(parent=parent, evaluation=vars(point), z=z)
            break
        plane = (
            tangent if local_control is None else np.asarray(local_control, dtype=float)
        )
        if abs(plane @ tangent) <= np.finfo(float).eps:
            reason = "nonmonotone control coordinate"
            break
        ds = controls.ds
        for cut in range(controls.event_bisections + 1):
            origin = z.copy()
            predicted = z + tangent * ds
            trial = predicted.copy()
            iterations = []
            measured_linear = None
            accepted_trial = False
            for it in range(controls.iterations + 1):
                try:
                    ev = evaluate(parent, trial)
                except TrialRejected as e:
                    cause = str(e)
                    ev = None
                    absolute = None
                    arc = float(plane @ (trial - predicted))
                    iterations.append(dict(z=trial.copy(), rejection=cause))
                    break
                arc = float(plane @ (trial - predicted))
                absolute = float(np.linalg.norm(ev.residual))
                jac = np.column_stack(
                    (ev.solver_k if ev.solver_k is not None else ev.k, ev.q)
                )
                bordered = np.vstack((jac, plane))
                rhs = -np.r_[ev.residual, arc]
                iterations.append(
                    dict(z=trial.copy(), matrix=bordered, rhs=rhs, evaluation=vars(ev))
                )
                if (
                    absolute <= controls.force_tol * ev.denominator
                    and abs(arc) <= controls.arc_tol
                ):
                    accepted_trial = ev.eligible
                    cause = (
                        "equilibrium"
                        if accepted_trial
                        else "support or constraint check"
                    )
                    break
                if it == controls.iterations:
                    cause = "Newton iteration limit"
                    break
                try:
                    delta = np.linalg.solve(bordered, rhs)
                except np.linalg.LinAlgError:
                    cause = "singular corrector"
                    break
                measured_linear = np.linalg.norm(
                    bordered @ delta - rhs
                ) / np.linalg.norm(rhs)
                iterations[-1]["true_linear_residual"] = measured_linear
                if (
                    not np.isfinite(measured_linear)
                    or measured_linear > controls.linear_tol
                ):
                    cause = "true linear residual"
                    break
                trial += delta
            switch = ev is not None and ev.event != incoming.event
            if accepted_trial and switch and cut < controls.event_bisections:
                cause, accepted_trial = "localize active-set change", False
            # No callback may update the accepted history on a rejected trial.
            if identity(parent) != parent_id:
                raise ValueError("accepted history changed during correction")
            index = len(rows) + 1
            row = dict(
                id=index,
                parent=parent_id,
                accepted=accepted_trial,
                reason=cause,
                ds=ds,
                s_parent=s,
                lambda_parent=z[-1],
                delta_lambda=trial[-1] - z[-1],
                z=trial.copy(),
                direction=tangent.copy(),
                iterations=it,
                event=ev.event if ev else None,
                absolute_force=absolute,
                denominator=ev.denominator if ev else None,
                relative_force=absolute / ev.denominator if ev else None,
                true_linear_residual=measured_linear,
                plane_residual=arc,
            )
            if ev and ev.details:
                row.update(ev.details.get("history", {}))
            row["trial_state_id"] = identity(ev.candidate) if ev else None
            row["result_checkpoint"] = (
                row["trial_state_id"] if accepted_trial else parent_id
            )
            rows.append(row)
            record = dict(
                origin=origin, parent=parent, predicted=predicted, iterations=iterations
            )
            # The two neighboring accepted points are included with every
            # critical interval. Regular interior iterations can be discarded.
            crossing = old[-1] * tangent[-1] < 0
            if (
                switch
                or crossing
                or not accepted_trial
                or index == 1
                or "next" in critical
            ):
                critical.pop("next", None)
                critical[index] = record
                if accepted:
                    critical[f"before-{index}"] = accepted[-1]
                if switch or crossing:
                    critical["next"] = True
            if accepted_trial:
                old, z, parent = tangent, trial.copy(), ev.candidate
                s += float(np.linalg.norm(z - origin))
                accepted.append(dict(z=z.copy(), state=parent, event=ev.event, s=s))
                # At zero incremental slip re-integration may select stick; the
                # oriented predictor and event comparison use the incoming side.
                # Replay only tests equilibrium and the original acceptance gates.
                incoming = ev
                point = evaluate(parent, z)
                if (
                    not point.eligible
                    or np.linalg.norm(point.residual)
                    > controls.force_tol * point.denominator
                ):
                    raise ValueError("accepted-state replay failed")
                break
            if cut == controls.event_bisections:
                reason = cause
                break
            ds *= 0.5
        if not accepted_trial:
            break
        if stop is not None and stop(accepted):
            reason = "requested interval reached"
            break
    critical.pop("next", None)
    if rows and "record" in locals():
        critical[rows[-1]["id"]] = record
    payload = dict(
        z=z,
        parent=identity(parent),
        previous_direction=old,
        k=incoming.k,
        q=incoming.q,
        event=incoming.event,
    )
    return dict(
        termination_reason=reason,
        physical_classification="D",
        accepted=accepted,
        restart=dict(payload=payload, sha256=identity(payload)),
        last_accepted=identity(parent),
        trials=rows,
        critical=critical,
    )


def simple_point(evaluate, parent, before, after, tol):
    """Locate R=0, Kv=0 and v.v=1 between two equilibria of a smooth fixed branch."""
    a, b = np.asarray(before), np.asarray(after)
    ea, eb = evaluate(parent, a), evaluate(parent, b)
    if (
        any(
            not ev.eligible or np.linalg.norm(ev.residual) > tol * ev.denominator
            for ev in (ea, eb)
        )
        or ea.event != eb.event
    ):
        return dict(kind="D", reason="unaccepted endpoint or event")
    chord = b - a
    if not np.any(chord):
        return dict(kind="D", reason="identical endpoints")
    try:
        ta, tb = direction(ea, chord), direction(eb, chord)
    except ArithmeticError:
        return dict(kind="D", reason="non-unique endpoint tangent")
    if ta[-1] * tb[-1] >= 0:
        return dict(kind="D", reason="no reversal on accepted branch")
    n = len(a) - 1
    guess = (a + b) / 2
    # Localize on equilibria before solving Moore--Spence. The chord plane
    # stays regular through a simple load extremum.
    plane = chord / np.linalg.norm(chord)
    left, right = a.copy(), b.copy()
    for _ in range(38):
        mid = (left + right) / 2

        def arc_equations(z):
            e = evaluate(parent, z)
            return np.r_[e.residual, plane @ (z - mid)]

        def arc_jac(z):
            e = evaluate(parent, z)
            return np.vstack((np.column_stack((e.k, e.q)), plane))

        refined = root(arc_equations, mid, jac=arc_jac, tol=min(tol, 1e-13))
        z = refined.x
        e = evaluate(parent, z)
        if (
            not e.eligible
            or e.event != ea.event
            or np.linalg.norm(e.residual) > tol * e.denominator
            or abs(plane @ (z - mid)) > tol
        ):
            break
        guess = z
        if direction(e, plane)[-1] * ta[-1] > 0:
            left = z
        else:
            right = z

    def equations(y):
        ev = evaluate(parent, y[: n + 1])
        v = y[n + 1 :]
        return np.r_[ev.residual, ev.k @ v, v @ v - 1]

    sol = root(
        equations,
        np.r_[guess, np.linalg.svd(evaluate(parent, guess).k)[2][-1]],
        tol=tol,
    )
    seed = np.r_[guess, np.linalg.svd(evaluate(parent, guess).k)[2][-1]]
    y = (
        sol.x
        if np.linalg.norm(equations(sol.x)) < np.linalg.norm(equations(seed))
        else seed
    )
    z, v = y[: n + 1], y[n + 1 :]
    ev = evaluate(parent, z)
    u, sv, _ = np.linalg.svd(ev.k)
    w = u[:, -1]
    h = np.cbrt(np.finfo(float).eps) * (1 + np.linalg.norm(z[:-1]))
    probes = [evaluate(parent, z + np.r_[sign * h * v, 0.0]) for sign in (-1, 1)]
    curvature = w @ (probes[1].k - probes[0].k) @ v / (2 * h)
    scale = max(np.linalg.norm(ea.k), np.linalg.norm(eb.k), 1.0)
    aug = np.vstack((np.column_stack((ev.k, ev.q)), ta))
    along = (z - a) @ chord / (chord @ chord)
    good = (
        np.linalg.norm(ev.residual) <= tol * ev.denominator
        and np.linalg.norm(ev.k @ v) <= tol * scale
        and abs(v @ v - 1) <= tol
        and ev.eligible
        and ev.event == ea.event
        and all(p.event == ev.event for p in probes)
        and 0 <= along <= 1
        and sum(sv <= tol * scale) == 1
        and abs(w @ ev.q) > tol * np.linalg.norm(ev.q)
        and abs(curvature) > tol * scale
        and np.linalg.matrix_rank(aug) == n + 1
    )
    return dict(
        kind="simple fold" if good else "D",
        z=z,
        R=ev.residual,
        wq=w @ ev.q,
        curvature=curvature,
        singular_values=sv,
        bordered_rank=np.linalg.matrix_rank(aug),
    )


def bifurcation(evaluate, parent, critical, epsilons, parent_at, tol):
    """Re-equilibrate both mode seeds at decreasing offsets on the original equations."""
    critical = np.asarray(critical)
    ev = evaluate(parent, critical)
    info = spectrum(ev.k, ev.q, ev.residual)
    if (
        not ev.eligible
        or np.linalg.norm(ev.residual) > tol * ev.denominator
        or info["nullity"] != 1
        or info["null_load_norm"] > tol
    ):
        return dict(kind="D", branches=[])
    v = info["right"][:, -1]
    if not np.all(np.asarray(epsilons) > 0) or not np.all(np.diff(epsilons) < 0):
        raise ValueError("mode offsets must decrease toward zero")
    branches = []
    for sign in (-1, 1):
        points = []
        for eps in epsilons:
            target = sign * eps

            def equations(z):
                e = evaluate(parent, z)
                return np.r_[e.residual, v @ (z[:-1] - parent_at(z[-1])) - target]

            result = root(equations, critical + np.r_[target * v, 0.0], tol=tol)
            z = result.x
            e = evaluate(parent, z)
            main = evaluate(parent, np.r_[parent_at(z[-1]), z[-1]])
            if (
                np.linalg.norm(equations(z)) > tol
                or not e.eligible
                or e.event != ev.event
                or not main.eligible
                or main.event != ev.event
                or np.linalg.norm(main.residual) > tol * main.denominator
            ):
                break
            points.append(z)
        branches.append(points)
    distances = [
        np.array([np.linalg.norm(z - critical) for z in group]) for group in branches
    ]
    connected = len(epsilons) >= 3 and all(
        len(d) == len(epsilons) and np.all(np.diff(d) < 0) and d[-1] <= 2 * epsilons[-1]
        for d in distances
    )
    return dict(kind="bifurcation" if connected else "D", branches=branches)


def assess_path(evaluate, points, tol):
    """Recheck both sides from the same parent before trying a local fold."""
    candidates = []
    for a, b, c in zip(points, points[1:], points[2:]):
        za, zb, zc = (p["z"] for p in (a, b, c))
        if (zb[-1] - za[-1]) * (zc[-1] - zb[-1]) >= 0:
            continue
        bracket = [za, zb, zc]
        if not a["event"] == b["event"] == c["event"]:
            candidates.append(
                dict(kind="D", reason="active-set change", bracket=bracket)
            )
            continue
        try:
            found = simple_point(evaluate, a["state"], za, zc, tol)
        except TrialRejected as error:
            found = dict(kind="D", reason=str(error))
        candidates.append(dict(bracket=bracket, **found))
    return candidates
