"""Read-only finite-precision branch diagnostics; no stability/acceptance verdict."""
import numpy as np


def component_bounds(holding, remainder):
    h, r = np.asarray(holding), np.asarray(remainder)
    if h.shape != r.shape or not np.all(np.isfinite(h)) or not np.all(np.isfinite(r)):
        raise ValueError('finite matching holding vectors required')
    allowance = 5e-14*np.maximum(1., np.abs(h))
    if (np.any(np.abs(r) > np.abs(h)+allowance)
            or np.any(r*h < -allowance)):
        raise ValueError('original signed holding component bounds exceeded')
    # The relative allowance must never authorize a load at an original zero.
    if np.any(r[h == 0.] != 0.):
        raise ValueError('original zero holding component changed')


def diagnose(k, load_direction, tangent, previous_tangent):
    k, load_direction = np.asarray(k), np.asarray(load_direction)
    left, singular, right = np.linalg.svd(k)
    augmented = np.linalg.svd(np.column_stack((k, -load_direction)), compute_uv=False)
    alignment = abs(float(left[:, -1]@load_direction))/max(float(np.linalg.norm(load_direction)), 1e-30)
    return dict(smallest_singular_kn_m=float(singular[-1]), largest_singular_kn_m=float(singular[0]),
                condition_number=float(singular[0]/singular[-1]) if singular[-1] else None,
                augmented_smallest_singular_kn_m=float(augmented[-1]),
                weak_left_load_alignment=alignment,
                tangent_coordinate_rate=float(tangent[-1]),
                turning_candidate=bool(tangent[-1]*previous_tangent[-1] < 0.),
                bifurcation_confirmed=False,
                weak_right_mode=right[-1].tolist(), weak_left_mode=left[:, -1].tolist(),
                singular_pair_right_residual=float(np.linalg.norm(k@right[-1]-singular[-1]*left[:, -1])),
                singular_pair_left_residual=float(np.linalg.norm(k.T@left[:, -1]-singular[-1]*right[-1])))


def probe_tangent(tangent, mode, amplitude):
    t = np.asarray(tangent, dtype=float)
    v = np.r_[np.asarray(mode, dtype=float), 0.]
    v -= t*(t@v)
    norm = np.linalg.norm(v)
    if norm < 1e-12:
        raise ValueError('weak mode is parallel to path tangent; transverse seed unavailable')
    v /= norm
    candidate = t+amplitude*v
    return candidate/np.linalg.norm(candidate)
