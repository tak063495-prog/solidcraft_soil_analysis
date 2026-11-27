"""Reconstruct the production Tet10 elastic seed tangent on pinned geometry.

GP tangent import is an explicit tensor snapshot, not a continuum-law surrogate
for the algorithmic plastic tangent. No nonlinear solve is performed here.
"""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix, save_npz, load_npz

ROOT = Path(__file__).resolve().parents[3]
# Existing reference utilities use sibling, top-level imports.
if str(ROOT / 'tools') not in sys.path:
    sys.path.insert(0, str(ROOT / 'tools'))
import fem_reference as fem
from prepare_slope_h_series import slope_mesh, GEOMETRY
from griffiths_2d import digest, GEOMETRY_DIGEST


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def free_dofs(mesh, case=None):
    """Only homogeneous Cartesian supports; reject unsupported couplings."""
    fixed = np.zeros((len(mesh['nodes']), 3), bool)
    if case is None:
        fixed[:, 1] = True
        fixed[mesh['node_sets']['fixed_x'], 0] = True
        fixed[mesh['node_sets']['fixed_z'], 2] = True
    else:
        if any(case.get(name) for name in ('joint_system', 'anchor_system', 'beam_system', 'plate_system',
                                           'structures', 'ties', 'initially_inactive_bodies')):
            raise ValueError('only an uncoupled displacement system is supported')
        if case['settings']['kinematics'] != 'explicit':
            raise ValueError('explicit Cartesian kinematics required')
        boundary = case['boundary']
        if boundary.get('nodal_loads') or boundary.get('surface_loads'):
            raise ValueError('gravity-only RHS required; import a frozen system for other loads')
        for c in boundary['constraints']:
            d = np.asarray(c['direction'], float)
            axes = np.flatnonzero(d)
            if (d.shape != (3,) or len(axes) != 1 or abs(d[axes[0]]) != 1
                    or c['value_m'] != 0 or np.any(c.get('gradient', [0, 0, 0]))):
                raise ValueError('homogeneous Cartesian constraints required')
            fixed[mesh['node_sets'][c['node_set']], axes[0]] = True
    return np.flatnonzero(~fixed.ravel())


def assemble(mesh, young=100000., poisson=.3, gamma=20., case=None, tangents=None):
    cells = fem.cells(mesh)
    if len(cells) > 1200:
        raise ValueError('reconstruction limited to 1200 elements (xy1); xy2 belongs on the PC')
    if not np.isfinite([young, poisson, gamma]).all() or young <= 0 or not -1 < poisson < .5 or gamma <= 0:
        raise ValueError('invalid elastic material or gravity')
    free = free_dofs(mesh, case)
    if not len(free):
        raise ValueError('no free degrees of freedom')
    lookup = np.full(3*len(mesh['nodes']), -1, int)
    lookup[free] = np.arange(len(free))
    geom = fem.Geometry(mesh)
    nq = len(fem.volume_rule(mesh)[1])
    if tangents is not None:
        tangents = np.asarray(tangents, float)
        if tangents.shape != (len(cells), nq, 6, 6) or not np.isfinite(tangents).all():
            raise ValueError('one finite 6x6 algorithmic tangent per element/GP required')
    elastic = fem.elastic_tensor(young, poisson)
    rows, cols, values = [], [], []
    force = np.zeros(3*len(mesh['nodes']))
    volume = 0.
    for ei, cell in enumerate(cells):
        q = geom.volume(ei)
        b, weight = q['b'], q['weight']
        if tangents is None:
            ke = np.einsum('qai,ab,qbj,q->ij', b, elastic, b, weight, optimize=True)
        else:
            ke = np.einsum('qai,qab,qbj,q->ij', b, tangents[ei], b, weight, optimize=True)
        dofs = (3*np.asarray(cell)[:, None]+np.arange(3)).ravel()
        force[dofs[2::3]] += -gamma*(weight @ q['shape'])
        local = lookup[dofs]
        keep = np.flatnonzero(local >= 0)
        indices = local[keep]
        rows.extend(np.repeat(indices, len(indices)))
        cols.extend(np.tile(indices, len(indices)))
        values.extend(ke[np.ix_(keep, keep)].ravel())
        volume += float(weight.sum())
    matrix = coo_matrix((values, (rows, cols)), shape=(len(free), len(free))).tocsr()
    matrix.sum_duplicates()
    # Keep structural zeros: the Rust ILU(0) pattern includes the element clique.
    return matrix, force[free], free, dict(
        elements=len(cells), nodes=len(mesh['nodes']), integration_points=len(cells)*nq,
        volume_m3=volume, young_kpa=young, poisson=poisson, gamma_kn_m3=gamma,
        tangent_kind='imported_algorithmic_GP_snapshot' if tangents is not None else 'elastic_zero_strain_seed',
        plastic_failure_tangent_reproduced=False)


def griffiths(depth=1):
    if depth != 1:
        raise ValueError('cloud profile uses xy1-ny1 only')
    geometry = read(GEOMETRY)
    if digest(geometry) != GEOMETRY_DIGEST:
        raise ValueError('pinned Griffiths corner geometry changed')
    mesh, _ = slope_mesh(geometry, 1, depth)
    return mesh


def save_system(folder, matrix, rhs, free, metadata):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=False)
    save_npz(folder/'matrix.npz', matrix)
    np.savez(folder/'vectors.npz', rhs=rhs, free_dofs=free)
    metadata = dict(metadata, matrix_sha256=sha(folder/'matrix.npz'), vectors_sha256=sha(folder/'vectors.npz'))
    (folder/'metadata.json').write_text(json.dumps(metadata, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def load_system(folder):
    folder = Path(folder)
    meta = read(folder/'metadata.json')
    for name in ('matrix', 'vectors'):
        if sha(folder/f'{name}.npz') != meta[f'{name}_sha256']:
            raise ValueError('frozen system bytes changed')
    a = csr_matrix(load_npz(folder/'matrix.npz'))
    with np.load(folder/'vectors.npz', allow_pickle=False) as data:
        b, free = data['rhs'], data['free_dofs']
    if (a.shape != (len(b), len(b)) or len(free) != len(b) or b.ndim != 1
            or free.ndim != 1 or not np.issubdtype(free.dtype, np.integer)
            or np.any(free < 0) or len(np.unique(free)) != len(free)
            or not np.isfinite(a.data).all() or not np.isfinite(b).all()):
        raise ValueError('invalid frozen sparse system')
    return a, b, free, meta
