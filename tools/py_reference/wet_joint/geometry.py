"""Affine Tet10/tri6 geometry and independent polynomial integration.

Only NumPy/SciPy and input JSON are used. Bulk elasticity uses a degree-two
four-point tetrahedron rule. Hydrostatic volume terms clip the tetrahedron
at the water plane and integrate each resulting simplex exactly to degree 2.
The contact product rule is the *input discretization*, not a continuum oracle.
"""
import copy
import itertools
import numpy as np
from scipy.special import roots_jacobi
from scipy.spatial import ConvexHull

EDGES = ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3))
FACES = ((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1))
TRI_EDGES = ((0, 1), (1, 2), (2, 0))
# Exact degree-two, positive tetrahedron rule (weights sum to 1).
A = (5 + 3*np.sqrt(5))/20
B = (5 - np.sqrt(5))/20
TET_L = np.full((4, 4), B)
np.fill_diagonal(TET_L, A)


def basis(lam):
    """Quadratic simplex basis, in vertex then cyclic-edge order."""
    lam = np.atleast_2d(lam)
    dim = lam.shape[1] - 1
    dl = np.vstack((-np.ones(dim), np.eye(dim)))
    edges = EDGES if dim == 3 else TRI_EDGES
    n = np.column_stack([lam*(2*lam-1)] + [4*lam[:, i]*lam[:, j] for i, j in edges])
    g = np.concatenate([(4*lam-1)[:, :, None]*dl,
                        np.stack([4*(lam[:, i, None]*dl[j]+lam[:, j, None]*dl[i])
                                  for i, j in edges], axis=1)], axis=1)
    return n, g


def strain_matrix(grad):
    b = np.zeros((len(grad), 6, 3*grad.shape[1]))
    for a in range(grad.shape[1]):
        x, y, z = grad[:, a].T
        k = 3*a
        b[:, 0, k] = x; b[:, 1, k+1] = y; b[:, 2, k+2] = z
        b[:, 3, k] = y; b[:, 3, k+1] = x
        b[:, 4, k+1] = z; b[:, 4, k+2] = y
        b[:, 5, k] = z; b[:, 5, k+2] = x
    return b


def tet_geometry(points, lam=TET_L):
    points = np.asarray(points)
    expected = np.array([(points[i]+points[j])/2 for i, j in EDGES])
    if not np.allclose(points[4:], expected, atol=1e-12, rtol=0):
        raise ValueError('only straight affine Tet10 is supported')
    jac = (points[1:4]-points[0]).T
    det = np.linalg.det(jac)
    if det <= 0:
        raise ValueError('nonpositive cell orientation')
    n, natural_grad = basis(lam)
    grad = natural_grad @ np.linalg.inv(jac)
    return n, strain_matrix(grad), lam @ points[:4], det/6


def wet_simplices(vertices, level):
    """Tessellate the clipped convex tetrahedron using its hull and centroid."""
    vertices = np.asarray(vertices)
    if np.all(vertices[:, 2] <= level):
        return [vertices]
    if np.all(vertices[:, 2] >= level):
        return []
    pts = [p for p in vertices if p[2] <= level]
    for i, j in itertools.combinations(range(4), 2):
        p, q = vertices[i], vertices[j]
        if (p[2]-level)*(q[2]-level) < 0:
            pts.append(p+(level-p[2])/(q[2]-p[2])*(q-p))
    pts = np.unique(np.asarray(pts), axis=0)
    center = pts.mean(axis=0)
    return [np.vstack((center, pts[face])) for face in ConvexHull(pts).simplices]


def hydro_volume(points, level, gamma):
    force = np.zeros(30)
    vertices = np.asarray(points)[:4]
    inv = np.linalg.inv((vertices[1:]-vertices[0]).T)
    for simplex in wet_simplices(vertices, level):
        xyz = TET_L @ simplex
        r = (xyz-vertices[0]) @ inv.T
        lam = np.column_stack((1-r.sum(axis=1), r))
        _, b, _, _ = tet_geometry(points, lam)
        volume = abs(np.linalg.det((simplex[1:]-simplex[0]).T))/6
        pressure = gamma*(level-xyz[:, 2])
        force -= np.einsum('qai,q->i', b[:, :3], pressure)*volume/4
    return force


def triangle_rule(order=4):
    axes = []
    for alpha in (1, 0):
        x, w = roots_jacobi(order, alpha, 0)
        axes.append(((x+1)/2, w/2**(alpha+1)))
    lam = []; weights = []
    for u, wu in zip(*axes[0]):
        for v, wv in zip(*axes[1]):
            lam.append(((1-u)*(1-v), u, (1-u)*v)); weights.append(wu*wv)
    return np.array(lam), np.array(weights)


def topology(mesh):
    result = {}
    for ei, cell in enumerate(mesh['tetrahedra10']):
        for local in FACES:
            corners = [cell[i] for i in local]
            nodes = corners + [cell[4+next(k for k, e in enumerate(EDGES)
                                                  if set(e) == {local[i], local[j]})]
                               for i, j in TRI_EDGES]
            result.setdefault(tuple(sorted(corners)), []).append((ei, nodes))
    return result


def face_geometry(xyz, nodes):
    p = xyz[nodes]
    expected = np.array([(p[i]+p[j])/2 for i, j in TRI_EDGES])
    if not np.allclose(p[3:], expected, atol=1e-12, rtol=0):
        raise ValueError('curved face is outside scope')
    area_vector = np.cross(p[1]-p[0], p[2]-p[0])
    return p, area_vector


def pressure_face(xyz, nodes, pressure):
    """Constant pressure on straight tri6, exact quadratic shape integral."""
    _, cross = face_geometry(xyz, nodes)
    # Vertex shape integrals vanish; each midside integrates to area/3.
    return np.vstack((np.zeros((3, 3)), np.tile(-pressure*cross/6, (3, 1))))


def refine(mesh):
    """Uniform 1-to-8 straight tetrahedron refinement, then Tet10 elevation.

    Node identity, not coordinate alone, keys the edges: paired mechanical
    faces remain separated. No history is transferred between mesh levels.
    """
    out = copy.deepcopy(mesh)
    xyz = [list(p) for p in mesh['nodes']]
    mids = {}
    for cell in mesh['tetrahedra10']:
        for k, (i, j) in enumerate(EDGES):
            mids[tuple(sorted((cell[i], cell[j])))] = cell[4+k]
    sets = {name: set(ids) for name, ids in mesh['node_sets'].items()}

    def midpoint(i, j):
        edge = tuple(sorted((i, j)))
        if edge not in mids:
            mids[edge] = len(xyz); xyz.append(((np.array(xyz[i])+xyz[j])/2).tolist())
        node = mids[edge]
        for ids in sets.values():
            if i in ids and j in ids:
                ids.add(node)
        return node

    cells = []; material = []; region = []
    for ei, original in enumerate(mesh['tetrahedra10']):
        a, b, c, d = original[:4]
        ab, bc, ca, ad, bd, cd = [midpoint(original[i], original[j]) for i, j in EDGES]
        children = ((a, ab, ca, ad), (ab, b, bc, bd), (ca, bc, c, cd), (ad, bd, cd, d),
                    (ab, cd, ca, ad), (ab, cd, ad, bd), (ab, cd, bd, bc), (ab, cd, bc, ca))
        for child in children:
            child = list(child)
            p = np.array(xyz)[child]
            if np.linalg.det((p[1:]-p[0]).T) < 0:
                child[1], child[2] = child[2], child[1]
            cells.append(child+[midpoint(child[i], child[j]) for i, j in EDGES])
            material.append(mesh['material_ids'][ei]); region.append(mesh['region_ids'][ei])
    faces = {}
    for name, items in mesh['face_sets'].items():
        faces[name] = []
        for a, b, c in items:
            ab, bc, ca = midpoint(a, b), midpoint(b, c), midpoint(c, a)
            faces[name].extend(((a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)))
    out.update(nodes=xyz, tetrahedra=[c[:4] for c in cells], tetrahedra10=cells,
               material_ids=material, region_ids=region, face_sets=faces,
               node_sets={name: sorted(ids) for name, ids in sets.items()})
    return out
