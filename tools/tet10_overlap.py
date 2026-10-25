"""Physical-volume displacement norms across affine Tet10 meshes.

Intersect tetrahedra as convex halfspace polyhedra. Integrate the squared
difference of the two quadratic fields on each intersection. A weighted QR
factor of a common polynomial basis avoids subtracting nearly equal energy
integrals. No geometry jitter, displacement clipping or GP stress transfer.
"""
import os
for _name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[_name] = '1'
from itertools import combinations
import numpy as np
import scipy
from scipy.spatial import ConvexHull, cKDTree
import fem_reference as ref


TRIPLES = np.asarray(list(combinations(range(8), 3)))


def monomials(x):
    a, b, c = x.T
    return np.column_stack((np.ones(len(x)), a, b, c, a*a, b*b, c*c, a*b, b*c, a*c))


def linear_product(a, b):
    result = np.zeros(10)
    result[0] = a[0]*b[0]
    result[1:4] = a[0]*b[1:]+b[0]*a[1:]
    result[4:7] = a[1:]*b[1:]
    result[7:] = [a[1]*b[2]+a[2]*b[1], a[2]*b[3]+a[3]*b[2], a[1]*b[3]+a[3]*b[1]]
    return result


class AffineMesh:
    def __init__(self, mesh):
        if mesh.get('schema') != 'soil-fem-mesh/3' or mesh.get('units') != 'm-kN-kPa':
            raise ValueError('soil-fem-mesh/3 in m-kN-kPa required')
        self.mesh = mesh
        self.xyz = np.asarray(mesh['nodes'], dtype=float)
        self.cells = np.asarray(mesh['tetrahedra10'], dtype=int)
        if self.cells.shape != (len(mesh['tetrahedra']), 10) or not np.isfinite(self.xyz).all():
            raise ValueError('invalid Tet10 geometry inventory')
        if not np.array_equal(self.cells[:, :4], mesh['tetrahedra']):
            raise ValueError('quadratic/corner connectivity mismatch')
        self.corners = self.xyz[self.cells[:, :4]]
        midpoint = np.stack([(self.corners[:, a]+self.corners[:, b])/2 for a, b in ref.EDGES], axis=1)
        self.maximum_affine_deviation_m = float(abs(self.xyz[self.cells[:, 4:]]-midpoint).max())
        self.scale = max(float(np.ptp(self.xyz, axis=0).max()), 1.)
        if self.maximum_affine_deviation_m > 1e-11*self.scale:
            raise ValueError('overlap integration requires affine Tet10 geometry; curved cells are unsupported')
        transform = np.transpose(self.corners[:, 1:]-self.corners[:, :1], (0, 2, 1))
        determinant = np.linalg.det(transform)
        if np.any(determinant <= 0):
            raise ValueError('nonpositive affine tetrahedron volume')
        self.volumes = determinant/6
        inverse = np.linalg.inv(transform)
        self.gradient = np.concatenate((-inverse.sum(axis=1)[:, None, :], inverse), axis=1)
        self.offset = -np.einsum('eai,ei->ea', self.gradient, self.corners[:, 0])
        self.offset[:, 0] += 1
        norms = np.linalg.norm(self.gradient, axis=2)
        self.normals = self.gradient/norms[:, :, None]
        self.plane_offset = self.offset/norms
        self.centers = self.corners.mean(axis=1)
        self.radii = np.linalg.norm(self.corners-self.centers[:, None, :], axis=2).max(axis=1)
        self.low, self.high = self.corners.min(axis=1), self.corners.max(axis=1)

    def coefficients(self, e, center, scale):
        linear = np.column_stack((self.offset[e]+self.gradient[e]@center, self.gradient[e]*scale))
        result = np.zeros((10, 10))
        for i in range(4):
            result[:, i] = 2*linear_product(linear[i], linear[i])
            result[:4, i] -= linear[i]
        for i, (a, b) in enumerate(ref.EDGES):
            result[:, i+4] = 4*linear_product(linear[a], linear[b])
        return result


def intersection(a, ea, b, eb):
    center = (a.centers[ea]+b.centers[eb])/2
    scale = max(a.radii[ea], b.radii[eb])
    normals = np.concatenate((a.normals[ea], b.normals[eb]))
    offsets = np.r_[a.plane_offset[ea], b.plane_offset[eb]]+normals@center
    matrices = normals[TRIPLES]
    valid = abs(np.linalg.det(matrices)) > 1e-13
    vertices = np.linalg.solve(matrices[valid], -offsets[TRIPLES[valid]][..., None])[..., 0]
    tolerance = 2e-11*scale
    vertices = vertices[np.all(vertices@normals.T+offsets >= -tolerance, axis=1)]
    unique = []
    for vertex in vertices:
        if not unique or np.linalg.norm(np.asarray(unique)-vertex, axis=1).min() > tolerance:
            unique.append(vertex)
    if len(unique) < 4:
        return None
    vertices = np.asarray(unique)
    if np.linalg.svd(vertices-vertices.mean(axis=0), compute_uv=False)[-1] <= tolerance:
        return None
    # Qhull degeneracy errors remain errors: do not perturb the physical mesh.
    hull = ConvexHull(vertices)
    interior = vertices[hull.vertices].mean(axis=0)
    return center, scale, vertices, interior, hull


class Overlay:
    def __init__(self, mesh_a, mesh_b, quadrature_order=3, progress=None):
        if quadrature_order not in (3, 4, 5):
            raise ValueError('overlap quadrature order must be 3, 4 or 5')
        if any(mesh_a.get(name) != mesh_b.get(name) for name in ('units', 'axis', 'crs', 'origin_m')):
            raise ValueError('coordinate metadata differ between overlap meshes')
        self.a, self.b = AffineMesh(mesh_a), AffineMesh(mesh_b)
        a, b = self.a, self.b
        tree = cKDTree(b.centers)
        rule, weights = ref.volume_rule(dict(tetrahedra10=[list(range(10))], quadrature_order=quadrature_order))
        records = []
        candidates_count = 0
        maximum_hull_volume_error = 0.
        for ea in range(len(a.cells)):
            possible = sorted(tree.query_ball_point(a.centers[ea], a.radii[ea]+b.radii.max()+1e-10))
            for eb in possible:
                if np.any(a.high[ea] < b.low[eb]-1e-10) or np.any(b.high[eb] < a.low[ea]-1e-10):
                    continue
                candidates_count += 1
                polyhedron = intersection(a, ea, b, eb)
                if polyhedron is None:
                    continue
                center, scale, vertices, interior, hull = polyhedron
                positions, quadrature_weights = [], []
                for triangle in hull.simplices:
                    transform = vertices[triangle]-interior
                    determinant = abs(float(np.linalg.det(transform)))
                    if determinant == 0:
                        continue
                    positions.append((interior+rule@transform)/scale)
                    quadrature_weights.append(weights*determinant)
                positions, w = np.concatenate(positions), np.concatenate(quadrature_weights)
                volume = float(w.sum())
                maximum_hull_volume_error = max(maximum_hull_volume_error, abs(volume-hull.volume))
                if abs(volume-hull.volume) > 1e-10*max(volume, 1.):
                    raise ValueError('intersection tetrahedralization differs from convex hull volume')
                polynomial = monomials(positions)
                _, factor = np.linalg.qr(np.sqrt(w)[:, None]*polynomial, mode='reduced')
                records.append((ea, eb, volume, factor,
                                a.coefficients(ea, center, scale), b.coefficients(eb, center, scale)))
            if progress and (ea+1) % 50 == 0:
                progress(ea+1, len(a.cells), len(records))
        if not records:
            raise ValueError('no positive-volume intersection')
        self.ea = np.asarray([r[0] for r in records], dtype=int)
        self.eb = np.asarray([r[1] for r in records], dtype=int)
        self.volume = np.asarray([r[2] for r in records])
        self.factor = np.stack([r[3] for r in records])
        self.ca, self.cb = [np.stack([r[index] for r in records]) for index in (4, 5)]
        self.metadata = dict(quadrature_order=quadrature_order, candidate_pairs=candidates_count,
                             positive_volume_pairs=len(records), total_intersection_volume_m3=float(self.volume.sum()),
                             maximum_hull_volume_consistency_error_m3=maximum_hull_volume_error,
                             maximum_affine_deviation_m=max(a.maximum_affine_deviation_m, b.maximum_affine_deviation_m),
                             scipy=scipy.__version__, numpy=np.__version__)
        self.coverage(np.ones(len(a.cells), dtype=bool), np.ones(len(b.cells), dtype=bool))

    def coverage(self, active_a, active_b):
        active_a, active_b = np.asarray(active_a, dtype=bool), np.asarray(active_b, dtype=bool)
        if active_a.shape != (len(self.a.cells),) or active_b.shape != (len(self.b.cells),):
            raise ValueError('active-domain inventory mismatch')
        selected = active_a[self.ea] & active_b[self.eb]
        cover_a = np.bincount(self.ea[selected], weights=self.volume[selected], minlength=len(active_a))
        cover_b = np.bincount(self.eb[selected], weights=self.volume[selected], minlength=len(active_b))
        errors = [abs(cover_a-self.a.volumes*active_a), abs(cover_b-self.b.volumes*active_b)]
        relative = [err/np.maximum(mesh.volumes, 1e-30) for err, mesh in zip(errors, (self.a, self.b))]
        maximum_relative = max(float(err.max()) for err in relative)
        maximum_absolute = max(float(err.max()) for err in errors)
        if maximum_relative > 1e-8 or maximum_absolute > 1e-8:
            raise ValueError(f'overlap does not cover the same active physical domains: {maximum_relative}, {maximum_absolute}')
        return selected, dict(active_volume_m3=float(self.volume[selected].sum()),
                              maximum_element_coverage_relative_error=maximum_relative,
                              maximum_element_coverage_error_m3=maximum_absolute)

    def compare(self, displacement_a, displacement_b, active_a, active_b):
        selected, coverage = self.coverage(active_a, active_b)
        ua, ub = np.asarray(displacement_a), np.asarray(displacement_b)
        if ua.shape != self.a.xyz.shape or ub.shape != self.b.xyz.shape or not np.isfinite(ua).all() or not np.isfinite(ub).all():
            raise ValueError('invalid displacement arrays')
        ca = np.einsum('pij,pja->pia', self.ca[selected], ua[self.a.cells[self.ea[selected]]])
        cb = np.einsum('pij,pja->pia', self.cb[selected], ub[self.b.cells[self.eb[selected]]])
        difference = np.einsum('pij,pja->pia', self.factor[selected], ca-cb)
        reference = np.einsum('pij,pja->pia', self.factor[selected], cb)
        error_squared = float(np.square(difference).sum())
        reference_squared = float(np.square(reference).sum())
        return dict(**coverage, difference_squared_integral_m5=error_squared,
                    reference_squared_integral_m5=reference_squared,
                    rms_difference_m=float(np.sqrt(error_squared/coverage['active_volume_m3'])),
                    relative_l2=float(np.sqrt(error_squared/reference_squared)) if reference_squared > 0 else None)
