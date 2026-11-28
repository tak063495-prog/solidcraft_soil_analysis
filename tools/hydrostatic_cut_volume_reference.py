"""Input-only hydrostatic volume force on affine Tet4/Tet10.

Verification only: clip each existing tetrahedron at z=H. On every resulting
tetrahedron, p is linear and grad(N) is at most linear, so barycentric moments
integrate -grad(N)*p exactly (apart from floating point). Curved Tet10 are
explicitly outside this reference; no production quadrature/GP is changed.
"""
import numpy as np
import fem_reference as fem


def affine_cell(mesh, element):
    nodes = fem.cells(mesh)[element]
    points = np.asarray(mesh['nodes'], dtype=float)[nodes]
    assert np.isfinite(points).all(), 'nonfinite geometry'
    vertices = points[:4]
    jacobian = (vertices[1:]-vertices[0]).T
    assert np.linalg.det(jacobian) > 0, 'positive affine tetrahedron required'
    if len(nodes) == 10:
        mids = np.asarray([(vertices[i]+vertices[j])/2 for i, j in fem.EDGES])
        assert np.allclose(points[4:], mids, rtol=0, atol=1e-13), 'curved Tet10 outside cut reference'
    return nodes, vertices, np.linalg.inv(jacobian)


def wet_tetrahedra(vertices, level):
    """Clip four polygon faces and the horizontal cap, then fan from inside."""
    vertices = np.asarray(vertices, dtype=float)
    assert np.isfinite(level), 'finite water level required'
    if level <= vertices[:, 2].min():
        return []
    if level >= vertices[:, 2].max():
        return [vertices.copy()]
    tolerance = 1e-13*max(np.ptp(vertices, axis=0).max(), 1.)

    def unique(points):
        result = []
        for point in points:
            if not any(np.linalg.norm(point-old) <= tolerance for old in result):
                result.append(point)
        return result

    def crossing(a, b):
        value = a+(b-a)*(level-a[2])/(b[2]-a[2])
        value[2] = level
        return value

    polygons = []
    for face in fem.FACES:
        polygon = vertices[list(face)]
        clipped = []
        for a, b in zip(polygon, np.roll(polygon, -1, axis=0)):
            if a[2] <= level:
                clipped.append(a)
            if (a[2] < level < b[2]) or (b[2] < level < a[2]):
                clipped.append(crossing(a, b))
        clipped = unique(clipped)
        if len(clipped) >= 3:
            polygons.append(clipped)
    cap = unique([v for v in vertices if v[2] == level]+[
        crossing(vertices[i], vertices[j]) for i, j in fem.EDGES
        if min(vertices[i, 2], vertices[j, 2]) < level < max(vertices[i, 2], vertices[j, 2])])
    if len(cap) >= 3:
        center = np.mean(cap, axis=0)
        cap.sort(key=lambda p: np.arctan2(p[1]-center[1], p[0]-center[0]))
        polygons.append(cap)
    center = np.mean(unique([p for polygon in polygons for p in polygon]), axis=0)
    cells = []
    for polygon in polygons:
        for i in range(1, len(polygon)-1):
            tet = np.asarray([center, polygon[0], polygon[i], polygon[i+1]])
            determinant = np.linalg.det((tet[1:]-tet[0]).T)
            if determinant != 0:
                cells.append(tet)
    return cells


def force(mesh, element, level, gamma, check_order=None):
    """Return element force and wet/pressure moments, without any solve output.

    check_order uses independent GL/Duffy on the clipped cells as a second
    integration algorithm. The default uses analytic barycentric moments.
    """
    assert np.isfinite(gamma) and gamma >= 0, 'finite nonnegative unit weight required'
    nodes, vertices, inverse = affine_cell(mesh, element)
    gradients = np.vstack((-inverse.sum(axis=0), inverse))
    load = np.zeros((len(nodes), 3))
    volume = 0.; first = np.zeros(3); pressure_integral = 0.; pressure_first = np.zeros(3)
    if check_order is not None:
        assert check_order >= 3
        natural, weights = fem.volume_rule({'tetrahedra10': [()]}, True, check_order)
    for tet in wet_tetrahedra(vertices, level):
        determinant = abs(np.linalg.det((tet[1:]-tet[0]).T))
        size = determinant/6
        p = gamma*(level-tet[:, 2])
        r = (tet-vertices[0])@inverse.T
        barycentric = np.column_stack((1-r.sum(axis=1), r))
        # Integral of two affine fields: V/20 * [(sum a)(sum b)+sum(a*b)].
        integral_p = size*p.sum()/4
        integral_pl = size/20*(p.sum()*barycentric.sum(axis=0)+p@barycentric)
        volume += size; first += size*tet.mean(axis=0)
        pressure_integral += integral_p
        pressure_first += size/20*(p.sum()*tet.sum(axis=0)+p@tet)
        if check_order is None:
            if len(nodes) == 4:
                load -= gradients*integral_p
            else:
                load[:4] -= (4*integral_pl-integral_p)[:, None]*gradients
                for k, (i, j) in enumerate(fem.EDGES):
                    load[k+4] -= 4*(integral_pl[i]*gradients[j]+integral_pl[j]*gradients[i])
        else:
            xyz = tet[0]+natural@(tet[1:]-tet[0])
            r = (xyz-vertices[0])@inverse.T
            _, derivative = fem.basis(r, len(nodes) == 10)
            gradient = derivative@inverse
            pressure = gamma*(level-xyz[:, 2])
            load -= np.einsum('qai,q,q->ai', gradient, pressure, determinant*weights)
    return dict(nodes=nodes, force=load, wet_volume=volume, wet_first_moment=first,
                pressure_integral=pressure_integral, pressure_first_moment=pressure_first)
