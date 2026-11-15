"""The bounded crossed matrix, frozen before observing the solutions."""
import math
import numpy as np
from b01follow_common import EDGES, mm, rigid, rotate

SHAPES = ('regular', 'trapezoid-0.2', 'trapezoid-0.4', 'parallelogram-0.2',
          'parallelogram-0.4', 'interior-0.06', 'interior-0.12', 'faceted-hinge-0.2')
FIELDS = ('rigid-x', 'rigid-y', 'rigid-z', 'membrane', 'bending', 'load')


def rotations():
    def axis(a, degrees):
        t = math.radians(degrees); c, s = math.cos(t), math.sin(t)
        if a == 0: return np.array([[1,0,0],[0,c,-s],[0,s,c]], float)
        if a == 1: return np.array([[c,0,s],[0,1,0],[-s,0,c]], float)
        return np.array([[c,-s,0],[s,c,0],[0,0,1]], float)
    return [('identity', np.eye(3)), ('rx37', axis(0,37)), ('rz61', axis(2,61)),
            ('general', mm(mm(axis(2,29),axis(1,-31)),axis(0,17)))]


def patch(kind):
    if kind not in SHAPES:
        raise ValueError('fixed eight-shape matrix only')
    uv = np.array([(i/4,j/4) for i in range(5) for j in range(5)])
    x = np.column_stack((uv, np.zeros(25)))
    for i in range(0,5,2):
        for j in range(0,5,2):
            p = 5*i+j; u,v = uv[p]
            if kind.startswith('trapezoid'):
                x[p,0] = (1-float(kind.split('-')[1])*v)*u
            elif kind.startswith('parallelogram'):
                x[p,0] = u+float(kind.split('-')[1])*v
            elif kind.startswith('interior') and i == j == 2:
                delta = float(kind.split('-')[1]); x[p,:2] += [delta,-2*delta/3]
            elif kind.startswith('faceted'):
                x[p,2] = .2*abs(u-.5)
    elements = []
    for i in range(2):
        for j in range(2):
            for tri in [((2*i,2*j),(2*i+2,2*j),(2*i+2,2*j+2)),
                        ((2*i,2*j),(2*i+2,2*j+2),(2*i,2*j+2))]:
                corners = [5*a+b for a,b in tri]
                mids = [5*((tri[a][0]+tri[b][0])//2)+(tri[a][1]+tri[b][1])//2 for a,b in EDGES]
                for mid,(a,b) in zip(mids,EDGES):
                    x[mid] = (x[corners[a]]+x[corners[b]])*.5
                elements.append(corners+mids)
    boundary = [i for i,p in enumerate(uv) if any(v in (0.,1.) for v in p)]
    return x, elements, boundary


def target_fields(x, q):
    a,b = np.asarray(x)[:,:2].T
    membrane = np.zeros((len(x),6)); membrane[:,0] = .001*a+.0002*b
    membrane[:,1] = -.0001*a+.0007*b; membrane[:,5] = -.00015
    bending = np.zeros_like(membrane)
    bending[:,2] = .5*(.001*a*a+.0006*a*b-.0007*b*b)
    bending[:,3] = .0003*a-.0007*b; bending[:,4] = -.001*a-.0003*b
    r = rigid(x)
    fields = {name:rotate(r[:,i+3].reshape(-1,6)*.001,q).ravel() for i,name in enumerate(FIELDS[:3])}
    fields.update(membrane=rotate(membrane,q).ravel(), bending=rotate(bending,q).ravel())
    return fields


def rejection_controls():
    x, elements, _ = patch('regular'); p = x[elements[0]].copy()
    rows = []
    for name, points in [('curved-edge-inplane',p.copy()), ('outplane-midpoint',p.copy())]:
        points[3,1 if name.startswith('curved') else 2] += .02
        rows.append((name,points))
    u,v = p[:,:2].T
    rows.append(('cylinder',np.column_stack((u,2*np.sin(v/2),2*(1-np.cos(v/2))))))
    rows.append(('sphere',np.column_stack((2*np.sin(u/2)*np.cos(v/2),2*np.sin(v/2),2*np.cos(u/2)*np.cos(v/2)-2))))
    return rows


def specs():
    return dict(passed=False,b01_complete=False,diagnostic_only=True,shapes=list(SHAPES),
                rotations=[dict(id=n,matrix=q.tolist()) for n,q in rotations()],
                rotation_order='column vectors: Rx17 then Ry(-31) then Rz29 for general',
                patches=32,elements_per_patch=8,nodes_per_patch=25,dofs_per_patch=150,
                fields=list(FIELDS),new_global_solves=192,
                midpoint_policy='rebuild arithmetic midpoints after corner distortion; no independent finite movement',
                scope_expanded=False,original_zero_mode_observation_threshold=1e-10)
