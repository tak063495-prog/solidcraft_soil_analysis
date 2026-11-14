"""Contact-only common-point operators retaining original parent face ancestry."""
import numpy as np
from c02diag_transfer import transfer,apply,dot


def natural(vertices,points):
    v=np.asarray(vertices,dtype=float);p=np.asarray(points,dtype=float)
    d=v[1:]-v[0]
    axes=max(((0,1),(0,2),(1,2)),key=lambda ab:abs(d[0,ab[0]]*d[1,ab[1]]-d[0,ab[1]]*d[1,ab[0]]))
    a,b=axes;det=d[0,a]*d[1,b]-d[0,b]*d[1,a]
    if abs(det)<1e-15:raise ValueError('degenerate face')
    q=p-v[0]
    r=(q[:,a]*d[1,b]-q[:,b]*d[1,a])/det
    s=(d[0,a]*q[:,b]-d[0,b]*q[:,a])/det
    x=np.column_stack((r,s))
    if np.max(abs(v[0]+r[:,None]*d[0]+s[:,None]*d[1]-p),initial=0)>1e-10:
        raise ValueError('point off face')
    return x


class ContactPlan:
    def __init__(self,mesh,system,h,positions,roots,method):
        self.size=len(system.points);self.target_size=len(positions);self.groups=[]
        owners=np.full(len(positions),-1)
        xyz=np.asarray(mesh['nodes'])
        groups={}
        for i,p in enumerate(system.points):
            groups.setdefault((p['joint_id'],p['face_index']),[]).append(i)
        for (_,face),indices in sorted(groups.items()):
            indices=np.asarray(indices);p=system.points[indices[0]]
            v=xyz[np.asarray(p['nodes_a'])[:3]]
            ids=np.flatnonzero((roots==face//4**h)&(owners<0))
            x=natural(v,positions[ids]);inside=(x.min(axis=1)>=-1e-10)&(x.sum(axis=1)<=1+1e-10)
            ids=ids[inside];x=x[inside]
            if not len(ids):continue
            owners[ids]=face
            source=np.asarray([system.points[i]['position_m'] for i in indices])
            q=natural(v,source);weights=np.asarray([system.points[i]['weight_m2'] for i in indices])
            if method=='nearest':
                matrix=np.zeros((len(ids),len(indices)))
                for row,point in enumerate(positions[ids]):
                    j=min(range(len(indices)),key=lambda k:(dot(source[k]-point,source[k]-point),int(indices[k])))
                    matrix[row,j]=1.
            else:matrix=transfer(q,weights,x,method)
            self.groups.append((ids,indices,matrix))
        if np.any(owners<0):raise ValueError('unmapped common contact point')
        self.owners=owners

    def apply(self,values):
        values=np.asarray(values,dtype=float)
        if values.shape[0]!=self.size:raise ValueError('GP count mismatch')
        output=np.empty((self.target_size,)+values.shape[1:])
        for ids,indices,matrix in self.groups:output[ids]=apply(matrix,values[indices])
        return output

    def history(self,values):
        values=np.asarray(values,dtype=float)
        if values.shape[1]!=self.size:raise ValueError('history GP count mismatch')
        output=np.empty((len(values),self.target_size)+values.shape[2:])
        for ids,indices,matrix in self.groups:
            output[:,ids]=np.moveaxis(apply(matrix,np.moveaxis(values[:,indices],1,0)),0,1)
        return output
