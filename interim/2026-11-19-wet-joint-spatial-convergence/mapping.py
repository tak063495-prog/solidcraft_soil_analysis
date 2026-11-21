"""Comparison operators only. Never transfer solved histories into a solve.

Contact: each 16-point collapsed Tri6 rule defines a tensor cubic in
u=r, v=s/(1-r). Material: elementwise total degree two least squares over
all 27 GPs. Nodal: original Tet10 shape functions. All traces retain ancestry.
"""
import numpy as np
from scipy.sparse import csr_matrix
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import study_c02 as c
fem=c.fem

CONTACT=('jump','traction','plastic_slip','accumulated','elastic_energy','plastic_dissipation','numerical_dissipation','endpoint_work','physical_pressure','gap')
MATERIAL=('stress','strain','void','pore','plastic_strain')
NODAL=('u',) # Nodal force vectors are dual functionals, not interpolated fields.

def monomials(x,degree=2):
    powers=[(i,j,k) for i in range(degree+1) for j in range(degree+1-i) for k in range(degree+1-i-j)]
    return np.array([x[:,0]**i*x[:,1]**j*x[:,2]**k for i,j,k in powers]).T

def tensor(x):
    u=x[:,0];v=x[:,1]/(1-u)
    return np.array([u**i*v**j for i in range(4) for j in range(4)]).T

def natural(vertices,points):
    delta=(vertices[1:]-vertices[0]).T
    x=np.linalg.lstsq(delta,(points-vertices[0]).T,rcond=None)[0].T
    assert np.max(abs(vertices[0]+x@delta.T-points))<1e-10
    return x

def contains(x):return (x.min(axis=1)>=-1e-10)&(x.sum(axis=1)<=1+1e-10)

def sparse(rows,columns,weights,shape):
    return csr_matrix((np.concatenate(weights),(np.concatenate(rows),np.concatenate(columns))),shape=shape)

class Maps:
    def __init__(self,meshes,case):
        self.meshes=meshes;self.case=case
        self.ground=c.fem.normalize_ground(c.read(c.ROOT/'examples/ground-joint-water-mc.json'))
        self.systems=[c.old.Water(m,c.read(c.ROOT/'examples/ground-joint-water-mc.json'),case).contact for m in meshes]
        self.material_position=[];self.material_root=[];self.contact_position=[];self.contact_root=[]
        for h,(m,s) in enumerate(zip(meshes,self.systems)):
            self.material_position.extend(q['position'] for q in (fem.volume(m,e) for e in range(len(m['tetrahedra']))))
            self.material_root.extend(np.repeat(np.arange(len(m['tetrahedra']))//8**h,27))
            self.contact_position.extend([p['position_m'] for p in s.points])
            self.contact_root.extend([p['face_index']//4**h for p in s.points])
        self.material_position=np.concatenate(self.material_position)
        self.material_root=np.array(self.material_root)
        self.contact_position=np.asarray(self.contact_position);self.contact_root=np.asarray(self.contact_root)
        fine=meshes[-1];self.node_position=np.asarray(fine['nodes']);self.node_root=np.full(len(self.node_position),-1)
        for e,cell in enumerate(fine['tetrahedra10']):
            for node in cell:
                if self.node_root[node]<0:self.node_root[node]=e//8**(len(meshes)-1)
        assert (self.node_root>=0).all()
        self.operators=[];self.audit=[]
        for h,(m,s) in enumerate(zip(meshes,self.systems)):
            nodal,_,_=self.volume_map(h,self.node_position,self.node_root,True)
            material,nearest,owners=self.volume_map(h,self.material_position,self.material_root,False)
            contact,contact_nearest=self.contact_map(h)
            self.operators.append(dict(u=nodal,material=material,contact=contact,nearest_material=nearest,nearest_contact=contact_nearest,material_owners=owners))
            self.audit.append(dict(level=h,elements=len(m['tetrahedra']),nodes=len(m['nodes']),material_gp=27*len(m['tetrahedra']),contact_gp=len(s.points),
                                   node_geometry_error_m=float(np.max(abs(nodal@np.asarray(m['nodes'])-self.node_position))),
                                   material_geometry_error_m=float(np.max(abs(material@np.concatenate([fem.volume(m,e)['position'] for e in range(len(m['tetrahedra']))])-self.material_position))),
                                   contact_geometry_error_m=float(np.max(abs(contact@np.array([p['position_m'] for p in s.points])-self.contact_position))),
                                   maximum_contact_weight_absolute=float(abs(contact.data).max()),contact_row_sum_error=float(abs(np.asarray(contact.sum(axis=1)).ravel()-1).max()),
                                   material_row_sum_error=float(abs(np.asarray(material.sum(axis=1)).ravel()-1).max())))

    def volume_map(self,h,points,roots,nodal):
        m=self.meshes[h];xyz=np.asarray(m['nodes']);cells=np.array(m['tetrahedra10']);rows=[];cols=[];weights=[];nr=[];nc=[];owners=np.full(len(points),-1)
        qnatural=fem.volume_rule(m)[0];pinv=np.linalg.pinv(monomials(qnatural))
        for e,cell in enumerate(cells):
            ids=np.flatnonzero((roots==e//8**h)&(owners<0))
            if not len(ids):continue
            x=natural(xyz[cell[:4]],points[ids]);inside=contains(x);ids=ids[inside];x=x[inside]
            if not len(ids):continue
            owners[ids]=e
            if nodal:
                w=fem.basis(x,True)[0];columns=np.broadcast_to(cell,w.shape)
            else:
                w=monomials(x)@pinv;columns=np.broadcast_to(e*27+np.arange(27),w.shape)
                nearest=((x[:,None]-qnatural[None])**2).sum(axis=2).argmin(axis=1)
                nr.append(ids);nc.append(e*27+nearest)
            rows.append(np.broadcast_to(ids[:,None],w.shape).ravel());cols.append(columns.ravel());weights.append(w.ravel())
        assert (owners>=0).all(),'unmapped physical volume point'
        result=sparse(rows,cols,weights,(len(points),len(xyz) if nodal else len(cells)*27))
        nearest=None if nodal else csr_matrix((np.ones(len(points)),(np.concatenate(nr),np.concatenate(nc))),shape=result.shape)
        return result,nearest,owners

    def contact_map(self,h):
        s=self.systems[h];xyz=np.asarray(self.meshes[h]['nodes']);rows=[];cols=[];weights=[];nr=[];nc=[];owners=np.full(len(self.contact_position),-1)
        for f in range(len(s.points)//16):
            indices=np.arange(16*f,16*(f+1));p=s.points[16*f];vertices=xyz[np.array(p['nodes_a'])[:3]]
            ids=np.flatnonzero((self.contact_root==f//4**h)&(owners<0))
            x=natural(vertices,self.contact_position[ids]);inside=contains(x);ids=ids[inside];x=x[inside]
            if not len(ids):continue
            owners[ids]=f
            q=natural(vertices,np.array([s.points[i]['position_m'] for i in indices]));w=tensor(x)@np.linalg.inv(tensor(q))
            rows.append(np.broadcast_to(ids[:,None],w.shape).ravel());cols.append(np.broadcast_to(indices,w.shape).ravel());weights.append(w.ravel())
            nearest=((x[:,None]-q[None])**2).sum(axis=2).argmin(axis=1);nr.append(ids);nc.append(indices[nearest])
        assert (owners>=0).all(),'unmapped physical contact point'
        shape=(len(self.contact_position),len(s.points))
        return sparse(rows,cols,weights,shape),csr_matrix((np.ones(len(self.contact_position)),(np.concatenate(nr),np.concatenate(nc))),shape=shape)

    def apply(self,h,fields,nearest=False):
        ops=self.operators[h];result={}
        for key in NODAL:result[key]=ops[key]@fields[key]
        for key in MATERIAL:result[key]=ops['nearest_material' if nearest else 'material']@fields[key]
        for key in CONTACT:result[key]=ops['nearest_contact' if nearest else 'contact']@fields[key]
        # On this audited elastic branch these two evaluations are exact FE/input
        # evaluations, so their difference from GP reconstruction is measurable.
        ground=self.ground;seeds=ground['initialization']['geostatic']['seed']['states']
        mesh=self.meshes[h];region={r['id']:r['body'] for r in mesh['regions']}
        e0=np.array([seeds[region[mesh['region_ids'][e]]]['void_ratio'] for e in ops['material_owners']])
        result['void_from_strain']=(1+e0)*np.exp(result['strain'][:,:3].sum(axis=1))-1
        result['pore_physical']=ground['initialization']['water_unit_weight_kn_m3']*np.maximum(ground['initialization']['water_level_m']-self.material_position[:,2],0)
        result['open']=(result['gap']>0).astype(float)
        return result

    def coordinates(self):
        return dict(contact=dict(position_m=self.contact_position,root_face=self.contact_root),material=dict(position_m=self.material_position,root_element=self.material_root),node=dict(position_m=self.node_position,root_element=self.node_root),levels=self.audit)
