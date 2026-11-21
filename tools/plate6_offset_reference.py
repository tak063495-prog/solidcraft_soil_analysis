"""NEW reference scope: straight flat Tri6 with a neutral-plane offset.

The existing plate_reference, plate6_reference and beam_reference scopes are
unchanged. This module explicitly composes the planar polynomial reference
with a rigid-arm virtual-work map; it never consumes Rust element matrices.
"""
import copy
import numpy as np
from plate6_reference import Element as Planar
from beam_reference import System as BaseSystem

class Element:
    def __init__(self,points,section,e=0.):
        points=np.asarray(points,float)
        if not np.isfinite(e) or abs(e)>section['thickness_m']/2:
            raise ValueError('offset reference scope: finite |e| <= thickness/2')
        planar=Planar(points,section)  # retains the old straight/flat rejection
        self.e=float(e);self.normal=planar.axes[2];self.area=planar.area;self.axes=planar.axes
        self.neutral=Planar(points-e*self.normal,section) if e else planar
        self.arm=e*self.normal
        self.map=np.eye(36)
        for i in range(6):
            for a in range(3):self.map[6*i:6*i+3,6*i+3+a]=-np.cross(np.eye(3)[a],self.arm)
        self.stiffness=self.map.T@self.neutral.stiffness@self.map if e else planar.stiffness.copy()

    def load(self,traction):
        f=self.neutral.load(traction)
        return self.map.T@f if self.e else f

    def recover(self,surface):
        return self.neutral.recover(self.map@surface if self.e else surface)

    def attachment_translation(self,l,surface):
        # Independently evaluate the three displacement polynomials at x_s;
        # reconstruct from neutral kinematics and the signed rigid arm.
        s=np.asarray(surface).reshape(6,6);p=(self.map@s.ravel()).reshape(6,6)
        n,*_=self.neutral.shapes(*(np.asarray(l)@self.neutral.xy[:3]))
        return n@p[:,:3]+np.cross(n@s[:,3:],self.arm)

class System(BaseSystem):
    def __init__(self,mesh,cfg):
        if not cfg.get('plates') or any(p.get('formulation')!='mindlin_tri6_offset' for p in cfg['plates']):
            raise ValueError('new offset reference requires explicit offset Tri6 members only')
        stripped=copy.deepcopy(cfg);stripped['plates']=[]
        super().__init__(mesh,stripped)
        self.cfg=cfg;self.plates=cfg['plates'];self.plate_data=[]
        sections={s['id']:s for s in cfg['plate_sections']}
        shifts={}
        for p in self.plates:
            nodes=[self.node_ids[n] for n in p['nodes']];section=sections[p['section_id']]
            kernel=Element(self.positions[nodes],section,p.get('attachment_offset_m',0.))
            for n in nodes:
                if n in shifts and np.linalg.norm(shifts[n]-kernel.arm)>64*np.finfo(float).eps*section['thickness_m']:
                    raise ValueError('offset steps/folds unverified')
                shifts[n]=kernel.arm
            self.plate_data.append((nodes,section,kernel))

    def assemble(self,state,factor):
        result=super().assemble(state,factor)
        # The old planar assembly has zero plate load moments. The new scope
        # must retain the rigid-arm moments in the external resultant too.
        for i,(nodes,section,kernel) in enumerate(self.plate_data):
            if not state['active_plates'][i]:continue
            p=self.plates[i]
            traction=np.array(p.get('constant_traction_kpa',[0.]*3),float)+factor*np.asarray(p.get('boundary_traction_kpa',[0.]*3))
            traction[2]-=section['unit_weight_kn_m3']*section['thickness_m']
            result[5][:]+=kernel.load(traction).reshape(6,6)[:,3:].sum(axis=0)
        return result
