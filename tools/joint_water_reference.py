"""Independent pressure sources, wet-face integration and dense Darcy FE.

Only input geometry/policies generate hydraulic references. Stored mechanical
stresses/contact histories are audited, not an independent plastic 3D solve.
"""
import copy
import numpy as np
import fem_reference as fem
import darcy_reference
from hydrostatic_reference import face_force
from joint_reference import System


class Water:
    def __init__(self,mesh,ground,case,layout=None):
        self.mesh=mesh;self.ground=ground;self.case=case;self.contact=System(mesh,case['joint_system'],layout)
        self.definitions={d['id']:d for d in case['joint_system']['joints']};self.gamma=ground['initialization']['water_unit_weight_kn_m3']
        roles=set(case['joint_system'].get('structure_bodies',[]));regions={r['id']:r['body'] for r in mesh['regions']};self.soil=np.asarray([regions[r] not in roles for r in mesh['region_ids']]);self.faces=list({(p['joint_id'],p['face_index']):p for p in self.contact.points}.values())

    def level(self,p,current):
        policy=self.definitions[p['joint_id']]['water']
        return policy.get('level_m') if policy['kind']=='prescribed' and policy.get('level_m') is not None else current

    def package(self,pressure,forces,active):
        active=np.asarray(active);load=np.zeros_like(self.mesh['nodes'],dtype=float)
        for p,force in zip(self.faces,forces):
            if active[p['owner_a']] and active[p['owner_b']]:load[p['nodes_a']]+=force;load[p['nodes_b']]-=force
        return dict(pressure=np.asarray(pressure),forces=forces,load=load)

    def hydrostatic(self,active,current):
        active=np.asarray(active);xyz=np.asarray(self.mesh['nodes']);forces=[]
        for p in self.faces:
            h=self.level(p,current)
            force=face_force(xyz[p['nodes_a']],h,self.gamma)[0] if h is not None and active[p['owner_a']] and active[p['owner_b']] else np.zeros((len(p['nodes_a']),3))
            forces.append(force)
        pressure=[self.gamma*max(self.level(p,current)-p['position_m'][2],0.) if self.level(p,current) is not None and active[p['owner_a']] and active[p['owner_b']] else 0. for p in self.contact.points]
        return self.package(pressure,forces,active)

    def retain(self,field,active):
        pressure=[v if active[p['owner_a']] and active[p['owner_b']] else 0. for p,v in zip(self.contact.points,field['pressure'])]
        forces=[v if active[p['owner_a']] and active[p['owner_b']] else np.zeros_like(v) for p,v in zip(self.faces,field['forces'])]
        return self.package(pressure,forces,active)

    def steady(self,case,active,material_ids,current):
        """Build head equivalence from input faces, then independently assemble.

        Names used for prescribed faces differ from the Rust adapter. Neither
        its hydraulic mesh nor its map/head solution generates this reference.
        """
        active=np.asarray(active,dtype=bool);soil_active=active&self.soil;xyz=np.asarray(self.mesh['nodes']);parent=list(range(len(xyz)));flow=copy.deepcopy(case);extra={}
        def root(n):
            while parent[n]!=n:n=parent[n]
            return n
        for p in self.faces:
            if not(active[p['owner_a']] and active[p['owner_b']]):continue
            policy=self.definitions[p['joint_id']]['water']
            if policy['kind']=='continuous':
                for a,b in zip(p['nodes_a'],p['nodes_b']):ra,rb=root(a),root(b);parent[max(ra,rb)]=min(ra,rb)
            if policy['kind']=='prescribed':
                h=self.level(p,current);assert h is not None
                for side in ('a','b'):
                    if self.soil[p['owner_'+side]]:
                        name='reference:'+p['joint_id']+':'+side
                        extra.setdefault(name,(h,[]))[1].append(p['nodes_'+side][:3])
        roots=np.asarray([root(n) for n in range(len(xyz))]);unique,mapping=np.unique(roots,return_inverse=True);mesh=copy.deepcopy(self.mesh);mesh['nodes']=xyz[unique].tolist()
        for key in ('tetrahedra','tetrahedra10'):
            if mesh.get(key):mesh[key]=mapping[np.asarray(mesh[key])].tolist()
        mesh['node_sets']={n:np.unique(mapping[s]).tolist() for n,s in mesh['node_sets'].items()}
        mesh['face_sets']={n:[list(f) for f in sorted({tuple(sorted(mapping[f])) for f in s})] for n,s in mesh['face_sets'].items()}
        for n,(h,faces) in extra.items():mesh['face_sets'][n]=mapping[faces].tolist();flow['heads'].append(dict(face_set=n,head_m=h))
        solved=darcy_reference.solve(mesh,flow,soil_active.tolist(),material_ids)
        field={n:solved[n][mapping] for n in ('head','relative_head','datum')};field['defined']=np.isin(mapping,solved['used'])
        return dict(field=field,mapping=mapping,hydraulic_mesh=mesh,case=flow,solution=solved,soil_active=soil_active)

    def from_heads(self,active,field,current):
        active=np.asarray(active);xyz=np.asarray(self.mesh['nodes'])
        def pressure(p,shape,z):
            policy=self.definitions[p['joint_id']]['water'];values={}
            for side in ('a','b'):
                if self.soil[p['owner_'+side]]:
                    nodes=np.asarray(p['nodes_'+side]);h=field['relative_head'][nodes];assert field['defined'][nodes].all()
                    values[side]=self.gamma*((field['datum'][nodes[0]]-z)+h[0]+shape@(h-h[0]))
            if policy['kind']=='prescribed':v=self.gamma*max(self.level(p,current)-z,0.);assert all(abs(v-p)<1e-8 for p in values.values())
            elif policy['kind']=='soil':v=values[policy['soil_side']]
            else:assert abs(values['a']-values['b'])<1e-8;v=values['a']
            assert v>=-1e-8;return max(v,0.)
        forces=[]
        for p in self.faces:
            nodes=p['nodes_a'];value=np.zeros((len(nodes),3))
            if active[p['owner_a']] and active[p['owner_b']]:
                if not(self.soil[p['owner_a']] or self.soil[p['owner_b']]):value=face_force(xyz[nodes],self.level(p,current),self.gamma)[0]
                else:
                    # Different volume/surface quadrature implementation from Rust.
                    q=fem.face(self.mesh,p['owner_a'],nodes[:3])
                    assert q['nodes'].tolist()==nodes
                    ps=np.asarray([pressure(p,n,x[2]) for n,x in zip(q['shape'],q['position'])]);value=-np.einsum('qa,qb,q->ab',q['shape'],q['area_vector'],ps)
            forces.append(value)
        ps=[pressure(p,np.asarray(p['shape']),p['position_m'][2]) if active[p['owner_a']] and active[p['owner_b']] else 0. for p in self.contact.points]
        return self.package(ps,forces,active)

    def audit(self,state,report,expected):
        native=state['joint_water'];active=state['active_elements'];contact_force=self.contact.audit(state,report)
        assert report['water']==native
        assert len(native['face_nodal_force_a_kn'])==len(expected['forces']) and len(native['physical_pressure_kpa'])==len(expected['pressure'])
        assert np.allclose(native['nodal_load_kn'],expected['load'].ravel(),rtol=1e-9,atol=1e-9)
        assert np.allclose(native['physical_pressure_kpa'],expected['pressure'],rtol=1e-10,atol=1e-10)
        assert np.allclose(report['water']['nodal_load_kn'],native['nodal_load_kn'],rtol=0,atol=0)
        for native_force,reference in zip(native['face_nodal_force_a_kn'],expected['forces']):assert np.allclose(native_force,reference,rtol=1e-9,atol=1e-9)
        assert np.allclose(native['resultant_force_kn'],expected['load'].sum(axis=0),rtol=1e-9,atol=1e-8)
        assert np.allclose(native['resultant_moment_kn_m'],np.cross(np.asarray(self.mesh['nodes']),expected['load']).sum(axis=0),rtol=1e-9,atol=1e-8)
        for i,(p,r) in enumerate(zip(self.contact.points,report['points'])):
            w=r['hydraulic'];effective=r['normal_pressure_kpa'] if r['active'] else 0.;policy=self.definitions[p['joint_id']]['water'];assert w['source']==policy['kind']==native['sources'][i]
            assert abs(w['physical_pore_pressure_kpa']-expected['pressure'][i])<1e-9
            assert abs(w['effective_contact_pressure_kpa']-effective)<1e-9
            assert abs(w['total_normal_pressure_kpa']-effective-expected['pressure'][i])<1e-9
            assert abs(w['friction_limit_kpa']-self.contact.mu[i]*max(effective,0.))<1e-9
        return contact_force
