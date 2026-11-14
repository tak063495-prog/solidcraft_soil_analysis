"""Input-only independent 3D elasticity for the EXISTING hydrostatic BIRTH scope.

No native histories, accepted fractions or results enter solve(). Hydraulic
storage, flux and BIRTH water transfer are undefined here, never filled with 0.
"""
import copy
import numpy as np
import fem_reference as fem
from joint_water_reference import Water
from joint_ext_reference import constraint_space, elastic_operators
from verify_high_order import physical
from mc_reference import normalized_yield
from wet_joint_cycle_reference import contact_work
from types import SimpleNamespace
from fields import row


class Model:
    def __init__(self, mesh, ground, case):
        self.mesh, self.ground, self.case = mesh, ground, case
        self.xyz = np.asarray(mesh['nodes']); self.nd = self.xyz.size
        if self.nd > 350: raise ValueError('small dense independent reference budget')
        self.level = ground['initialization']['water_level_m']
        self.gamma = ground['initialization']['water_unit_weight_kn_m3']
        definition = case['joint_system']['joints'][0]
        if (len(case['joint_system']['joints']) != 1 or definition.get('crack_tip')
            or definition['birth_policy'] != 'stress_free' or definition['water'] != {'kind':'continuous'}
            or case['initially_inactive_bodies'] != ['UpperClay']
            or case.get('beam_system') or case.get('anchor_system')
            or [s['kind'] for s in case['stages']] != ['birth','joint_birth','pressure']):
            raise ValueError('existing fully submerged single MC body / CONTINUOUS B only')
        birth = case['stages'][0]
        assert birth['placement']=='cad_new_nodes' and birth['bodies']==['UpperClay']
        assert birth['hydraulic_initialization']==dict(kind='hydrostatic',level_m=self.level)
        assert self.level >= self.xyz[:,2].max()
        self.water = Water(mesh,ground,case); assert self.water.soil.all()
        self.system = self.water.contact
        self.ops = elastic_operators(mesh,ground)
        self.sample = [fem.volume(mesh,e) for e in range(len(self.ops))]
        regions = {r['id']:r['body'] for r in mesh['regions']}
        self.names = [regions[r] for r in mesh['region_ids']]
        self.initial_active = np.asarray([b!='UpperClay' for b in self.names])
        self.states = ground['initialization']['geostatic']['seed']['states']
        assert ground['initialization']['geostatic']['seed']['kind']=='body_states'
        assert ground['initialization']['geostatic'].get('initial_boundary_factor',0)==0
        self.faces = sorted(fem.topology(mesh).items())
        self.pore = np.concatenate([self.gamma*(self.level-q['position'][:,2]) for q in self.sample])
        self.tensor = np.asarray([kt*np.eye(3)+(kn-kt)*np.outer(n,n)
                                for kt,kn,n in zip(self.system.kt,self.system.kn,self.system.normal)])
        self.cj = np.einsum('qag,qab,qbh,q->gh',self.system.c,self.tensor,self.system.c,self.system.weight,optimize=True)

    def external(self, active, pressure):
        dry = copy.deepcopy(self.ground); dry['initialization']['water_level_m']=None
        base, _ = physical(self.mesh,dry,active,0.,True,joint_system=self.system)
        loaded, _ = physical(self.mesh,dry,active,pressure,True,joint_system=self.system)
        wet, _ = physical(self.mesh,self.ground,active,0.,True,joint_system=self.system)
        field = self.water.hydrostatic(active,self.level)
        gravity = np.zeros_like(self.xyz)
        for e in np.flatnonzero(active):
            q = self.ops[e][0]
            mid = self.mesh['material_ids'][e]
            gamma = next(m['unit_weight_kn_m3'] for m in self.ground['materials'] if m['id']==mid)
            gravity[q['nodes'],2] -= gamma*(q['weight']@q['shape'])
        exterior=[];water_faces=[]
        cfg=self.ground['initialization']['geostatic']
        hydro={tuple(sorted(f)) for name in cfg.get('hydrostatic_face_sets',[]) for f in self.mesh['face_sets'][name]}
        for face,owners in self.faces:
            live=[e for e in owners if active[e]];opposite=self.system.counterparts.get(face)
            exposed=len(live)==1 and not (opposite is not None and active[opposite])
            exterior.append(exposed);water_faces.append(exposed and (len(owners)==2 or opposite is not None or face in hydro))
        return dict(gravity_load=gravity,base_load=base-gravity,surcharge_load=loaded-base,
                    exterior_water_load=wet-base,joint_water_load=field['load'],
                    joint_water_pressure=field['pressure'],joint_face_water_force=np.asarray(field['forces']),
                    exterior_faces=np.asarray(exterior),wet_exterior_faces=np.asarray(water_faces)), (wet+loaded-base+field['load']).ravel()

    def operators(self, active, origins, seeds):
        k=np.zeros((self.nd,self.nd)); effective=np.zeros(self.nd);pore=np.zeros(self.nd)
        for e,(q,d,ids,ke) in enumerate(self.ops):
            if active[e]:
                k[np.ix_(ids,ids)]+=ke
                effective[ids]-=np.einsum('qai,a,q->i',q['b'],seeds[e],q['weight'])+ke@origins[e,ids]
                ps=np.zeros((len(q['weight']),6));ps[:,:3]=self.gamma*(self.level-q['position'][:,2,None])
                pore[ids]-=np.einsum('qai,qa,q->i',q['b'],ps,q['weight'])
        return k,effective,pore

    def fields(self,u,active,origins,seeds,e0,traction,reference,born,load,holding,external):
        stress=[];strain=[];void=[];eff=np.zeros(self.nd);pore=np.zeros(self.nd)
        for e,(q,d,ids,ke) in enumerate(self.ops):
            op=self.sample[e]
            eps=np.einsum('qai,i->qa',op['b'],u[ids]-origins[e,ids]) if active[e] else np.zeros((len(op['weight']),6))
            sig=seeds[e]-eps@d.T
            mid=self.mesh['material_ids'][e];mat=next(m['mohr_coulomb'] for m in self.ground['materials'] if m['id']==mid)
            if np.max(normalized_yield(mat,sig))>1e-8: raise ValueError('independent elastic scope reached MC yield')
            stress.extend(sig);strain.extend(eps);void.extend((1+e0[e])*np.exp(eps[:,:3].sum(axis=1))-1)
            if active[e]:
                # GL125 force assembly; native sample quadrature is only an output grid.
                eps5=np.einsum('qai,i->qa',q['b'],u[ids]-origins[e,ids]);sig5=seeds[e]-eps5@d.T
                eff[ids]-=np.einsum('qai,qa,q->i',q['b'],sig5,q['weight'])
                ps=np.zeros_like(sig5);ps[:,:3]=self.gamma*(self.level-q['position'][:,2,None])
                pore[ids]-=np.einsum('qai,qa,q->i',q['b'],ps,q['weight'])
        joint_force=np.einsum('qag,qa,q->g',self.system.c,traction,self.system.weight) if born else np.zeros(self.nd)
        internal=eff+pore+joint_force
        qn=np.einsum('qa,qa->q',traction,self.system.normal);tau=traction-qn[:,None]*self.system.normal
        jump=np.einsum('qag,g->qa',self.system.c,u) if born else np.zeros_like(traction)
        gap=np.einsum('qa,qa->q',jump-reference,self.system.normal)
        fields=dict(u=u.reshape(-1,3),stress=stress,strain=strain,void=void,
            total_stress=np.asarray(stress)+self.pore[:,None]*np.array([1,1,1,0,0,0]),
            pore_pressure=self.pore,head=np.full(len(self.pore),self.level),
            plastic_strain=np.zeros_like(stress),equivalent_plastic_strain=np.zeros(len(stress)),plastic_work=np.zeros(len(stress)),
            material_return=np.zeros(len(stress),int),active_elements=active.copy(),
            jump=jump,traction=traction,reference_jump=reference,initial_traction=np.zeros_like(traction),
            plastic_slip=np.zeros_like(traction),accumulated=np.zeros(len(traction)),plastic_dissipation=np.zeros(len(traction)),
            numerical_dissipation=np.zeros(len(traction)),endpoint_work=np.zeros(len(traction)),
            elastic_offset=reference,initial_energy=np.zeros(len(traction)),
            contact_status=np.ones(len(traction),int),contact_open=gap>0,opening=np.maximum(gap,0),gap=gap,
            normal_pressure=-qn,shear_traction=tau,elastic_energy=qn**2/(2*self.system.kn)+(tau*tau).sum(axis=1)/(2*self.system.kt),
            joint_active=[born],joint_tied=[not born],generation=[int(born)],archived_generation_count=[0],
            effective_internal=eff.reshape(-1,3),pore_internal=pore.reshape(-1,3),contact_internal=joint_force.reshape(-1,3),
            internal=internal.reshape(-1,3),applied_load=load.reshape(-1,3),reaction=(internal-load).reshape(-1,3),holding_load=holding.reshape(-1,3))
        fields.update(external)
        return {k:np.asarray(v) for k,v in fields.items()}


def solve(m,steps,sink,prepared_sink):
    """Every solve input is mesh/ground/case and an a-priori increment count."""
    assert steps>=2 and steps%2==0
    active=m.initial_active.copy();origins=np.zeros((len(m.ops),m.nd))
    seeds=np.asarray([m.states[b]['effective_stress_compression_positive_kpa'] for b in m.names],float)
    e0=np.asarray([m.states[b]['void_ratio'] for b in m.names]);u=np.zeros(m.nd)
    born=False;reference=np.zeros((len(m.system.points),3));traction=reference.copy()
    pressure=m.ground['initialization']['surface_surcharge_kpa'];rows=[];max_residual=0.
    work=np.zeros(len(traction));numerical=work.copy();previous_elastic=work.copy()
    for phase,stage in enumerate([dict(kind='geostatic',name='geostatic')]+m.case['stages']):
        kind=stage['kind'];start_pressure=pressure;target=stage.get('target_pressure_kpa',pressure)
        if kind=='birth':
            placed=u.reshape(-1,3).copy();existing={int(n) for e in np.flatnonzero(active) for n in fem.cells(m.mesh)[e]}
            for n in set(range(len(placed)))-existing:placed[n]=0.
            for p in m.system.points:
                for a,b in zip(p['nodes_a'],p['nodes_b']):
                    if a in existing and b not in existing:placed[b]=placed[a]
                    elif b in existing and a not in existing:placed[a]=placed[b]
            fresh=stage['states']['UpperClay']
            assert fresh['plastic_strain_compression_positive']==[0.]*6 and fresh['equivalent_plastic_strain']==fresh['plastic_work_kpa']==0. and fresh['last_return']=='elastic'
            for e,b in enumerate(m.names):
                if not active[e]:origins[e]=placed.ravel();seeds[e]=fresh['effective_stress_compression_positive_kpa'];e0[e]=fresh['void_ratio']
            active[:]=True;u=placed.ravel()
        elif kind=='joint_birth':
            born=True;reference=np.einsum('qag,g->qa',m.system.c,u);traction=np.zeros_like(reference)
        t,offset=constraint_space(m.mesh,m.case,active,m.system,not born)
        # All boundary factors in this bounded force-driven fixture are 0.
        assert m.case['stages'][-1]['kind']=='pressure'
        k,seed,pore=m.operators(active,origins,seeds)
        external,physical=m.external(active,pressure)
        initial_internal=k@u+seed+pore
        if born:initial_internal+=np.einsum('qag,qa,q->g',m.system.c,traction,m.system.weight)
        holding=initial_internal-physical if kind in ('geostatic','birth','joint_birth') else np.zeros(m.nd)
        if kind=='birth':
            prepared_sink(row(stage['name'],0.),m.fields(u,active,origins,seeds,e0,traction,reference,born,physical+holding,holding,external))
            sink.placement(u.reshape(-1,3))
        matrix=k+(m.cj if born else 0.)
        constant=seed+pore-(np.einsum('qag,qab,qb,q->g',m.system.c,m.tensor,reference,m.system.weight) if born else 0.)
        for index in range(1,steps+1):
            fraction=index/steps;pressure=start_pressure+(target-start_pressure)*fraction
            external,physical=m.external(active,pressure);load=physical+(1-fraction)*holding
            reduced=t.T@matrix@t;rhs=t.T@(load-constant);v=np.linalg.solve(reduced,rhs);u=t@v
            traction=np.einsum('qab,qb->qa',m.tensor,np.einsum('qag,g->qa',m.system.c,u)-reference) if born else np.zeros_like(reference)
            if born:
                pn=-np.einsum('qa,qa->q',m.system.normal,traction);tau=traction+pn[:,None]*m.system.normal
                assert pn.min()>0 and np.all(np.linalg.norm(tau,axis=1)<m.system.mu*pn),'outside closed-stick reference scope'
            fields=m.fields(u,active,origins,seeds,e0,traction,reference,born,load,(1-fraction)*holding,external)
            if born:
                elastic=fields['elastic_energy'];jump=fields['jump']-reference
                old=previous_jump if index>1 or kind!='joint_birth' else np.zeros_like(jump)
                counts=np.maximum(1,np.ceil(np.linalg.norm(jump-old,axis=1)/m.case['joint_system']['maximum_increment_m']).astype(int))
                for count in np.unique(counts):
                    mask=counts==count
                    sub=SimpleNamespace(**{key:getattr(m.system,key)[mask] for key in ('normal','kn','kt','mu')})
                    oldq=np.einsum('qab,qb->qa',m.tensor[mask],old[mask])
                    pd,num,increment=contact_work(sub,old[mask],oldq,jump[mask],int(count))
                    assert np.max(np.abs(pd),initial=0)<1e-10
                    work[mask]+=increment;numerical[mask]+=num
                previous_jump=jump.copy();previous_elastic=elastic.copy()
                fields['endpoint_work']=work.copy();fields['numerical_dissipation']=numerical.copy()
            residual=float(np.linalg.norm(t.T@(fields['internal'].ravel()-load))/max(np.linalg.norm(load),1.));max_residual=max(max_residual,residual)
            assert residual<1e-10
            sink(row(stage['name'],fraction),fields)
            rows.append(dict(stage=stage['name'],progress=fraction,independent_relative_force_residual=residual))
    return dict(reference='fresh input-only GL125/Tri6/SVD 3D elastic equilibrium',accepted=len(rows),rows=rows,
                maximum_independent_relative_force_residual=max_residual,water_balance=dict(flux=None,storage=None,birth_water_transfer=None,acceptance=None,reason='prescribed hydrostatic field; no water conservation equation is solved'))
