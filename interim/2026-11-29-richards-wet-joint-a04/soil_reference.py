"""Input-only soil Richards/construction reference; no production solves/results.

Reuse the separate C-04 full RT0 cell/trace Newton, not Rust's condensation.
3D Tet4 elasticity and every water/load integral are assembled here from inputs.
Only existing dry skeleton / frozen porosity / MC elastic branch is represented.
"""
import copy
import importlib.util
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
import fem_reference as fem
from verify_high_order import physical, constraints
from mc_reference import normalized_yield
spec=importlib.util.spec_from_file_location('a04_independent_rt0',ROOT/'interim/2026-11-22-wet-joint-transient-water/hydraulic_reference.py')
hyd=importlib.util.module_from_spec(spec);spec.loader.exec_module(hyd)

SAMPLE_KEYS=('pressure_head_m','physical_pore_pressure_kpa','bishop_chi','effective_pore_stress_kpa','effective_saturation','water_saturation','water_content','water_equivalent_mass_per_volume')


def samples(h,z,law,gamma):
    psi=np.asarray(h)-z;n=law['n'];m=1-1/n
    se=(1+(law['alpha_per_m']*np.maximum(-psi,0))**n)**(-m)
    theta=law['residual_water_content']+(law['saturated_water_content']-law['residual_water_content'])*se
    p=gamma*psi
    return np.stack((psi,p,se,se*p,se,theta/law['saturated_water_content'],theta,theta+law['saturated_specific_storage_per_m']*np.maximum(psi,0)),axis=-1)


class Model:
    def __init__(self,mesh,ground,case):
        assert not mesh.get('tetrahedra10') and len(mesh['nodes'])<=32
        assert case['excavations'][0]['bodies']==['UpperClay'] and len(case['excavations'])==1
        self.mesh,self.ground,self.case=mesh,ground,case;self.base=case['one_way'];self.flow=self.base['flow']
        assert self.base['bishop']=={'kind':'effective_saturation'} and not self.base['mechanics'].get('joint_system')
        assert not self.base['initial_exterior_water']
        assert all(not p['reservoirs'] for p in self.base['exterior_periods'][:2])
        for period,exterior in zip(self.flow['analysis']['periods'],self.base['exterior_periods']):
            for reservoir in exterior['reservoirs']:
                bc=next(b['condition'] for b in period['boundaries'] if b['face_set']==reservoir['face_set'])
                if bc['kind']!='head' or bc['value'].get('gradient',[0,0,0])!=[0,0,0] or bc['value']['head_m']!=reservoir['level_m']:
                    raise ValueError('declared reservoir and hydraulic head boundary are inconsistent')
        assert self.base['mechanics']['settings']['kinematics']=='side_rollers'
        self.xyz=np.asarray(mesh['nodes']);self.cells=np.asarray(mesh['tetrahedra']);self.nd=self.xyz.size
        self.nc=len(self.cells);self.gamma=ground['initialization']['water_unit_weight_kn_m3']
        self.mats={m['id']:m for m in ground['materials']};self.laws={m['id']:m['retention'] for m in self.flow['materials']}
        self.q=[fem.volume(mesh,e) for e in range(self.nc)]
        self.d=[fem.elastic_tensor(self.mats[i]['mohr_coulomb']['young_kpa'],self.mats[i]['mohr_coulomb']['poisson']) for i in mesh['material_ids']]
        self.ids=[(3*c[:,None]+np.arange(3)).ravel() for c in self.cells]
        self.seed=np.array([fem.initial_fields(mesh,ground,e,self.q[e]['position'])[0][0] for e in range(self.nc)])
        self.e0=np.array([fem.initial_fields(mesh,ground,e,self.q[e]['position'])[3]['void_ratio']['value'] for e in range(self.nc)])
        self.after=np.array([r!=next(r['id'] for r in mesh['regions'] if r['body']=='UpperClay') for r in mesh['region_ids']])
        self.full=self.hydraulic(np.ones(self.nc,bool));self.cut=self.hydraulic(self.after)

    def hydraulic(self,active):
        elems=np.flatnonzero(active);small=copy.deepcopy(self.mesh)
        small['tetrahedra']=[self.mesh['tetrahedra'][e] for e in elems]
        small['material_ids']=[self.mesh['material_ids'][e] for e in elems]
        small['region_ids']=[self.mesh['region_ids'][e] for e in elems]
        model=hyd.Model(small,self.flow);model.elements=elems
        return model

    def initial_hydraulic(self):
        m=self.full;h=m.seed.astype(np.longdouble);a=self.flow['initial_head']
        trace=(a['head_m']-m.datum+(m.face_centres-a.get('reference_point_m',[0.,0.,0.]))@np.array(a.get('gradient',[0.,0.,0.]))).astype(np.longdouble)
        storage,mass=m.local(h);flux=hyd.local_solve(mass,h[:,None]-trace[m.cf])
        return dict(h=h,trace=trace,storage=storage,flux=flux)

    def force(self,stress,active):
        result=np.zeros(self.nd)
        for e in np.flatnonzero(active):
            q=self.q[e];result[self.ids[e]]-=q['b'][0].T@stress[e]*q['weight'].sum()
        return result

    def stiffness(self,active):
        k=np.zeros((self.nd,self.nd))
        for e in np.flatnonzero(active):
            b=self.q[e]['b'][0];k[np.ix_(self.ids[e],self.ids[e])]+=b.T@self.d[e]@b*self.q[e]['weight'].sum()
        return k

    def fluid(self,m,v,reservoirs=()):
        active=np.zeros(self.nc,bool);active[m.elements]=True
        h=v['h']+m.datum;body=np.zeros((len(self.xyz),3));internal=np.zeros_like(body);raw=[];mp=[];positions=[];weights=[]
        for local,e in enumerate(m.elements):
            law=self.laws[self.mesh['material_ids'][e]]
            value=samples(h[local],m.position[local,:,2],law,self.gamma)
            raw.extend(value);positions.extend(m.position[local]);weights.extend(m.weight[local])
            # Each Tet4 gradient is constant, hydraulic quadrature is independent.
            vertices=self.xyz[self.cells[e]];gradient=np.vstack((-np.ones(3),np.eye(3)))@np.linalg.inv((vertices[1:]-vertices[0]).T)
            natural=(m.position[local]-vertices[0])@np.linalg.inv(vertices[1:]-vertices[0])
            bary=np.column_stack((1-natural.sum(axis=1),natural))
            body[self.cells[e],2]-=self.gamma*(m.weight[local]*value[:,7])@bary
            internal[self.cells[e]]-=gradient*(m.weight[local]@value[:,3])
            mp.append(samples(h[local],self.q[e]['position'][0,2],law,self.gamma))
        dry=physical(self.mesh,self.ground,active,self.ground['initialization']['surface_surcharge_kpa'],True)[0]
        exterior=np.zeros_like(body)
        if not active.all():
            for reservoir in reservoirs:
                for face in self.mesh['face_sets'][reservoir['face_set']]:
                    owners=fem.topology(self.mesh)[tuple(sorted(face))];live=[e for e in owners if active[e]]
                    assert len(live)==1
                    q=fem.face(self.mesh,live[0],face)
                    exterior[q['nodes']]-=q['shape'].T@(self.gamma*np.maximum(reservoir['level_m']-q['position'][:,2],0)[:,None]*q['area_vector'])
        physical_load=dry+body+exterior;rhs=physical_load-internal
        theta_s=np.array([self.laws[self.mesh['material_ids'][e]]['saturated_water_content'] for e in m.elements])
        water=float(np.sum(theta_s*m.volume+v['storage']))
        return dict(hydraulic_samples=np.asarray(raw),mechanical_samples=np.asarray(mp),hydraulic_positions=np.asarray(positions),hydraulic_weights=np.asarray(weights),
            fluid_internal=internal,water_body=body,exterior=exterior,physical=physical_load,rhs=rhs,water=water,dry=dry)

    def mechanical(self,active,load,old_u,old_stress,old_strain,old_void):
        _,fixed=constraints(self.mesh,active,'side_rollers');free=np.flatnonzero(~fixed.ravel())
        used=np.unique(self.cells[active]);k=self.stiffness(active);seed_force=self.force(self.seed,active)
        u=old_u.copy();u[used]=0.;a=k[np.ix_(free,free)];b=(load.ravel()-seed_force)[free]
        x=np.linalg.solve(a,b);linear=float(np.linalg.norm(a@x-b)/max(np.linalg.norm(b),1e-30))
        assert linear<=1e-10,linear
        u.ravel()[free]=x;strain=old_strain.copy();stress=old_stress.copy();void=old_void.copy()
        for e in np.flatnonzero(active):
            strain[e]=self.q[e]['b'][0]@u.ravel()[self.ids[e]]
            stress[e]=self.seed[e]-self.d[e]@strain[e]
            void[e]=(1+self.e0[e])*np.exp(strain[e,:3].sum())-1
            assert normalized_yield(self.mats[self.mesh['material_ids'][e]]['mohr_coulomb'],stress[e:e+1]).max()<=1e-8
        internal=self.force(stress,active).reshape(-1,3)
        residual=float(np.linalg.norm((internal-load).ravel()[free])/max(np.linalg.norm(load),1))
        assert residual<=1e-10,residual
        return u,stress,strain,void,internal,residual

    def fields(self,m,v,fluid,u,stress,strain,void,active,load,holding,removed,inflow,initial,physical_internal=None):
        hp=fluid['hydraulic_samples'];mp=fluid['mechanical_samples']
        total=stress.copy();total[m.elements,:3]+=mp[:,3,None]
        effective=self.force(stress,active).reshape(-1,3)
        fields=dict(cell_ids=m.elements,face_keys=np.array(m.keys),active_elements=active,
            cell_head_m=v['h']+m.datum,face_head_m=v['trace']+m.datum,outflow_m3_s=v['flux'],storage_offset_m3=v['storage'],
            hydraulic_positions=fluid['hydraulic_positions'],hydraulic_weights=fluid['hydraulic_weights'],
            hydraulic_darcy_flux_m_s=np.einsum('eqia,ei->eqa',m.basis,v['flux']).reshape(-1,3),
            u=u,stress=stress,strain=strain,void=void,total_stress=total,
            plastic_strain=np.zeros_like(stress),equivalent_plastic_strain=np.zeros(self.nc),plastic_work=np.zeros(self.nc),material_return=np.zeros(self.nc,int),
            effective_internal=effective,fluid_internal=fluid['fluid_internal'],water_body=fluid['water_body'],exterior=fluid['exterior'],
            physical_load=fluid['physical'],equivalent_rhs=fluid['rhs'],applied_load=load,holding_load=holding,
            reaction=effective-load,physical_reaction=effective+fluid['fluid_internal']-fluid['physical']-holding,
            water_m3=np.array(fluid['water']),removed_water_m3=np.array(removed),net_inflow_m3=np.array(inflow),mass_defect_m3=np.array(fluid['water']+removed-initial-inflow))
        for j,key in enumerate(SAMPLE_KEYS):fields['hydraulic_'+key]=hp[:,j];fields['mechanical_'+key]=mp[:,j]
        return fields

    def solve(self,sink):
        m=self.full;v=self.initial_hydraulic();fluid=self.fluid(m,v);active=np.ones(self.nc,bool)
        u=np.zeros_like(self.xyz);stress=self.seed.copy();strain=np.zeros_like(stress);void=self.e0.copy()
        seed_force=self.force(stress,active).reshape(-1,3);removed=0.;inflow=0.;initial=fluid['water'];rows=[]
        for i in range(1,5):
            progress=i/4;holding=(seed_force-fluid['rhs'])*(1-progress);load=fluid['rhs']+holding
            u,stress,strain,void,internal,res=self.mechanical(active,load,u,stress,strain,void)
            row=dict(phase='initialization',name='initial_equilibrium',time_s=0.,progress=progress)
            sink(row,self.fields(m,v,fluid,u,stress,strain,void,active,load,holding,removed,inflow,initial));rows.append(dict(**row,force_residual=res))
        dt=self.flow['settings']['maximum_time_step_s'];t=0.;old_load=load.copy()
        for index,period in enumerate(self.flow['analysis']['periods']):
            if index==2:
                old_m=m;old_v=v;m=self.cut;active=self.after.copy();lookup={k:i for i,k in enumerate(old_m.keys)}
                local=[list(old_m.elements).index(e) for e in m.elements]
                v=dict(h=old_v['h'][local].copy(),storage=old_v['storage'][local].copy(),flux=old_v['flux'][local].copy(),trace=old_v['trace'][[lookup[k] for k in m.keys]].copy())
                fluid=self.fluid(m,v,self.base['exterior_periods'][index]['reservoirs']);removed=initial+inflow-fluid['water']
                # Independently integrate removed absolute water, not ledger closure.
                dead=np.flatnonzero(~active);theta=np.array([self.laws[self.mesh['material_ids'][e]]['saturated_water_content'] for e in dead])
                removed=float(np.sum(theta*old_m.volume[dead]+old_v['storage'][dead]))
                removed_force=self.force(stress,~active).reshape(-1,3)
                holding0=old_load-fluid['rhs']-removed_force
                used=np.unique(self.cells[active]);unused=sorted(set(range(len(u)))-set(used));holding0[unused]=0.
                for j in range(1,5):
                    holding=holding0*(1-j/4);load=fluid['rhs']+holding
                    u,stress,strain,void,internal,res=self.mechanical(active,load,u,stress,strain,void)
                    row=dict(phase='excavation_release',name='cut',time_s=t,progress=j/4)
                    sink(row,self.fields(m,v,fluid,u,stress,strain,void,active,load,holding,removed,inflow,initial));rows.append(dict(**row,force_residual=res))
                old_load=load.copy()
            while t<period['end_time_s']:
                step=min(dt,period['end_time_s']-t)
                result=m.step(v['h'],v['trace'],v['storage'],period,step)
                t+=step;v=dict(h=result['cell_head_m']-m.datum,trace=result['face_head_m']-m.datum,storage=result['storage_offset_m3'],flux=result['outflow_m3_s'])
                inflow+=step*(result['source_m3_s']-result['boundary_outflow_m3_s'])
                fluid=self.fluid(m,v,self.base['exterior_periods'][index]['reservoirs']);holding=np.zeros_like(u);load=fluid['rhs']
                u,stress,strain,void,internal,res=self.mechanical(active,load,u,stress,strain,void)
                row=dict(phase='transient',name=period['name'],time_s=t,progress=1.)
                sink(row,self.fields(m,v,fluid,u,stress,strain,void,active,load,holding,removed,inflow,initial));rows.append(dict(**row,force_residual=res,hydraulic_scaled_residual=result['maximum_scaled_residual'],true_linear_residual=result['maximum_true_linear_residual']))
                old_load=load.copy()
        return dict(complete=True,input_only=True,independent_implementation='full RT0/trace Newton + separately assembled 3D elasticity',rows=rows,
            primary_nodal_head_unknown=False,contact=None,contact_reason='soil-only; JOINT is rejected by production')
