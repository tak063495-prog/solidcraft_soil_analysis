"""Independent dry 3D MCC rate FEM for verification, not the production solver.

Uses analytic continuum elastoplastic tangents and explicit midpoint integration
with yield/unloading events. No Rust stress update, algorithmic tangent, trial
displacement, holding vector or material history generates this solution.
Stress and engineering strain are compression positive, ordered xx,yy,zz,xy,yz,xz.
Source equations: Itasca MCC 15,18,19,17 (see docs/MCC_RATE_FEM.ja.md).
"""
import os
for _name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[_name]='1'
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import fem_reference as ref
from verify_high_order import physical, constraints


def terms(state, parameters):
    """State columns: stress[6], ln(pc), ln(v), plastic engineering strain[6]."""
    stress=state[:,:6];pc=np.exp(state[:,6]);v=np.exp(state[:,7]);p=stress[:,:3].mean(axis=1)
    if not np.isfinite(state).all() or np.any(p<=0) or np.any(v<=1):
        raise ValueError('independent MCC rate path left p>0, v>1')
    m2,lam,kap,nu=parameters.T
    dev=stress.copy();dev[:,:3]-=p[:,None]
    q2=1.5*np.square(dev[:,:3]).sum(axis=1)+3*np.square(dev[:,3:]).sum(axis=1)
    yield_value=(q2+m2*p*(p-pc))/(m2*pc**2)
    normal=dev*np.array([3.,3.,3.,6.,6.,6.])
    trace=m2*(2*p-pc);normal[:,:3]+=trace[:,None]/3
    bulk=v*p/kap;shear=3*bulk*(1-2*nu)/(2*(1+nu))
    elastic=np.zeros((len(state),6,6));elastic[:,:3,:3]=(bulk-2*shear/3)[:,None,None]
    for a in range(3):elastic[:,a,a]+=2*shear;elastic[:,a+3,a+3]=shear
    dn=np.einsum('qab,qb->qa',elastic,normal)
    hardening=m2*p*pc*v/(lam-kap)*trace
    denominator=np.einsum('qa,qa->q',normal,dn)+hardening
    return yield_value,elastic,normal,dn,denominator


class System:
    constitutive_terms=staticmethod(terms)

    def __init__(self, mesh, ground, active, case=None):
        self.mesh=mesh;self.ground=ground;self.active=np.asarray(active,dtype=bool)
        self.element_ids=np.flatnonzero(self.active);self.xyz=np.asarray(mesh['nodes'])
        kinematics=(case or {}).get('settings',{}).get('kinematics','side_rollers')
        if kinematics not in ('side_rollers','explicit'):raise ValueError('unsupported reference kinematics')
        self.used,self.fixed=constraints(mesh,self.active,'side_rollers')
        self.rotation=None
        boundary=(case or {}).get('boundary') or {'constraints':[]}
        if boundary.get('constraints') or kinematics=='explicit':
            for c in boundary.get('constraints',[]):
                if c['value_m']!=0 or np.any(c.get('gradient',[0.,0.,0.])):
                    raise ValueError('independent reference currently supports homogeneous directional constraints')
            model={**case,'boundary':boundary}
            bases,error=ref.staged_constraint_space(mesh,model,dict(active_elements=self.active,displacement_m=np.zeros_like(self.xyz)))
            if error:raise ValueError('nonhomogeneous reference constraints')
            self.rotation=np.zeros((len(self.xyz),3,3));self.fixed[:]=True
            for node,basis in bases.items():
                self.rotation[node,:,:basis.shape[1]]=basis;self.fixed[node,:basis.shape[1]]=False
        self.free=np.flatnonzero(~self.fixed.reshape(-1))
        if len(self.free)>6500:raise ValueError('independent verification is limited to 6500 free DOFs')
        q=[ref.volume(mesh,int(e)) for e in self.element_ids]
        self.b=np.stack([v['b'] for v in q]);self.weight=np.stack([v['weight'] for v in q])
        self.dofs=np.stack([(v['nodes'][:,None]*3+np.arange(3)).reshape(-1) for v in q])
        self.rb=self.b
        if self.rotation is not None:
            local_rotation=self.rotation[self.dofs[:,::3]//3]
            self.rb=np.einsum('eqani,enij->eqanj',self.b.reshape(*self.b.shape[:3],-1,3),local_rotation).reshape(self.b.shape)
        lookup=np.full(self.xyz.size,-1);lookup[self.free]=np.arange(len(self.free));reduced=lookup[self.dofs]
        shape=(len(q),self.dofs.shape[1],self.dofs.shape[1])
        rows=np.broadcast_to(reduced[:,:,None],shape).reshape(-1);cols=np.broadcast_to(reduced[:,None,:],shape).reshape(-1)
        self.mask=(rows>=0)&(cols>=0);self.rows=rows[self.mask];self.cols=cols[self.mask]
        unique,self.inverse=np.unique(self.rows*len(self.free)+self.cols,return_inverse=True)
        self.sparse_rows=unique//len(self.free);self.sparse_cols=unique%len(self.free)
        materials={m['id']:m['mcc'] for m in ground['materials']}
        self.parameters=np.array([[materials[mesh['material_ids'][e]][key] for key in ('critical_state_m','lambda','kappa','poisson')]
                                  for e in self.element_ids]);self.parameters[:,0]**=2
        self.parameters=np.repeat(self.parameters,self.b.shape[1],axis=0)
        self.calls=0
        self.linear_iterations=0;self.dense_fallbacks=0

    def reduce(self, force):
        if self.rotation is not None:force=np.einsum('nij,ni->nj',self.rotation,force.reshape(-1,3)).reshape(-1)
        return force[self.free]

    def lift(self, coordinates):
        if self.rotation is None:return coordinates
        return np.einsum('nij,nj->ni',self.rotation,coordinates.reshape(-1,3)).reshape(-1)

    def linear_solve(self, local_stiffness, rhs):
        """Independent COO accumulation/diagonal PCG, audited by true residual."""
        n=len(self.free)
        values=np.bincount(self.inverse,weights=local_stiffness.reshape(-1)[self.mask],minlength=len(self.sparse_rows))
        def product(x):return np.bincount(self.sparse_rows,weights=values*x[self.sparse_cols],minlength=n)
        diagonal=np.bincount(self.sparse_rows,weights=values*(self.sparse_rows==self.sparse_cols),minlength=n)
        x=np.zeros(n);r=rhs.copy();scale=max(np.linalg.norm(rhs),1e-30)
        if scale<=1e-25:return x
        if np.all(diagonal>0):
            z=r/diagonal;p=z.copy();rz=float(r@z)
            for iteration in range(max(100,4*n)):
                ap=product(p);curvature=float(p@ap)
                if curvature<=0:break
                alpha=rz/curvature;x+=alpha*p;r-=alpha*ap
                if np.linalg.norm(r)<2e-12*scale:
                    actual=rhs-product(x)
                    if np.linalg.norm(actual)<1e-11*scale:
                        self.linear_iterations+=iteration+1;return x
                    r=actual
                z=r/diagonal;next_rz=float(r@z);p=z+(next_rz/rz)*p;rz=next_rz
        # Softening may make the continuum tangent indefinite. This bounded
        # numerical fallback makes no physical stability/uniqueness claim.
        if n>1500:raise ValueError('reference PCG failed; indefinite/ill-conditioned tangent needs separate investigation')
        matrix=np.zeros((n,n));matrix[self.sparse_rows,self.sparse_cols]=values
        x=np.linalg.solve(matrix,rhs);self.dense_fallbacks+=1
        if np.linalg.norm(rhs-product(x))>1e-10*scale:raise ValueError('independent linear true residual exceeded tolerance')
        return x

    def force(self, state):
        stress=state[:,:6].reshape(*self.weight.shape,6)
        force=-np.einsum('eqai,eqa,eq->ei',self.b,stress,self.weight,optimize=True)
        out=np.zeros(self.xyz.size);np.add.at(out,self.dofs.reshape(-1),force.reshape(-1));return out

    def rate(self,state,force_rate,plastic):
        self.calls+=1
        yf,d,n,dn,den=terms(state,self.parameters)
        if np.any(den[plastic]<=0):raise ValueError('nonpositive active MCC consistency denominator')
        ep=d.copy();ep[plastic]-=np.einsum('qa,qb,q->qab',dn[plastic],dn[plastic],1/den[plastic])
        ep=ep.reshape(*self.weight.shape,6,6)
        ke=np.einsum('eqai,eqab,eqbj,eq->eij',self.rb,ep,self.rb,self.weight,optimize=True)
        coordinates=np.zeros(self.xyz.size);coordinates[self.free]=self.linear_solve(ke,self.reduce(force_rate))
        du=self.lift(coordinates)
        de=-np.einsum('eqai,ei->eqa',self.b,du[self.dofs]).reshape(-1,6)
        # At the centre of an elastic MCC ellipse (p=pc/2, q=0) both
        # gradient and consistency denominator vanish. No plastic law is used
        # there. For inactive points only the loading sign is relevant.
        divisor=np.where(plastic,den,np.maximum(np.abs(den),1.))
        potential=np.einsum('qa,qa->q',dn,de)/divisor
        multiplier=potential*plastic
        flow=multiplier[:,None]*n
        dy=np.zeros_like(state);dy[:,:6]=np.einsum('qab,qb->qa',d,de-flow)
        dy[:,6]=np.exp(state[:,7])/(self.parameters[:,1]-self.parameters[:,2])*flow[:,:3].sum(axis=1)
        dy[:,7]=-de[:,:3].sum(axis=1);dy[:,8:]=flow
        return du,dy,potential*np.linalg.norm(n,axis=1),yf

    def select(self,state,force_rate,plastic):
        eligible=plastic|(terms(state,self.parameters)[0]>=-2e-9)
        plastic=plastic.copy();seen=set()
        for _ in range(80):
            key=plastic.tobytes()
            if key in seen:raise ValueError('independent continuum active set cycled; no displacement accepted')
            seen.add(key);_,_,potential,_=self.rate(state,force_rate,plastic)
            new=eligible&(potential>1e-12)
            if np.array_equal(new,plastic):return plastic
            plastic=new
        raise ValueError('independent continuum active set iteration limit')

    def advance(self,state,displacement,force_rate,count,progress=None):
        """Reuse only the identical step-start rate inside event localization."""
        from mcc_rate_midpoint import advance_cached
        return advance_cached(self,state,displacement,force_rate,count,progress)

    def advance_uncached(self,state,displacement,force_rate,count,progress=None):
        """Integrate load progress 0..1, localize switches; never project stresses."""
        t=0.;state=state.copy();u=displacement.copy();plastic=np.zeros(len(state),dtype=bool)
        maximum_drift=0.;events=0;accepted=0;minimum_step=1.;step=1/count
        def midpoint(s,h,mask):
            _,dy,_,_=self.rate(s,force_rate,mask)
            du,dy,_,_=self.rate(s+dy*h/2,force_rate,mask)
            return s+dy*h,du*h
        while t<1-1e-13:
            plastic=self.select(state,force_rate,plastic)
            h=min(step,1-t);candidate,du=midpoint(state,h,plastic)
            yf=terms(candidate,self.parameters)[0]
            _,_,potential,_=self.rate(candidate,force_rate,plastic)
            crossed=(~plastic)&(yf>1e-9);unloaded=plastic&(potential<-1e-11);next_plastic=plastic.copy()
            if np.any(crossed)|np.any(unloaded):
                lo=0.;hi=h
                for _ in range(28):
                    mid=(lo+hi)/2;trial,_=midpoint(state,mid,plastic)
                    fy=terms(trial,self.parameters)[0]
                    _,_,ld,_=self.rate(trial,force_rate,plastic)
                    if np.any((~plastic)&(fy>1e-9)) or np.any(plastic&(ld<-1e-11)):hi=mid
                    else:lo=mid
                h=hi;candidate,du=midpoint(state,h,plastic)
                fy=terms(candidate,self.parameters)[0];_,_,ld,_=self.rate(candidate,force_rate,plastic)
                next_plastic|=(~plastic)&(fy>=0.999e-9);next_plastic&=ld>=-0.999e-11
                events+=1
            # The switch may have shortened h; assess the accepted candidate.
            actual_yield=terms(candidate,self.parameters)[0]
            drift=float(np.max(np.maximum(actual_yield,0.)))
            if np.any(next_plastic):drift=max(drift,float(np.max(np.abs(actual_yield[next_plastic]))))
            if drift>2e-6:
                step=h/2
                if step<1e-9:raise ValueError('independent yield drift cannot be controlled without a projection')
                continue
            if h<1e-11:raise ValueError('independent event integration stalled')
            state=candidate;u+=du;plastic=next_plastic;t+=h;accepted+=1;minimum_step=min(minimum_step,h);maximum_drift=max(maximum_drift,drift)
            if progress is not None and accepted%100==0:progress(t,accepted,events)
        return state,u,dict(accepted_steps=accepted,events=events,maximum_normalized_yield_drift=maximum_drift,
                           minimum_progress_step=minimum_step,rate_assemblies=self.calls,
                           independent_pcg_iterations=self.linear_iterations,dense_fallbacks=self.dense_fallbacks)


def seed(mesh,ground,identity=None):
    out=[]
    cfg=ground['initialization'].get('geostatic')
    regions={r['id']:r['body'] for r in mesh['regions']};point_index=0
    if cfg and cfg['seed']['kind']=='integration_points':
        if identity is None or identity.get('schema')!='soil-fem-mesh-identity/1' or identity['mesh_sha256']!=cfg['seed']['mesh_sha256']:
            raise ValueError('integration-point seeds require matching mesh identity metadata')
    for ei in range(len(mesh['tetrahedra'])):
        if cfg:
            q=ref.volume(mesh,ei);rows=[]
            for qi in range(len(q['weight'])):
                if cfg['seed']['kind']=='body_states':row=cfg['seed']['states'][regions[mesh['region_ids'][ei]]]
                elif cfg['seed']['kind']=='integration_points':
                    source=cfg['seed']['points'][point_index];meta=identity['integration_points'][point_index]
                    if source['element_index']!=ei or source['local_index']!=qi or source['material_id']!=mesh['material_ids'][ei]:
                        raise ValueError('seed point association differs from the reference mesh')
                    if meta['element_index']!=ei or meta['local_index']!=qi or not np.isclose(meta['weight_m3'],q['weight'][qi],rtol=1e-10,atol=1e-12):
                        raise ValueError('mesh identity integration inventory differs from independent quadrature')
                    if not np.allclose(source['position_m'],q['position'][qi],rtol=0.,atol=1e-8) or not np.allclose(meta['position_m'],q['position'][qi],rtol=0.,atol=1e-8):
                        raise ValueError('seed point position differs from independent geometry')
                    row=source['state']
                else:raise ValueError('unsupported geostatic seed kind')
                if row['preconsolidation_pressure_kpa']<=0 or row['void_ratio']<=0:raise ValueError('inadmissible MCC seed pc or void ratio')
                rows.append([*row['effective_stress_compression_positive_kpa'],np.log(row['preconsolidation_pressure_kpa']),
                             np.log1p(row['void_ratio']),*row['plastic_strain_compression_positive']]);point_index+=1
            out.append(np.asarray(rows));continue
        q=ref.volume(mesh,ei);s,pore,material,layer=ref.initial_fields(mesh,ground,ei,q['position'])
        if np.any(pore):raise ValueError('this reference supports dry mechanics only')
        m=material['mcc'];mean=s[:,:3].mean(axis=1)
        if layer['consolidation']['definition']!='mcc_pc_over_mean_effective':raise ValueError('explicit mean-effective OCR required')
        if layer['void_ratio']['mode']!='from_mcc_reference':raise ValueError('horizontal reference requires explicit MCC reference void-ratio mode')
        pc=mean*layer['consolidation']['ratio']
        v=m['specific_volume_ncl_at_reference']-m['lambda']*np.log(pc/m['reference_pressure_kpa'])+m['kappa']*np.log(pc/mean)
        state=np.zeros((len(s),14));state[:,:6]=s;state[:,6]=np.log(pc);state[:,7]=np.log(v);out.append(state)
    if cfg and cfg['seed']['kind']=='integration_points' and (point_index!=len(cfg['seed']['points']) or point_index!=len(identity['integration_points'])):
        raise ValueError('seed or identity point count differs from the reference mesh')
    return np.stack(out)


def solve(mesh,ground,case,count,progress=None,identity=None):
    if ground['initialization']['method'] not in ('horizontal_k0','equilibrate') or ground['initialization']['water_level_m'] is not None:
        raise ValueError('independent reference requires dry K0 or explicit geostatic initialization')
    if any(case.get(k) for k in ('joint_system','beam_system','anchor_system','initially_inactive_bodies')):
        raise ValueError('unsupported reference boundary/member configuration')
    active=np.ones(len(mesh['tetrahedra']),dtype=bool);state=seed(mesh,ground,identity);u=np.zeros(np.asarray(mesh['nodes']).size)
    pressure=ground['initialization']['surface_surcharge_kpa'];original=case.get('pressure_on_original_surface_only',False)
    cfg=ground['initialization'].get('geostatic');factor=cfg['initial_boundary_factor'] if cfg else 0.
    def load():
        result=physical(mesh,ground,active,pressure,original)[0]
        result+=factor*ref.conditions_load(mesh,case.get('boundary') or {},active,skip_removed_surfaces=True)
        return result.reshape(-1)
    frames=[]
    def frame(name,system,diagnostics):
        local=state[active].reshape(-1,14);physical_load=load()
        residual=np.linalg.norm(system.reduce(physical_load-system.force(local)))/max(np.linalg.norm(physical_load),1.)
        frames.append(dict(name=name,displacement_m=u.reshape(-1,3).tolist(),active_elements=active.tolist(),
            material_state=state.tolist(),relative_free_force_residual=float(residual),diagnostics=diagnostics))
    system=System(mesh,ground,active,case)
    if np.max(terms(state.reshape(-1,14),system.parameters)[0])>1e-8:raise ValueError('initial seed lies outside MCC yield surface')
    geostatic=None
    if cfg:
        original_seed=state.copy();initial_force=system.force(state.reshape(-1,14))
        local,u,diagnostics=system.advance(state.reshape(-1,14),u,load()-initial_force,count,progress)
        state=local.reshape(state.shape)
        maximum_strain=float(np.max(np.abs(np.einsum('eqai,ei->eqa',system.b,u[system.dofs]))))
        if maximum_strain>cfg['maximum_initial_engineering_strain']:raise ValueError('independent initial strain exceeded the input limit')
        geostatic=dict(seed_material_state=original_seed.tolist(),maximum_engineering_strain=maximum_strain,diagnostics=diagnostics)
        frame('initial_equilibrated',system,diagnostics)
    else:frame('initial_k0',system,{})
    reference_displacement=u.reshape(-1,3).copy()
    if frames[-1]['relative_free_force_residual']>1e-8:raise ValueError('reference initial field is not in discrete equilibrium')
    for stage in case['stages']:
        before=load()
        if stage['kind']=='pressure':
            pressure=stage['target_pressure_kpa'];system=System(mesh,ground,active,case)
        elif stage['kind']=='boundary':
            factor=stage['target_factor'];system=System(mesh,ground,active,case)
        elif stage['kind']=='excavate':
            groups=[]
            for body in stage['bodies']:
                if body in mesh.get('excavation_sets',{}):
                    groups.append(mesh['excavation_sets'][body])
                else:
                    regions={r['id'] for r in mesh['regions'] if r['body']==body}
                    indices=[i for i,r in enumerate(mesh['region_ids']) if r in regions]
                    if not indices:raise ValueError('unknown or empty reference excavation body: '+body)
                    groups.append(indices)
            remove=np.unique(np.concatenate(groups))
            if not active[remove].all():raise ValueError('overlapping reference excavation')
            active[remove]=False;pressure=stage['remaining_surface_pressure_kpa'];system=System(mesh,ground,active,case)
            before=system.force(state[active].reshape(-1,14))
        else:raise ValueError('reference supports pressure/boundary/excavation stages only')
        after=load()
        local,u,diag=system.advance(state[active].reshape(-1,14),u,after-before,count,progress)
        state[active]=local.reshape(state[active].shape);frame(stage['name'],system,diag)
    return dict(schema='soil-fem-independent-mcc-rate-fem/1',frames=frames,geostatic=geostatic,reference_displacement_m=reference_displacement.tolist(),
                method='analytic continuum tangent; independent sparse NumPy PCG/global assembly and nodal SVD; explicit midpoint with loading/unloading events; no stress projection',
                limitations='dry small-strain MCC; homogeneous directional constraints; initial/pressure/boundary/excavation stages; acceptance requires independent refinement and field comparison')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mesh',type=Path);p.add_argument('ground',type=Path);p.add_argument('case',type=Path)
    p.add_argument('output',type=Path);p.add_argument('--steps',type=int,default=200);p.add_argument('--identity',type=Path);a=p.parse_args()
    if a.output.exists():raise ValueError('refusing to overwrite independent reference output')
    if a.steps<2:raise ValueError('at least two reference steps required')
    inputs=[a.mesh,a.ground,a.case,Path(__file__),Path(__file__).with_name('mcc_rate_midpoint.py'),
            Path(ref.__file__),Path(__file__).with_name('verify_high_order.py')]+([a.identity] if a.identity else [])
    fingerprints={str(path):hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs}
    report=solve(*(json.loads(path.read_text()) for path in (a.mesh,a.ground,a.case)),a.steps,
        progress=lambda t,n,e:print(f'Independent progress {t:.6f}, steps {n}, events {e}',flush=True),identity=json.loads(a.identity.read_text()) if a.identity else None)
    if fingerprints!={str(path):hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs}:
        raise ValueError('reference inputs or source files changed during the computation')
    report['input_sha256']=fingerprints
    report['reference_steps']=a.steps
    report['integration_implementation']='cached step-start midpoint; unchanged events and yield-drift limits'
    with a.output.open('x') as f:json.dump(report,f)
    print('Independent MCC 3D rate FEM completed',flush=True)
