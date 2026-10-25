"""Independent transfer/force audit and elastic 3D PDE given hydraulic heads.

The mechanical solution uses only inputs and hydraulic primary heads, never
native displacement/stress updates. Hydraulic DOFs are an input to this
ONE_WAY comparison; an independent coupled flow solution is not implied.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import copy
import numpy as np
import fem_reference as fem
from richards_column_reference import retention
from hydrostatic_reference import face_force
from verify_high_order import physical
from verify_richards import audit as hydraulic_audit


def geometry(mesh,e,independent=False,order=5):
    ids=fem.cells(mesh)[e];xyz=np.asarray(mesh['nodes'])[ids]
    r,w=fem.volume_rule({'tetrahedra10':[True],'quadrature_order':mesh.get('quadrature_order',3)},independent,order)
    shape,gradient=fem.basis(r,len(ids)==10);jac=np.einsum('ai,qaj->qij',xyz-xyz[0],gradient);det=np.linalg.det(jac);assert np.all(det>0)
    grad=np.einsum('qai,qij->qaj',gradient,np.linalg.inv(jac))
    return dict(nodes=ids,position=xyz[0]+shape@(xyz-xyz[0]),weight=w*det,shape=shape,gradient=grad)


def constraint_space(mesh,case,ground,active=None):
    xyz=np.asarray(mesh['nodes']);size=xyz.size;kind=case['settings']['kinematics'];offset=np.zeros(size)
    live=np.ones(len(xyz),dtype=bool)
    if active is not None:
        live[:]=False;live[np.unique(fem.cells(mesh)[np.asarray(active,dtype=bool)])]=True
    if kind=='oedometer':
        planes=sorted(set(np.round(xyz[live,2],8)));q=np.zeros((size,len(planes)-1))
        for n,p in enumerate(xyz):
            if not live[n]:continue
            j=planes.index(round(p[2],8))
            if j:q[3*n+2,j-1]=1.
        return q,offset
    assert kind in ('side_rollers','explicit')
    factor=ground['initialization'].get('geostatic',{}).get('initial_boundary_factor',0.)
    blocks=[];count=0
    for n,p in enumerate(xyz):
        if not live[n]:
            blocks.append((n,count,np.zeros((3,0))));continue
        rows=[];values=[]
        if kind=='side_rollers':
            for a in range(3):
                fixed=abs(p[2]-xyz[:,2].min())<1e-8 if a==2 else min(abs(p[a]-xyz[:,a].min()),abs(p[a]-xyz[:,a].max()))<1e-8
                if fixed:rows.append(np.eye(3)[a]);values.append(0.)
        for c in case.get('boundary',{}).get('constraints',[]):
            if n not in mesh['node_sets'][c['node_set']]:continue
            d=np.asarray(c['direction']);norm=np.linalg.norm(d);rows.append(d/norm)
            values.append(factor*(c['value_m']+np.asarray(c.get('gradient',[0.,0.,0.]))@(p-c.get('reference_point_m',[0.,0.,0.])))/norm)
        if rows:
            a=np.asarray(rows);b=np.asarray(values);_,s,vh=np.linalg.svd(a,full_matrices=True);rank=int((s>1e-10*max(s.max(),1)).sum());block=vh[rank:].T
            offset[3*n:3*n+3]=np.linalg.lstsq(a,b,rcond=1e-10)[0];assert np.max(np.abs(a@offset[3*n:3*n+3]-b))<1e-10
        else:block=np.eye(3)
        blocks.append((n,count,block));count+=block.shape[1]
    q=np.zeros((size,count))
    for n,j,block in blocks:q[3*n:3*n+3,j:j+block.shape[1]]=block
    return q,offset


def base_load(mesh,ground,mechanics,active=None):
    active=np.ones(len(mesh['tetrahedra']),dtype=bool) if active is None else np.asarray(active,dtype=bool)
    out,_=physical(mesh,ground,active,ground['initialization']['surface_surcharge_kpa'],mechanics.get('pressure_on_original_surface_only',False))
    factor=ground['initialization'].get('geostatic',{}).get('initial_boundary_factor',0.)
    if factor:out+=factor*fem.conditions_load(mesh,mechanics.get('boundary',{}),active)
    return out


def fluid_load(mesh,ground,case,hydraulic,reservoirs,geos,base):
    law={m['id']:m['retention'] for m in case['flow']['materials']};gamma=ground['initialization']['water_unit_weight_kn_m3'];datum=case['flow']['head_datum_m']
    internal=np.zeros_like(base);body=np.zeros_like(base);exterior=np.zeros_like(base);mass=0.;raw=[]
    for e,g in enumerate(geos):
        if hydraulic['relative_cell_head_m'][e] is None:raw.append(None);continue
        psi=(datum-g['position'][:,2])+hydraulic['relative_cell_head_m'][e];p=law[mesh['material_ids'][e]];theta,_,_=retention(psi,p)
        se=(theta-p['residual_water_content'])/(p['saturated_water_content']-p['residual_water_content']);sw=theta/p['saturated_water_content']
        closure=case['bishop'];chi=se if closure['kind']=='effective_saturation' else sw if closure['kind']=='water_saturation' else np.full(len(psi),closure['value'])
        u=gamma*psi;effective=chi*u;water=theta+p['saturated_specific_storage_per_m']*np.maximum(psi,0)
        internal[g['nodes']]-=np.einsum('qai,q,q->ai',g['gradient'],effective,g['weight'])
        body[g['nodes'],2]-=gamma*np.einsum('qa,q,q->a',g['shape'],water,g['weight']);mass+=g['weight']@water
        raw.append(dict(pressure_head_m=psi,physical_pore_pressure_kpa=u,bishop_chi=chi,effective_pore_stress_kpa=effective,water_content=theta,effective_saturation=se,water_saturation=sw,water_equivalent_mass_per_volume=water))
    owners=fem.topology(mesh)
    for reservoir in reservoirs:
        for face in mesh['face_sets'][reservoir['face_set']]:
            live=[e for e in owners[tuple(sorted(face))] if hydraulic['relative_cell_head_m'][e] is not None];assert len(live)==1
            e=live[0];nodes=fem.face_nodes(mesh,e,face);force,_=face_force(np.asarray(mesh['nodes'])[nodes],reservoir['level_m'],gamma)
            exterior[nodes]+=force
    actual=base+body+exterior
    return dict(fluid_internal_nodal_force_kn=internal.ravel(),water_body_nodal_load_kn=body.ravel(),exterior_nodal_load_kn=exterior.ravel(),physical_nodal_load_kn=actual.ravel(),equivalent_mechanical_rhs_kn=(actual-internal).ravel(),water_equivalent_mass_volume_m3=float(mass),raw=raw)


def fluid_reservoirs(run,frame):
    if frame['phase']!='transient':return run['case']['initial_exterior_water']
    return next(p['reservoirs'] for p in run['case']['exterior_periods'] if p['name']==frame['name'])


def reconstructed_samples(mesh,run,hydraulic):
    """Derived output only; the primary material history/heads remain native."""
    hyd=[geometry(mesh,e) for e in range(len(mesh['tetrahedra']))];mechanical=[]
    for e in range(len(hyd)):
        g=fem.volume(mesh,e);g['gradient']=np.stack([g['b'][:,0,0::3],g['b'][:,1,1::3],g['b'][:,2,2::3]],axis=-1);mechanical.append(g)
    result=[]
    for geos in (hyd,mechanical):
        values=fluid_load(mesh,run['ground'],run['case'],hydraulic,[],geos,np.zeros_like(mesh['nodes'],dtype=float))['raw'];points=[]
        for e,(g,data) in enumerate(zip(geos,values)):
            if data is None:continue
            for i in range(len(g['weight'])):points.append(dict(element_index=e,local_index=i,position_m=g['position'][i].tolist(),weight_m3=float(g['weight'][i]),**{k:float(v[i]) for k,v in data.items()}))
        result.append(points)
    return result


def hydraulic_view(mesh,run):
    first=dict(name='initial',state=run['seed']['hydraulic'],points=[],balance=None)
    frames=[first]+[dict(name=f['name'],state=f['state']['hydraulic'],points=[],balance=f['hydraulic_balance']) for f in run['frames'] if f['phase']=='transient']
    attempts=[dict(period=a['period'],accepted=True,time_step_s=a['time_step_s'],start_time_s=a['start_time_s']) for a in run['attempts'] if a['phase']=='transient' and a['accepted']]
    return dict(schema='soil-fem-richards-result/1',complete=True,active_elements=[True]*len(mesh['tetrahedra']),cell_element_indices=list(range(len(mesh['tetrahedra']))),material_ids=mesh['material_ids'],case=run['case']['flow'],faces=run['hydraulic_faces'],frames=frames,attempts=attempts)


def audit(mesh,run):
    assert run['complete'] and run['schema']=='soil-fem-one-way-result/1'
    ground=run['ground'];case=run['case'];geos=[geometry(mesh,e) for e in range(len(mesh['tetrahedra']))];mechanical=[fem.volume(mesh,e) for e in range(len(mesh['tetrahedra']))]
    base=base_load(mesh,ground,case['mechanics']);q,offset=constraint_space(mesh,case['mechanics'],ground);maximum_force_error=maximum_field_error=maximum_residual=0.;rows=[]
    for f in run['frames']:
        state=f['state'];fluid=fluid_load(mesh,ground,case,state['hydraulic'],fluid_reservoirs(run,f),geos,base);native=state['fluid']
        for name in ('fluid_internal_nodal_force_kn','water_body_nodal_load_kn','exterior_nodal_load_kn','physical_nodal_load_kn','equivalent_mechanical_rhs_kn'):
            error=float(np.max(np.abs(fluid[name]-native[name])));maximum_force_error=max(maximum_force_error,error)
            assert error<1e-8*max(np.linalg.norm(fluid[name]),1.),(f['name'],name,error)
        assert abs(fluid['water_equivalent_mass_volume_m3']-native['water_equivalent_mass_volume_m3'])<1e-11*max(fluid['water_equivalent_mass_volume_m3'],1.)
        for key in (() if not native['hydraulic_points'] else ('physical_pore_pressure_kpa','bishop_chi','effective_pore_stress_kpa','water_content','effective_saturation','water_saturation','water_equivalent_mass_per_volume')):
            expected=np.concatenate([a[key] for a in fluid['raw']]);stored=np.array([p[key] for p in native['hydraulic_points']]);error=float(np.max(np.abs(expected-stored)));maximum_field_error=max(maximum_field_error,error)
            assert np.allclose(expected,stored,rtol=1e-10,atol=1e-12),(f['name'],key,error)
        if f['initialization_progress']<1.:continue
        # Material stresses are audited, not regenerated by a second plastic
        # return map. The independent elastic PDE below is a separate check.
        internal=np.zeros_like(base);start=0
        for g in mechanical:
            count=len(g['weight']);states=state['mechanical']['material_states'][start:start+count];start+=count
            sigma=np.array([s['effective_stress_compression_positive_kpa'] for s in states]);internal[g['nodes']]-=np.einsum('qai,qa,q->i',g['b'],sigma,g['weight']).reshape(-1,3)
        residual=q.T@(fluid['equivalent_mechanical_rhs_kn']-internal.ravel());relative=float(np.linalg.norm(residual)/max(np.linalg.norm(fluid['equivalent_mechanical_rhs_kn']),1.));maximum_residual=max(maximum_residual,relative)
        assert relative<1.05*case['mechanics']['settings']['relative_force_tolerance']+1e-10,(f['name'],relative)
        displacement=np.asarray(state['mechanical']['displacement_m']).ravel();free=np.linalg.lstsq(q,displacement-offset,rcond=None)[0]
        assert np.max(np.abs(q@free+offset-displacement))<1e-10
        rows.append(dict(name=f['name'],time_s=state['hydraulic']['time_s'],independent_relative_mechanical_force_residual=relative,water_mass_volume_m3=fluid['water_equivalent_mass_volume_m3']))
    flow=hydraulic_audit(mesh,hydraulic_view(mesh,run))
    return dict(schema='soil-fem-one-way-independent-force-audit/1',passed=True,maximum_independent_force_component_error_kn=maximum_force_error,maximum_raw_field_difference=maximum_field_error,maximum_relative_mechanical_force_residual=maximum_residual,hydraulic_audit=flow,frames=rows,
        scope='Independent constitutive/weight/traction/force reassembly. Stored material stresses are audited; this is not an independent plastic mechanics or coupled PDE solution.')


class ElasticReference:
    def __init__(self,mesh,ground,case,order=5,active=None):
        assert all('mohr_coulomb' in m for m in ground['materials'])
        self.active=np.ones(len(mesh['tetrahedra']),dtype=bool) if active is None else np.asarray(active,dtype=bool)
        self.mesh=mesh;self.ground=ground;self.case=case;self.q,self.offset=constraint_space(mesh,case['mechanics'],ground,active)
        size=len(mesh['nodes'])*3;self.base=base_load(mesh,ground,case['mechanics'],active);self.geos=[geometry(mesh,e,True,order) for e in range(len(mesh['tetrahedra']))]
        self.internal=np.zeros(size);full=np.zeros((size,size));materials={m['id']:m for m in ground['materials']};regions={r['id']:r for r in mesh['regions']}
        for e in range(len(mesh['tetrahedra'])):
            if not self.active[e]:continue
            g=fem.volume(mesh,e,True,order);ids=(g['nodes'][:,None]*3+np.arange(3)).ravel();mat=materials[mesh['material_ids'][e]]['mohr_coulomb'];d=fem.elastic_tensor(mat['young_kpa'],mat['poisson'])
            if ground['initialization']['method']=='equilibrate':
                cfg=ground['initialization']['geostatic'];assert cfg['seed']['kind']=='body_states';body=regions[mesh['region_ids'][e]]['body'];s=cfg['seed']['states'][body]['effective_stress_compression_positive_kpa'];sigma=np.broadcast_to(s,(len(g['weight']),6))
            else:sigma=fem.initial_fields(mesh,ground,e,g['position'])[0]
            self.internal[ids]-=np.einsum('qai,qa,q->i',g['b'],sigma,g['weight'])
            full[np.ix_(ids,ids)]+=np.einsum('qai,ab,qbj,q->ij',g['b'],d,g['b'],g['weight'],optimize=True)
        self.full=full;self.k=self.q.T@full@self.q
        # The Cholesky is independent of native tangents and material updates.
        self.factor=np.linalg.cholesky(self.k)

    def solve(self,hydraulic,reservoirs,holding=None):
        fluid=fluid_load(self.mesh,self.ground,self.case,hydraulic,reservoirs,self.geos,self.base)
        load=fluid['equivalent_mechanical_rhs_kn']+(np.zeros_like(self.internal) if holding is None else np.asarray(holding))
        rhs=self.q.T@(load-self.internal-self.full@self.offset)
        # NumPy direct solve on the original reduced operator.
        free=np.linalg.solve(self.k,rhs);residual=np.linalg.norm(self.k@free-rhs)/max(np.linalg.norm(rhs),1.);assert residual<1e-10
        return dict(displacement_m=(self.offset+self.q@free).reshape(-1,3),relative_linear_residual=float(residual))
