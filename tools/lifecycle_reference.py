"""Independent incremental elastic FE solution and lifecycle state-transfer audit.

The reference owns its stresses and displacement. Hydraulic primary heads are
inputs, so this is neither an independent Richards solve nor plastic 3D proof.
"""
import copy
import numpy as np
import fem_reference as fem
import one_way_reference as ref
from verify_high_order import internal, physical
from verify_richards import audit as flow_audit
from verify_one_way_construction import domain


def ids(mesh,state):
    return state['mechanical'].get('element_material_ids') or mesh['material_ids']


def load(mesh,ground,case,active,material_ids):
    out,_=physical(mesh,ground,active,ground['initialization']['surface_surcharge_kpa'],case['mechanics'].get('pressure_on_original_surface_only',False),material_ids=material_ids)
    factor=ground['initialization'].get('geostatic',{}).get('initial_boundary_factor',0.)
    if factor:out+=factor*fem.conditions_load(mesh,case['mechanics'].get('boundary',{}),active)
    return out


def reservoirs(run,frame):
    case=run['case'];name=frame['name']
    if frame['phase']!='transient':
        name=next(e['before_period'] for e in case.get('excavations',[])+case.get('transitions',[]) if e['name']==name)
    return next(p['reservoirs'] for p in case['one_way']['exterior_periods'] if p['name']==name)


def audit(mesh,run):
    assert run['complete'];ground=run['ground'];case=run['case']['one_way']
    geos=[ref.geometry(mesh,e) for e in range(len(mesh['tetrahedra']))]
    mechanical=[fem.volume(mesh,e) for e in range(len(geos))];offset=np.cumsum([0]+[len(g['weight']) for g in mechanical])
    previous=run['initialization']['last_committed'];active=previous['mechanical']['active_elements'];material_ids=ids(mesh,previous)
    segments=[];segment=[dict(name='initial',state=previous['hydraulic'],balance=None,points=[])];faces=None
    maximum_force=maximum_residual=maximum_budget=maximum_transfer=0.;rows=[];seen=set()
    events={r['name']:r for r in run.get('transitions',[])};cuts={r['name']:r for r in run['excavations']}
    def mass(e,h,mid):
        p=next(m['retention'] for m in case['flow']['materials'] if m['id']==mid);g=geos[e]
        psi=case['flow']['head_datum_m']-g['position'][:,2]+h;theta,_,_=ref.retention(psi,p)
        return float(g['weight']@(theta+p['saturated_specific_storage_per_m']*np.maximum(psi,0.)))
    for f in run['frames']:
        state=f['state'];now_active=state['mechanical']['active_elements'];now_ids=ids(mesh,state)
        jump=f['phase']!='transient' and f['name'] not in seen
        if jump:
            seen.add(f['name']);segments.append((active,material_ids,faces,segment))
            r=(events if f['phase']=='material_transition' else cuts)[f['name']]
            water=r['water'] if f['phase']=='material_transition' else r['removal'];origin=previous['hydraulic']
            if f['phase']=='material_transition':
                selected={c['element_index'] for c in water['changed_cells']}
                for c in water['changed_cells']:
                    e=c['element_index'];after=c['after'];new=mass(e,after['relative_head_m'],c['new_material_id']);maximum_transfer=max(maximum_transfer,abs(new-after['water_equivalent_volume_m3']))
                    assert abs(new-after['water_equivalent_volume_m3'])<1e-11*max(new,1.)
                    if c['before'] is not None:
                        old=mass(e,c['before']['relative_head_m'],c['old_material_id']);assert abs(old-c['before']['water_equivalent_volume_m3'])<1e-11*max(old,1.)
                action=next(t['action'] for t in run['case']['transitions'] if t['name']==f['name'])
                if action['kind']=='material_set' and action['water']['kind']=='preserve_water':
                    assert water['water_added_m3']==water['water_removed_m3']==0.
                    for c in water['changed_cells']:
                        assert abs(c['before']['water_equivalent_volume_m3']-c['after']['water_equivalent_volume_m3'])<1e-12
                archived=r['old_material_points']
            else:
                selected={c['element_index'] for c in water['removed_cells']};total=0.
                for c in water['removed_cells']:
                    e=c['element_index'];actual=mass(e,c['relative_head_m'],material_ids[e]);assert abs(actual-c['water_equivalent_volume_m3'])<1e-11*max(actual,1.);total+=actual
                assert abs(total-water['water_removed_m3'])<1e-11*max(total,1.)
                archived=r.get('removed_material_points',[])
            for p in archived:
                j=p['integration_point_index'];assert p['state']==previous['mechanical']['material_states'][j]
                assert p['material_id']==material_ids[run['mechanical_integration_points'][j]['element_index']]
            for e,live in enumerate(active):
                if live and e not in selected:
                    for key in ('relative_cell_head_m','cell_storage_offset_m3','cell_outflow_m3_s'):assert water['state'][key][e]==origin[key][e]
            for i,j in enumerate(water['new_to_old_face']):
                if j is not None:
                    assert water['new_faces'][i]['corners']==water['old_faces'][j]['corners']
                    assert water['state']['relative_face_head_m'][i]==origin['relative_face_head_m'][j]
            segment=[dict(name='boundary_jump',state=water['state'],balance=None,points=[])];faces=f['hydraulic_faces']
        active=now_active;material_ids=now_ids;faces=f['hydraulic_faces']
        if f['phase']=='transient':segment.append(dict(name=f['name'],state=state['hydraulic'],balance=f['hydraulic_balance'],points=[]))
        current=domain(mesh,active);current['material_ids']=material_ids
        fluid=ref.fluid_load(current,ground,case,state['hydraulic'],reservoirs(run,f),geos,load(mesh,ground,case,active,material_ids))
        for key in ('fluid_internal_nodal_force_kn','water_body_nodal_load_kn','exterior_nodal_load_kn','physical_nodal_load_kn','equivalent_mechanical_rhs_kn'):
            error=float(np.max(np.abs(fluid[key]-state['fluid'][key])));maximum_force=max(maximum_force,error);assert error<1e-8*max(np.linalg.norm(fluid[key]),1.)
        force,_=internal(mesh,ground,state['mechanical']['material_states'],active,mechanical,offset,material_ids=material_ids,pore_pressure_kpa=[0.]*int(offset[-1]))
        q,affine=ref.constraint_space(mesh,case['mechanics'],ground,active);rhs=fluid['equivalent_mechanical_rhs_kn']+np.asarray(f['holding_load_kn'])
        residual=float(np.linalg.norm(q.T@(rhs-force.ravel()))/max(np.linalg.norm(rhs),1.));maximum_residual=max(maximum_residual,residual)
        assert residual<1.05*case['mechanics']['settings']['relative_force_tolerance']+1e-10
        error=fluid['water_equivalent_mass_volume_m3']+f['cumulative_removed_water_m3']-f.get('cumulative_added_water_m3',0.)-run['initial_water_m3']-f['cumulative_net_water_inflow_m3']
        maximum_budget=max(maximum_budget,abs(error));assert abs(error-f['signed_global_water_budget_error_m3'])<1e-10*max(run['initial_water_m3'],1.)
        rows.append(dict(name=f['name'],phase=f['phase'],relative_force_residual=residual,water_budget_error_m3=error));previous=state
    segments.append((active,material_ids,faces,segment));hydraulic=[]
    for active,material_ids,faces,frames in segments:
        if len(frames)<2:continue
        pairs={(f['name'],f['state']['time_s']) for f in frames[1:]}
        attempts=[dict(period=a['name'],accepted=True,time_step_s=a['dt_s'],start_time_s=a['time_s']) for a in run['attempts'] if a['accepted'] and a['phase']=='transient' and (a['name'],a['time_s']+a['dt_s']) in pairs]
        water_run=dict(schema='soil-fem-richards-result/1',complete=True,active_elements=active,cell_element_indices=np.flatnonzero(active).tolist(),material_ids=material_ids,case=case['flow'],faces=faces,frames=frames,attempts=attempts)
        hydraulic.append(flow_audit(domain(mesh,active),water_run))
    return dict(passed=True,maximum_force_component_difference_kn=maximum_force,maximum_relative_mechanical_residual=maximum_residual,maximum_global_water_budget_error_m3=maximum_budget,maximum_transfer_water_difference_m3=maximum_transfer,hydraulic_segments=hydraulic,frames=rows,scope='Independent water-law integration, unchanged primary/trace transfer, archived histories, stored-stress force and RT0 balance reassembly; not an independent plastic or flow PDE.')


class IncrementalElastic:
    def __init__(self,mesh,run,order):
        self.mesh=copy.deepcopy(mesh);self.ground=run['ground'];self.case=run['case']['one_way'];self.materials={m['id']:m for m in self.ground['materials']}
        assert all('mohr_coulomb' in m for m in self.materials.values())
        inactive=set(run['case'].get('initially_inactive_bodies',[]));regions={r['id']:r for r in mesh['regions']}
        self.bodies=[regions[r]['body'] for r in mesh['region_ids']];self.active=np.asarray([b not in inactive for b in self.bodies]);self.material_ids=mesh['material_ids'].copy()
        self.geo=[fem.volume(mesh,e,True,order) for e in range(len(self.active))];self.water_geo=[ref.geometry(mesh,e,True,order) for e in range(len(self.active))]
        self.stress=[]
        for e,g in enumerate(self.geo):
            if self.ground['initialization']['method']=='horizontal_k0':sigma=fem.initial_fields(mesh,self.ground,e,g['position'])[0]
            else:sigma=np.broadcast_to(self.ground['initialization']['geostatic']['seed']['states'][self.bodies[e]]['effective_stress_compression_positive_kpa'],(len(g['weight']),6)).copy()
            self.stress.append(sigma)
        self.u=np.zeros_like(mesh['nodes'],dtype=float)

    def change(self,event):
        selected=np.asarray([b in event['bodies'] for b in self.bodies]);kind=event['kind']
        if kind=='cut':self.active[selected]=False
        elif kind=='birth':
            old_nodes=np.unique(fem.cells(self.mesh)[self.active]);self.active[selected]=True
            if event['placement']=='cad_new_nodes':
                new_nodes=np.setdiff1d(np.unique(fem.cells(self.mesh)[self.active]),old_nodes);self.u[new_nodes]=0.
            for e in np.flatnonzero(selected):self.stress[e][:]=event['states'][self.bodies[e]]['effective_stress_compression_positive_kpa']
        else:
            assert kind=='material_set'
            for e in np.flatnonzero(selected):
                self.material_ids[e]=event['material_id']
                if event['history']['mode']=='replace':self.stress[e][:]=event['history']['states'][self.bodies[e]]['effective_stress_compression_positive_kpa']

    def solve(self,hydraulic,reservoir):
        q,offset=ref.constraint_space(self.mesh,self.case['mechanics'],self.ground,self.active);assert not np.any(offset)
        size=self.u.size;full=np.zeros((size,size));force=np.zeros(size);blocks=[]
        current=domain(self.mesh,self.active);current['material_ids']=self.material_ids
        fluid=ref.fluid_load(current,self.ground,self.case,hydraulic,reservoir,self.water_geo,load(self.mesh,self.ground,self.case,self.active,self.material_ids))
        for e in np.flatnonzero(self.active):
            g=self.geo[e];dofs=(g['nodes'][:,None]*3+np.arange(3)).ravel();p=self.materials[self.material_ids[e]]['mohr_coulomb'];d=fem.elastic_tensor(p['young_kpa'],p['poisson']);blocks.append((e,d,dofs))
            full[np.ix_(dofs,dofs)]+=np.einsum('qai,ab,qbj,q->ij',g['b'],d,g['b'],g['weight'],optimize=True)
            force[dofs]-=np.einsum('qai,qa,q->i',g['b'],self.stress[e],g['weight'])
        k=q.T@full@q;rhs=q.T@(fluid['equivalent_mechanical_rhs_kn']-force);increment=q@np.linalg.solve(k,rhs);self.u+=increment.reshape(-1,3)
        for e,d,dofs in blocks:self.stress[e]-=np.einsum('ab,qbi,i->qa',d,self.geo[e]['b'],increment[dofs])
        assert np.linalg.norm(k@np.linalg.lstsq(q,increment,rcond=None)[0]-rhs)<1e-9*max(np.linalg.norm(rhs),1.)
        return self.u.copy()


def elastic(mesh,run):
    assert all(not any(s['plastic_strain_compression_positive']) for f in run['frames'] for s in f['state']['mechanical']['material_states'])
    initial=run['initialization']['last_committed'];events={e['name']:dict(e,kind='cut') for e in run['case'].get('excavations',[])}
    events.update({e['name']:e['action'] for e in run['case'].get('transitions',[])})
    references=[IncrementalElastic(mesh,run,n) for n in (5,7,9)]
    old=[r.solve(initial['hydraulic'],run['case']['one_way']['initial_exterior_water']) for r in references];native_old=np.asarray(initial['mechanical']['displacement_m']);rows=[];seen=set()
    for f in run['frames']:
        if f['phase']!='transient' and f['progress']<1.:continue
        if f['phase']!='transient' and f['name'] not in seen:
            seen.add(f['name']);event=events[f['name']]
            old_active=references[0].active.copy()
            for r in references:r.change(event)
            if event['kind']=='birth' and event['placement']=='cad_new_nodes':
                new_nodes=np.setdiff1d(np.unique(fem.cells(mesh)[references[0].active]),np.unique(fem.cells(mesh)[old_active]));native_old=native_old.copy();native_old[new_nodes]=0.
                for v in old:v[new_nodes]=0.
        solved=[r.solve(f['state']['hydraulic'],reservoirs(run,f)) for r in references];nodes=np.unique(fem.cells(mesh)[references[0].active]);native=np.asarray(f['state']['mechanical']['displacement_m']);delta=[(v-o)[nodes] for v,o in zip(solved,old)]
        total=float(np.linalg.norm(native[nodes]-solved[-1][nodes])/max(np.linalg.norm(solved[-1][nodes]),1e-15));scale=np.linalg.norm(delta[-1]);native_delta=(native-native_old)[nodes]
        relative=float(np.linalg.norm(native_delta-delta[-1])/max(scale,1e-15));coarse=float(np.linalg.norm(delta[0]-delta[1])/max(scale,1e-15));fine=float(np.linalg.norm(delta[1]-delta[2])/max(scale,1e-15))
        assert total<3e-4,(f['name'],total)
        if scale>1e-8:assert relative<5e-4 and fine<1e-6,(f['name'],relative,coarse,fine)
        rows.append(dict(name=f['name'],phase=f['phase'],increment_norm_m=float(scale),total_relative_displacement_l2_difference=total,increment_relative_displacement_l2_difference=relative,independent_increment_quadrature_5_to_7_difference=coarse,independent_increment_quadrature_7_to_9_difference=fine,maximum_component_difference_m=float(np.max(np.abs(native[nodes]-solved[-1][nodes])))))
        old=solved;native_old=native
    return dict(passed=True,comparisons=rows,scope='Independent incremental elastic 3D FE using input K0 stresses, owned displacement/stresses and birth placement, current stiffness and preserved/replaced history. Given native hydraulic heads. No native mechanics generates the reference; no independent plastic or flow PDE is claimed.')
