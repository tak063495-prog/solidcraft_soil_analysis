"""Independent NumPy scalar FEM and ONE_WAY pressure audit (bounded dense).

The reference assembles from input geometry/conductivity/BCs, never Rust heads.
Nodal reaction outflow is -residual; raw gradient flux is a separate quantity.
"""
import numpy as np
import fem_reference as ref


def solve(mesh,case,active,material_ids=None):
    cells=ref.cells(mesh);xyz=np.asarray(mesh['nodes']);active=np.asarray(active,dtype=bool)
    used=np.unique(cells[active]);lookup=np.full(len(xyz),-1);lookup[used]=np.arange(len(used))
    if len(used)>2000:raise ValueError('independent dense Darcy reference limited to 2000 active heads')
    materials={m['id']:m for m in case['materials']};ids=material_ids or mesh['material_ids']
    matrix=np.zeros((len(used),len(used)));rhs=np.zeros(len(used));source=0.;neumann=0.;absolute_source=0.;absolute_neumann=0.
    quadrature={}
    for ei in np.flatnonzero(active):
        q=ref.volume(mesh,ei);quadrature[ei]=q;nodes=lookup[q['nodes']]
        g=np.stack([q['b'][:,a,a::3] for a in range(3)],axis=2)
        k=np.asarray(materials[ids[ei]]['hydraulic_conductivity_m_s']);k=(k+k.T)/2
        local=np.einsum('qia,ab,qjb,q->ij',g,k,g,q['weight'])
        matrix[np.ix_(nodes,nodes)]+=local
        rate=materials[ids[ei]].get('volume_source_per_s',0.)
        rhs[nodes]+=rate*(q['weight']@q['shape']);source+=rate*q['weight'].sum()
        absolute_source+=abs(rate)*q['weight'].sum()
    owners=ref.topology(mesh);heads={}
    for bc in case['heads']:
        for face in mesh['face_sets'][bc['face_set']]:
            adjacent=[e for e in owners[tuple(sorted(face))] if active[e]];assert len(adjacent)==1
            nodes=ref.face_nodes(mesh,adjacent[0],face)
            for n in nodes:
                offset=np.asarray(bc.get('gradient',[0.,0.,0.]))@(xyz[n]-bc.get('reference_point_m',[0.,0.,0.]))
                value=(bc['head_m'],offset)
                if n in heads:assert abs((heads[n][0]-value[0])+(heads[n][1]-value[1]))<1e-10
                heads[n]=value
    for bc in case.get('fluxes',[]):
        for face in mesh['face_sets'][bc['face_set']]:
            adjacent=[e for e in owners[tuple(sorted(face))] if active[e]];assert len(adjacent)==1
            q=ref.face(mesh,adjacent[0],face)
            if bc['kind']=='outward':out=bc['flux_m_s']*np.linalg.norm(q['area_vector'],axis=1)
            else:
                flux=np.asarray(bc['flux_m_s'])+(q['position']-bc.get('reference_point_m',[0.,0.,0.]))@np.asarray(bc.get('gradient_per_s',np.zeros((3,3)))).T
                out=np.einsum('qa,qa->q',flux,q['area_vector'])
            rhs[lookup[q['nodes']]]-=q['shape'].T@out;neumann+=out.sum()
            absolute_neumann+=abs(out).sum()
    head=np.zeros(len(used));datum=np.zeros(len(used));unseen=set(range(len(used)))
    while unseen:
        seed=min(unseen);unseen.remove(seed);group=[seed];queue=[seed]
        while queue:
            row=queue.pop()
            for node in np.flatnonzero(matrix[row]!=0):
                if node in unseen:unseen.remove(node);group.append(node);queue.append(node)
        anchors=[heads[used[node]][0] for node in sorted(group) if used[node] in heads]
        assert anchors,'independent domain is unanchored';datum[group]=anchors[0]
    fixed=np.asarray(sorted(lookup[list(heads)]));free=np.setdiff1d(np.arange(len(used)),fixed)
    for n,(base,offset) in heads.items():head[lookup[n]]=(base-datum[lookup[n]])+offset
    if len(free):head[free]=np.linalg.solve(matrix[np.ix_(free,free)],(rhs-matrix@head)[free])
    residual=matrix@head-rhs
    field=np.zeros(len(xyz));field[used]=head+datum
    relative=np.zeros(len(xyz));relative[used]=head;node_datum=np.zeros(len(xyz));node_datum[used]=datum
    reactions=np.zeros(len(xyz));reactions[used[fixed]]=-residual[fixed]
    den=absolute_source+absolute_neumann+abs(reactions).sum()
    balance=abs(reactions.sum()+neumann-source)/max(den,1e-30)
    assert balance<1e-9,balance
    return {'head':field,'relative_head':relative,'datum':node_datum,'used':used,'reactions':reactions,'water_balance':balance,'quadrature':quadrature}


def verify(mesh,case,solution,material_ids=None):
    assert solution['complete']
    independent=solve(mesh,case,solution['active_elements'],material_ids)
    actual=np.asarray(solution['field']['head_m']);used=independent['used']
    error=float(np.max(abs(actual[used]-independent['head'][used])))
    assert error<1e-8,error
    max_flux=0.;max_pore_head=0.;materials={m['id']:m for m in case['materials']}
    ids=material_ids or mesh['material_ids'];xyz=np.asarray(mesh['nodes'])
    for p in solution['points']:
        ei=p['element_index'];j=p['local_index'];q=independent['quadrature'][ei]
        nodes=q['nodes'];h=independent['relative_head'][nodes];d=h-h[0]
        gradient=np.asarray([q['b'][j,a,a::3]@d for a in range(3)])
        k=np.asarray(materials[ids[ei]]['hydraulic_conductivity_m_s']);k=(k+k.T)/2
        flux=-k@gradient;pressure=(independent['datum'][nodes[0]]-q['position'][j,2])+h[0]+q['shape'][j]@d
        max_flux=max(max_flux,float(np.max(abs(flux-p['darcy_flux_m_s']))))
        max_pore_head=max(max_pore_head,abs(pressure-p['pressure_head_m']))
        assert np.allclose(q['position'][j],p['position_m'],rtol=1e-12,atol=1e-12)
    assert max_flux<1e-12 and max_pore_head<1e-8,(max_flux,max_pore_head)
    for n,value in enumerate(solution['dirichlet_nodal_outflow_m3_s']):
        if value is not None:assert np.isclose(value,independent['reactions'][n],rtol=1e-7,atol=1e-14)
    return {'maximum_head_error_m':error,'maximum_flux_error_m_s':max_flux,'maximum_pressure_head_error_m':max_pore_head,'independent_water_balance':independent['water_balance']}


def stage_pressures(mesh,ground,run,q,offset):
    """Recompute every endpoint pore field, including frozen removed points."""
    xyz=np.asarray(mesh['nodes']);gamma=ground['initialization']['water_unit_weight_kn_m3']
    z=np.concatenate([v['position'][:,2] for v in q]);expected=None;reports=[]
    head=ground['initialization']['water_level_m']
    if head is None:head=float(xyz[:,2].min())
    audits=iter(run.get('flow_audits',[]));frame_pressures=[]
    for index,frame in enumerate(run['frames']):
        state=frame['state'];active=np.asarray(state['active_elements'])
        if index:
            stage=run['case']['stages'][index-1]
            if stage['kind']=='steady_seepage':
                audit=next(audits);assert audit['stage_index']==index-1 and audit['stage']==stage['name']
                assert audit['solution']['active_elements']==state['active_elements']
                ids=state.get('element_material_ids') or mesh['material_ids']
                assert audit['solution']['material_ids']==ids
                comparison=verify(mesh,stage['flow'],audit['solution'],ids);reports.append({'stage':stage['name'],**comparison})
                field=solve(mesh,stage['flow'],active,ids)
                if expected is None:expected=np.maximum(head-z,0)*gamma
                for ei in np.flatnonzero(active):
                    nodes=q[ei]['nodes'];h=field['relative_head'][nodes];expected[offset[ei]:offset[ei+1]]=((field['datum'][nodes[0]]-q[ei]['position'][:,2])+h[0]+q[ei]['shape']@(h-h[0]))*gamma
                head=stage['exterior_water_level_m']
            elif stage['kind']=='hydrostatic_water':
                head=stage['target_level_m']
                if expected is not None:
                    for ei in np.flatnonzero(active):expected[offset[ei]:offset[ei+1]]=np.maximum(head-q[ei]['position'][:,2],0)*gamma
        assert np.isclose(state['hydrostatic_water_level_m'],head,rtol=0,atol=1e-12)
        if expected is None:assert 'pore_pressure_kpa' not in state
        else:assert np.allclose(state['pore_pressure_kpa'],expected,rtol=1e-10,atol=1e-8)
        frame_pressures.append((np.maximum(head-z,0)*gamma if expected is None else expected).copy())
    assert next(audits,None) is None,'unexpected hydraulic audit'
    return frame_pressures,reports
