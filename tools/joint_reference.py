"""Independent paired-face geometry and mechanical audits of dry JOINT states.

NumPy basis/Jacobian construction checks the exported raw surface quadrature.
This module does not call the Rust constitutive update or read its tangent.
"""
import numpy as np
import fem_reference as ref


def declared_tip(mesh, definition, owners):
    """Independent combinatorial certificate, never inferred from coincidence."""
    declaration = definition.get('crack_tip')
    if declaration is None:
        return set(), {}
    assert mesh.get('tetrahedra10') and declaration['edges']
    faces = declaration['faces']
    a = [tuple(sorted(f['face_a'])) for f in faces]
    b = [tuple(sorted(f['face_b'])) for f in faces]
    assert len(a) == len(set(a)) and len(b) == len(set(b))
    assert set(a) == {tuple(sorted(f)) for f in mesh['face_sets'][definition['side_a']]}
    assert set(b) == {tuple(sorted(f)) for f in mesh['face_sets'][definition['side_b']]}
    correspondence = {}
    surface = {}
    for f, fa, fb in zip(faces, a, b):
        oa, ob = f['owner_a'], f['owner_b']
        assert oa != ob and owners[fa] == [oa] and owners[fb] == [ob]
        correspondence[fa] = (fb, oa, ob)
        nodes = ref.face_nodes(mesh, oa, fa)
        for k, (i, j) in enumerate(((0, 1), (1, 2), (2, 0))):
            edge = tuple(sorted((int(nodes[i]), int(nodes[j]))))
            surface.setdefault(edge, []).append(int(nodes[k+3]))
    assert all(len(m) <= 2 and len(set(m)) == 1 for m in surface.values())
    allowed, unique, degree = set(), set(), {}
    cells = ref.cells(mesh)
    for v, w, mid in declaration['edges']:
        edge = tuple(sorted((v, w)))
        assert len({v, w, mid}) == 3 and edge not in unique and surface.get(edge) == [mid]
        unique.add(edge)
        allowed.update((v, w, mid))
        for n in (v, w):
            degree[n] = degree.get(n, 0) + 1
        star = {e for e, cell in enumerate(cells) if v in cell[:4] and w in cell[:4]}
        links = {e: set() for e in star}
        for face, os in owners.items():
            if v in face and w in face and len(os) == 2:
                x, y = os
                links[x].add(y)
                links[y].add(x)
        for f in faces:
            if v in f['face_a'] and w in f['face_a']:
                assert v in f['face_b'] and w in f['face_b']
                x, y = f['owner_a'], f['owner_b']
                links[x].add(y)
                links[y].add(x)
        visited, pending = set(), [min(star)]
        while pending:
            e = pending.pop()
            if e not in visited:
                visited.add(e)
                pending.extend(links[e])
        assert visited == star and all(len(s) <= 2 for s in links.values())
    assert all(d <= 2 for d in degree.values())
    return allowed, correspondence


class System:
    def __init__(self,mesh,configuration,layout=None):
        self.mesh=mesh;self.configuration=configuration;xyz=np.asarray(mesh['nodes']);owners=ref.topology(mesh)
        front_edges, front_degree = set(), {}
        for d in configuration['joints']:
            for v, w, _ in (d.get('crack_tip') or {}).get('edges', []):
                edge = tuple(sorted((v, w)))
                assert edge not in front_edges, 'duplicate crack-tip edge across joints'
                front_edges.add(edge)
                for n in edge:front_degree[n] = front_degree.get(n, 0) + 1
        assert all(degree <= 2 for degree in front_degree.values()), 'branching crack tip across joints'
        assert all(len(os) <= 2 for os in owners.values())
        if any(d.get('crack_tip') is not None for d in configuration['joints']):
            cells = ref.cells(mesh); corners = set(cells[:, :4].flat)
            edges, middles = {}, {}
            assert cells.shape[1] == 10
            for cell in cells:
                assert len(set(cell)) == 10
                for k, (i, j) in enumerate(ref.EDGES):
                    edge = tuple(sorted((int(cell[i]), int(cell[j])))); mid = int(cell[k+4])
                    assert mid not in corners and edges.get(edge, mid) == mid and middles.get(mid, edge) == edge
                    edges[edge] = mid; middles[mid] = edge
        self.points=[];self.counterparts={};used=set()
        for d in configuration['joints']:
            allowed, correspondence = declared_tip(mesh, d, owners)
            for fi,face in enumerate(mesh['face_sets'][d['side_a']]):
                fa=tuple(sorted(face));oa,=owners[fa];na=ref.face_nodes(mesh,oa,fa);pa=xyz[na];matches=[]
                for fb in mesh['face_sets'][d['side_b']]:
                    fb=tuple(sorted(fb));ob,=owners[fb];nb=ref.face_nodes(mesh,ob,fb)
                    if len(na)!=len(nb):continue
                    distance=np.linalg.norm(pa[:,None]-xyz[nb][None,:],axis=2);permutation=distance.argmin(axis=1)
                    if np.max(distance[np.arange(len(na)),permutation])<=configuration['geometric_tolerance_m'] and len(set(permutation))==len(na):
                        if d.get('crack_tip') is not None:
                            assert np.all(permutation[:3] < 3)
                            assert sum(permutation[i] > permutation[j] for i in range(3) for j in range(i+1,3)) % 2 == 1
                            for k, (i, j) in enumerate(((0,1),(1,2),(2,0))):
                                edge = next(e for e,(v,w) in enumerate(((0,1),(1,2),(2,0))) if {v,w}=={int(permutation[i]),int(permutation[j])})
                                assert permutation[k+3] == edge+3
                        matches.append((fb,ob,nb[permutation]))
                assert len(matches)==1;fb,ob,nb=matches[0];assert fa not in used and fb not in used;used.update((fa,fb))
                shared = set(int(n) for n in na[na == nb])
                assert shared == set(int(n) for n in na) & allowed and not np.all(na == nb)
                if correspondence:
                    assert correspondence[fa] == (fb, oa, ob)
                self.counterparts[fa]=ob;self.counterparts[fb]=oa
                if len(na)==3:rule=[(np.array(r),1/6) for r in ((1/6,1/6),(2/3,1/6),(1/6,2/3))]
                else:rule=[(np.array([u,(1-u)*v]),wu*wv) for u,wu in zip(*ref.jacobi(4,1)) for v,wv in zip(*ref.jacobi(4,0))]
                for r,w in rule:
                    l=np.r_[1-r.sum(),r];dl=np.array([[-1.,-1.],[1.,0.],[0.,1.]])
                    shape=l.copy();gradient=dl.copy()
                    if len(na)==6:
                        shape=np.r_[l*(2*l-1),[4*l[a]*l[b] for a,b in ((0,1),(1,2),(2,0))]]
                        gradient=np.vstack(((4*l-1)[:,None]*dl,[4*(l[a]*dl[b]+l[b]*dl[a]) for a,b in ((0,1),(1,2),(2,0))]))
                    j=(pa-pa[0]).T@gradient;area=np.cross(j[:,0],j[:,1])*w;weight=np.linalg.norm(area);normal=area/weight
                    self.points.append({'joint_id':d['id'],'face_index':fi,'position_m':(pa[0]+shape@(pa-pa[0])).tolist(),
                        'normal':normal.tolist(),'weight_m2':float(weight),'nodes_a':na.tolist(),'nodes_b':nb.tolist(),'shape':shape.tolist(),'owner_a':oa,'owner_b':ob,'material':d['material']})
        self.normal=np.asarray([p['normal'] for p in self.points]);self.weight=np.asarray([p['weight_m2'] for p in self.points])
        self.kn=np.asarray([p['material']['normal_stiffness_kpa_m'] for p in self.points]);self.kt=np.asarray([p['material']['tangential_stiffness_kpa_m'] for p in self.points]);self.mu=np.asarray([p['material']['friction_coefficient'] for p in self.points])
        self.c=np.zeros((len(self.points),3,xyz.size))
        for i,p in enumerate(self.points):
            for a,b,w in zip(p['nodes_a'],p['nodes_b'],p['shape']):
                for k in range(3):self.c[i,k,3*a+k]-=w;self.c[i,k,3*b+k]+=w
        for p, c in zip(self.points, self.c):
            for node in set(p['nodes_a']) & set(p['nodes_b']):
                assert np.array_equal(c[:, 3*node:3*node+3], np.zeros((3, 3)))
        if layout is not None:
            assert len(layout)==len(self.points)
            for expected,actual in zip(self.points,layout):
                for key in ('joint_id','face_index','nodes_a','nodes_b','owner_a','owner_b'):assert expected[key]==actual[key],key
                for key in ('normal','shape','position_m','weight_m2'):assert np.allclose(expected[key],actual[key],atol=1e-12,rtol=1e-10),key

    def rigid_support(self,case,state,expected=None):
        """Independent dense SVD of constrained component rigid motions.

        This certificate covers JOINT and displacement constraints. Beam and
        anchor stiffness require their own independent assembly before use.
        """
        assert not case.get('beam_system') and not case.get('anchor_system')
        world=np.asarray(self.mesh['nodes']);xyz=world-world[0];cells=ref.cells(self.mesh);active=np.asarray(state['active_elements'])
        owners=ref.topology(self.mesh);neighbors={e:set() for e in np.flatnonzero(active)}
        for pair in owners.values():
            live=[e for e in pair if active[e]]
            if len(live)==2:
                a,b=live;neighbors[a].add(b);neighbors[b].add(a)
        components=[];pending=set(neighbors);component={}
        while pending:
            stack=[min(pending)];elements=set()
            while stack:
                e=stack.pop()
                if e not in pending:continue
                pending.remove(e);elements.add(e);stack.extend(neighbors[e])
            for e in elements:component[e]=len(components)
            components.append(np.unique(cells[sorted(elements)]))
        assert 0<len(components)<=500
        centers=[xyz[nodes].mean(axis=0) for nodes in components]
        radii=[np.linalg.norm(xyz[nodes]-center,axis=1).max() for nodes,center in zip(components,centers)]
        dimension=6*len(components)
        def motion(c,point):
            r=(point-centers[c])/radii[c];matrix=np.zeros((3,dimension))
            matrix[:,6*c:6*c+6]=np.c_[np.eye(3),np.cross(np.eye(3),r).T]
            return matrix
        rows=[];primary={}
        for c,nodes in enumerate(components):
            for n in nodes:
                if n in primary:rows.extend(motion(c,xyz[n])-motion(primary[n],xyz[n]))
                else:primary[n]=c
        local_case={**case,'boundary':case.get('boundary') or {'constraints':[]}}
        basis,_=ref.staged_constraint_space(self.mesh,local_case,state)
        for n,b in basis.items():rows.extend((np.eye(3)-b@b.T)@motion(primary[n],xyz[n]))
        stresses=np.asarray([s['effective_stress_compression_positive_kpa'] for s in state['material_states']])
        # Include only the active volume quadrature history, as in the force norm.
        live_points=np.repeat(active,len(ref.volume_rule(self.mesh)[1]))
        stress_scale=max(1.,float(np.sqrt((stresses[live_points]**2*np.array([1,1,1,2,2,2])).sum(axis=1)).max()))
        pressure_tolerance=max(1e-10,case['settings']['relative_force_tolerance']*stress_scale)
        for p,s in zip(self.points,state['joints']['points']):
            if not (active[p['owner_a']] and active[p['owner_b']]):continue
            n=np.asarray(p['normal']);q=np.asarray(s['traction_kpa']);pressure=-q@n
            if pressure<=pressure_tolerance:continue
            ra=sum(w*motion(component[p['owner_a']],xyz[n]) for n,w in zip(p['nodes_a'],p['shape']))
            rb=sum(w*motion(component[p['owner_b']],xyz[n]) for n,w in zip(p['nodes_b'],p['shape']))
            rows.append(n@(rb-ra));cap=p['material']['friction_coefficient']*pressure;tau=q+n*pressure
            shear_tolerance=pressure_tolerance*(1+p['material']['friction_coefficient'])
            if cap<=shear_tolerance:continue
            if np.linalg.norm(tau)<cap-max(cap*1e-8,shear_tolerance):
                t=np.cross(n,np.eye(3)[0 if abs(n[0])<.8 else 1]);t/=np.linalg.norm(t)
                rows.extend([t@(rb-ra),np.cross(n,t)@(rb-ra)])
            elif np.linalg.norm(tau)>0:
                t=np.cross(n,tau);t/=np.linalg.norm(t);rows.append(t@(rb-ra))
        matrix=np.asarray(rows).reshape(-1,dimension);length=np.linalg.norm(matrix,axis=1)
        matrix=matrix[length>1e-10]/length[length>1e-10,None]
        singular=np.linalg.svd(matrix,compute_uv=False)
        rank=int((singular>1e-9*max(float(singular.max(initial=0)),1)).sum())
        report={'continuum_components':len(components),'rigid_modes':dimension,'restrained_modes':rank,'unrestrained_modes':dimension-rank}
        if expected is not None:assert report==expected,(report,expected)
        return report

    def audit(self,state,report=None):
        states=state['joints']['points'];assert len(states)==len(self.points);u=np.asarray(state['displacement_m']);active=np.asarray(state['active_elements'])
        computed=np.einsum('qag,g->qa',self.c,u.reshape(-1));force=np.zeros(u.size);totals=np.zeros(4)
        for i,(p,s) in enumerate(zip(self.points,states)):
            n=self.normal[i];jump=np.asarray(s['jump_m']);q=np.asarray(s['traction_kpa']);gn=jump@n;qn=q@n;tau=q-qn*n
            slip=np.asarray(s['plastic_slip_m']);live=bool(active[p['owner_a']] and active[p['owner_b']])
            assert np.allclose(s['normal'],n,atol=1e-12)
            if live:assert np.allclose(jump,computed[i],atol=1e-10,rtol=0)
            assert np.isclose(qn,self.kn[i]*min(gn,0),atol=1e-7,rtol=1e-9)
            assert np.linalg.norm(tau)<=self.mu[i]*max(-qn,0)+1e-7
            assert abs(slip@n)<1e-10 and s['accumulated_plastic_slip_m']+1e-10>=np.linalg.norm(slip)
            assert s['plastic_dissipation_kpa_m']>=-1e-12 and s['numerical_dissipation_kpa_m']>=-1e-12
            if gn>1e-12:assert np.linalg.norm(q)<1e-7
            energy=qn**2/(2*self.kn[i])+tau@tau/(2*self.kt[i])
            assert np.isclose(s['endpoint_work_kpa_m'],energy+s['plastic_dissipation_kpa_m']+s['numerical_dissipation_kpa_m'],atol=1e-9,rtol=1e-8)
            if live:
                force+=self.c[i].T@q*self.weight[i]
                totals+=self.weight[i]*np.array([energy,s['plastic_dissipation_kpa_m'],s['numerical_dissipation_kpa_m'],s['endpoint_work_kpa_m']])
            if report is not None:
                r=report['points'][i];assert r['index']==i and r['active']==live
                for key,value in {'gap_m':gn,'opening_m':max(gn,0),'normal_pressure_kpa':-qn,'shear_traction_kpa':tau,'elastic_energy_kpa_m':energy,
                    'contact_reference_offset_m':jump-gn*n-tau/self.kt[i]-slip}.items():assert np.allclose(r[key],value,atol=1e-9,rtol=1e-8),key
        force=force.reshape(-1,3);moment=np.cross(np.asarray(self.mesh['nodes']),force).sum(axis=0)
        assert np.linalg.norm(force.sum(axis=0))<1e-7 and np.linalg.norm(moment)<1e-7
        if report is not None:
            for key,value in zip(('elastic_energy_kn_m','plastic_dissipation_kn_m','numerical_dissipation_kn_m','endpoint_work_kn_m'),totals):assert np.isclose(report[key],value,atol=1e-9,rtol=1e-8)
            assert np.allclose(report['resultant_internal_force_kn'],force.sum(axis=0),atol=1e-8)
            assert np.allclose(report['resultant_internal_moment_kn_m'],moment,atol=1e-8)
        return force


def coupled_forces(system,state,force,load,report=None):
    return force+system.audit(state,report),load
