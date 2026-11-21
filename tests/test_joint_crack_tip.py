"""Declared Tet10 edge tips: topology and independent paired-face operators."""
import copy
from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
import fem_reference as ref
from joint_reference import System
from joint_tip_fixture import conformal, specification, fixture
from split_joint_mesh import split


class JointCrackTip(unittest.TestCase):
    def test_multiple_faces_internal_tip_glued_region_geometry_and_material_are_exact(self):
        for same_body in (False, True):
            source, spec, mesh, configuration = fixture(same_body=same_body)
            declaration = configuration['joints'][0]['crack_tip']
            self.assertEqual(len(declaration['faces']), 2)
            self.assertEqual(len(declaration['edges']), 1)
            tips = set(declaration['edges'][0])
            original = np.asarray(source['nodes'])[ref.cells(source)]
            np.testing.assert_array_equal(np.asarray(mesh['nodes'])[ref.cells(mesh)], original)
            for key in ('material_ids', 'region_ids', 'regions'):
                self.assertEqual(mesh[key], source[key])
            selected = {tuple(sorted(f['nodes'])) for f in spec['interfaces'][0]['faces']}
            for face, owners in ref.topology(source).items():
                if len(owners) != 2:
                    continue
                a, b = owners
                na = ref.face_nodes(source, a, face)
                mapped = []
                for owner in (a, b):
                    mapped.append({int(n):int(ref.cells(mesh)[owner, list(ref.cells(source)[owner]).index(n)]) for n in na})
                if face in selected:
                    self.assertEqual({n for n in na if mapped[0][int(n)]==mapped[1][int(n)]}, set(na)&tips)
                else:
                    self.assertEqual(mapped[0], mapped[1])
            self.assertEqual(sum(ref.volume(source,e)['weight'].sum() for e in range(len(source['tetrahedra']))),
                             sum(ref.volume(mesh,e)['weight'].sum() for e in range(len(mesh['tetrahedra']))))

    def test_missing_partial_wrong_midpoint_interior_edge_and_duplicate_are_rejected(self):
        source, spec, _, _ = fixture()
        bad = copy.deepcopy(spec); bad['interfaces'][0].pop('tip_edges')
        with self.assertRaisesRegex(ValueError, 'undeclared shared'):
            split(source,bad)
        for mutation in ('mid', 'partial', 'duplicate', 'owner', 'interior', 'face'):
            bad = copy.deepcopy(spec); d=bad['interfaces'][0]
            if mutation == 'mid': d['tip_edges'][0][2]=d['tip_edges'][0][0]
            elif mutation == 'partial': d['tip_edges'][0]=d['tip_edges'][0][:2]
            elif mutation == 'duplicate': d['tip_edges'].append(d['tip_edges'][0].copy())
            elif mutation == 'owner': d['faces'][0]['owner_b']=d['faces'][0]['owner_a']
            elif mutation == 'face': d['faces'].append(copy.deepcopy(d['faces'][0]))
            else:
                faces=d['faces'];edge=sorted(set(faces[0]['nodes'])&set(faces[1]['nodes']))
                cell=ref.cells(source)[faces[0]['owner_a']]
                k=next(k for k,(i,j) in enumerate(ref.EDGES) if {int(cell[i]),int(cell[j])}==set(edge))
                d['tip_edges']=[edge+[int(cell[k+4])]]
            with self.subTest(mutation=mutation), self.assertRaises(ValueError): split(source,bad)

    def test_tri4_tip_nonconforming_midpoint_and_fully_tied_faces_are_rejected(self):
        source,spec,_,_=fixture()
        linear=copy.deepcopy(source);linear.pop('tetrahedra10');linear.pop('quadrature_order');linear['schema']='soil-fem-mesh/2'
        with self.assertRaises(ValueError):split(linear,spec)
        bad=copy.deepcopy(source);bad['tetrahedra10'][0][4]=bad['tetrahedra10'][0][5]
        with self.assertRaises(ValueError):split(bad,spec)
        # Declaring all boundary edges would erase jump DOFs on the first face.
        one=copy.deepcopy(spec);d=one['interfaces'][0];d['faces']=d['faces'][:1];ns=ref.face_nodes(source,d['faces'][0]['owner_a'],d['faces'][0]['nodes'])
        d['tip_edges']=[[int(ns[i]),int(ns[j]),int(ns[3+k])] for k,(i,j) in enumerate(((0,1),(1,2),(2,0)))]
        with self.assertRaisesRegex(ValueError, 'fully tied'):split(source,one)

    def test_independent_operator_exact_cancellation_virtual_work_and_rigid_motion(self):
        _,_,mesh,cfg=fixture()
        system=System(mesh,cfg)
        tips=cfg['joints'][0]['crack_tip']['edges'][0]
        np.testing.assert_array_equal(system.c[:,:,np.ravel([[3*n+k for k in range(3)] for n in tips])],0.)
        rng=np.random.default_rng(1701);u=rng.normal(size=(len(mesh['nodes']),3));du=rng.normal(size=u.shape);traction=rng.normal(size=(len(system.points),3))
        jump=np.array([np.asarray(p['shape'])@(u[p['nodes_b']]-u[p['nodes_a']]) for p in system.points])
        virtual_jump=np.array([np.asarray(p['shape'])@(du[p['nodes_b']]-du[p['nodes_a']]) for p in system.points])
        np.testing.assert_allclose(np.einsum('qag,g->qa',system.c,u.ravel()),jump,atol=1e-14,rtol=0)
        force=np.einsum('qag,qa,q->g',system.c,traction,system.weight).reshape(-1,3)
        self.assertLess(abs(float(np.sum(du*force)-np.einsum('qa,qa,q->',virtual_jump,traction,system.weight))),1e-12)
        np.testing.assert_allclose(force.sum(axis=0),0.,atol=1e-12,rtol=0)
        np.testing.assert_allclose(np.cross(mesh['nodes'],force).sum(axis=0),0.,atol=1e-12,rtol=0)
        xyz=np.asarray(mesh['nodes']);rigid=[.1,-.2,.3]+np.cross([.4,-.5,.6],xyz)
        np.testing.assert_allclose(np.einsum('qag,g->qa',system.c,rigid.ravel()),0.,atol=1e-12,rtol=0)

    def test_reference_refuses_undeclared_or_tampered_tip_topology(self):
        _,_,mesh,cfg=fixture()
        bad=copy.deepcopy(cfg);bad['joints'][0].pop('crack_tip')
        with self.assertRaises(AssertionError):System(mesh,bad)
        for key in ('edges','owner_a','face_b'):
            bad=copy.deepcopy(cfg);d=bad['joints'][0]['crack_tip']
            if key=='edges':d['edges'][0][2]=d['edges'][0][0]
            elif key=='owner_a':d['faces'][0][key]=d['faces'][0]['owner_b']
            else:d['faces'][0][key]=d['faces'][0]['face_a']
            with self.subTest(key=key), self.assertRaises(AssertionError):System(mesh,bad)

    def test_multi_edge_front_and_nondeclared_boundary_are_distinct(self):
        _,_,mesh,cfg=fixture(y=(0.,.5,1.))
        self.assertEqual(len(cfg['joints'][0]['crack_tip']['edges']),2)
        System(mesh,cfg)
        source=conformal(y=(0.,.5,1.));spec=specification(source)
        spec['interfaces'][0]['tip_edges'].pop()
        with self.assertRaisesRegex(ValueError,'undeclared shared'):split(source,spec)

    def test_branched_front_disconnected_edge_star_and_face_interior_are_rejected(self):
        source = conformal(y=(0.,1.,2.))
        spec = specification(source, full=True)
        faces = spec['interfaces'][0]['faces']
        pair = next((a,b) for a in faces for b in faces if len(set(a['nodes']) & set(b['nodes']))==1)
        d = spec['interfaces'][0]; d['faces'] = list(pair); d['tip_edges'] = []
        for f in pair:
            ns = ref.face_nodes(source,f['owner_a'],f['nodes'])
            d['tip_edges'].extend([[int(ns[i]),int(ns[j]),int(ns[3+k])] for k,(i,j) in enumerate(((0,1),(1,2),(2,0)))])
        with self.assertRaisesRegex(ValueError,'branching'):split(source,spec)
        source,spec,_,_=fixture()
        bad=copy.deepcopy(spec);bad['interfaces'][0]['faces'][0]['nodes'][0]=bad['interfaces'][0]['tip_edges'][0][2]
        with self.assertRaisesRegex(ValueError,'face/owner'):split(source,bad)
        # A second volume fan touching only the tip edge is nonmanifold.
        v,w,mid=spec['interfaces'][0]['tip_edges'][0]
        new=len(source['nodes']);source['nodes'].extend([[1.3,.3,.2],[1.4,.6,.8]])
        tet=[v,w,new,new+1];cell=tet.copy()
        for i,j in ref.EDGES:
            if {tet[i],tet[j]}=={v,w}:cell.append(mid)
            else:
                cell.append(len(source['nodes']))
                source['nodes'].append(((np.asarray(source['nodes'][tet[i]])+source['nodes'][tet[j]])/2).tolist())
        source['tetrahedra'].append(tet);source['tetrahedra10'].append(cell)
        source['material_ids'].append(1);source['region_ids'].append(1)
        with self.assertRaisesRegex(ValueError,'nonmanifold tip edge star'):split(source,spec)

    def test_branching_cannot_be_hidden_by_multiple_interface_ids(self):
        source=conformal(same_body=True)
        nodes=[[0.,0.,0.]];cells=[];interfaces=[]
        for k,(x,y) in enumerate(((1.,0.),(0.,1.),(-1.,0.))):
            n=len(nodes);nodes.extend([[x,y,0.],[-y,x,0.],[0.,0.,-1.],[0.,0.,1.]])
            cells.extend([[0,n+1,n,n+2],[0,n,n+1,n+3]])
            interfaces.append(dict(id=f'J{k}',faces=[dict(nodes=[0,n,n+1],owner_a=2*k,owner_b=2*k+1)]))
        edges={};quadratic=[]
        for tet in cells:
            cell=tet.copy()
            for i,j in ref.EDGES:
                edge=tuple(sorted((tet[i],tet[j])))
                if edge not in edges:
                    edges[edge]=len(nodes);nodes.append(((np.asarray(nodes[edge[0]])+nodes[edge[1]])/2).tolist())
                cell.append(edges[edge])
            quadratic.append(cell)
        source.update(nodes=nodes,tetrahedra=cells,tetrahedra10=quadratic,material_ids=[1]*6,region_ids=[1]*6,node_sets={'all':list(range(len(nodes)))},face_sets={})
        for k,d in enumerate(interfaces):
            n=1+4*k;d['tip_edges']=[[0,n,edges[0,n]]]
        with self.assertRaisesRegex(ValueError,'branching crack tip across interfaces'):
            split(source,dict(schema='soil-fem-joint-split/1',interfaces=interfaces))
        # The independent constructor must reject the complete network before
        # individual side matching can obscure this three-arm declaration.
        template=fixture()[3]['joints'][0];cfg=fixture()[3];cfg['joints']=[]
        for d in interfaces:
            joint=copy.deepcopy(template);joint['id']=d['id'];joint['crack_tip']['edges']=d['tip_edges'];cfg['joints'].append(joint)
        with self.assertRaisesRegex(AssertionError,'branching crack tip across joints'):System(source,cfg)


if __name__=='__main__':unittest.main()
