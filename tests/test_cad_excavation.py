"""Exercise the real OpenCASCADE STEP/Boolean/Tet4 path, without CAD mocks."""
import contextlib
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
import subprocess
import gmsh
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from mesh_layers import generate


def source(folder,name,kind,values,index):
    gmsh.initialize()
    try:
        gmsh.option.setNumber('General.Terminal',0)
        tag=(gmsh.model.occ.addCylinder if kind=='cylinder' else gmsh.model.occ.addBox)(*values)
        gmsh.model.occ.synchronize();path=folder/f'{name}.step';gmsh.write(str(path))
        bounds=gmsh.model.occ.getBoundingBox(3,tag)
        return {'id':index,'body':name,'step':path.name,'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                'bounds_m':(np.array(bounds).reshape(2,3)*.001).tolist(),'volume_m3':gmsh.model.occ.getMass(3,tag)*1e-9}
    finally:gmsh.finalize()


class CadExcavation(unittest.TestCase):
    def test_body_plane_sets_use_all_quadratic_nodes_and_preserve_existing_data_atomically(self):
        import copy
        from body_boundary_sets import add_body_boundary_sets
        points=[[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],[0.,0.,-1.]]
        corners=[[0,1,2,3],[0,2,1,4]];middle={};quadratic=[]
        for cell in corners:
            nodes=cell.copy()
            for a,b in ((0,1),(1,2),(2,0),(0,3),(1,3),(2,3)):
                edge=tuple(sorted((cell[a],cell[b])))
                if edge not in middle:
                    middle[edge]=len(points);points.append(((np.asarray(points[edge[0]])+points[edge[1]])/2).tolist())
                nodes.append(middle[edge])
            quadratic.append(nodes)
        mesh=dict(nodes=points,tetrahedra=corners,tetrahedra10=quadratic,region_ids=[1,2],material_ids=[1,2],
            regions=[dict(id=1,body='Upper',bounds_m=[[0.,0.,0.],[1.,1.,1.]]),
                dict(id=2,body='Lower',bounds_m=[[0.,0.,-1.],[1.,1.,0.]])],
            node_sets={'existing':[0,1]},face_sets={'existing':[[0,1,2]]})
        original=copy.deepcopy(mesh);add_body_boundary_sets(mesh)
        self.assertEqual(mesh['face_sets']['body:Upper:zmin'],[[0,1,2]])
        self.assertEqual(mesh['face_sets']['body:Lower:zmax'],[[0,1,2]])
        expected=sorted([0,1,2,middle[(0,1)],middle[(1,2)],middle[(0,2)]])
        self.assertEqual(mesh['node_sets']['body:Upper:zmin'],expected)
        self.assertNotIn([0,1,2],mesh['face_sets']['body:Upper:exterior'])
        self.assertNotIn('body:Upper:zmax',mesh['face_sets'])
        pruned=copy.deepcopy(mesh)
        for key in ('node_sets','face_sets'):pruned[key]={k:v for k,v in pruned[key].items() if not k.startswith('body:')}
        self.assertEqual(pruned,original)
        once=copy.deepcopy(mesh);add_body_boundary_sets(mesh);self.assertEqual(mesh,once)
        collision=copy.deepcopy(original);collision['node_sets']['body:Upper:zmin']=[0]
        before=copy.deepcopy(collision)
        with self.assertRaisesRegex(ValueError,'collides'):add_body_boundary_sets(collision)
        self.assertEqual(collision,before)
        curved=copy.deepcopy(original);curved['nodes'][middle[(0,1)]][2]=.1
        add_body_boundary_sets(curved)
        self.assertNotIn('body:Upper:zmin',curved['face_sets'])
        self.assertNotIn('body:Lower:zmax',curved['face_sets'])

    def test_an_empty_owned_fragment_is_rejected_before_any_mesh_is_exported(self):
        from unittest.mock import patch
        with tempfile.TemporaryDirectory() as d:
            path,_=self.prepare(Path(d));output=path.parent/'lost-fragment.json'
            real=gmsh.model.mesh.getElements;calls=[]
            def missing(dim=-1,tag=-1):
                if dim==3:
                    calls.append(tag)
                    if len(calls)==2:return np.array([],dtype=int),[],[]
                return real(dim,tag)
            with patch.object(gmsh.model.mesh,'getElements',side_effect=missing),contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaisesRegex(ValueError,'has no tetrahedral elements'):
                    generate(path,output,2.)
            self.assertEqual(len(calls),2)
            self.assertFalse(output.exists())

    def prepare(self,folder,kind='box',values=(3000,3000,6000,4000,4000,3000)):
        regions=[]
        for i,z in enumerate((0,4000)):
            r=source(folder,f'Clay{i}', 'box',(0,0,z,10000,10000,4000),i+1);r['material_id']=i+1;regions.append(r)
        tool=source(folder,'Pit',kind,values,3)
        manifest={'schema':'soil-fem-cad-bodies/1','units':'m-kN-kPa','ground_sha256':'a'*64,'regions':regions,'excavation_tools':[tool]}
        path=folder/'bodies.json';path.write_text(json.dumps(manifest));return path,manifest

    def mesh(self,path):
        output=path.parent/'mesh.json'
        with contextlib.redirect_stdout(io.StringIO()):generate(path,output,1.5)
        mesh=json.loads(output.read_text());p=np.array(mesh['nodes'])[mesh['tetrahedra']]
        volumes=np.linalg.det((p[:,1:]-p[:,:1]).transpose(0,2,1))/6
        self.assertTrue(np.all(volumes>0));self.assertAlmostEqual(volumes.sum(),800.,places=8)
        self.assertEqual(set(mesh['material_ids']),{1,2})
        self.assertTrue(mesh['node_sets']['Pit:floor']);self.assertTrue(mesh['node_sets']['Pit:surrounding_surface'])
        return mesh,p.mean(axis=1),volumes

    def test_clipped_tool_does_not_add_soil_or_change_material_ownership(self):
        with tempfile.TemporaryDirectory() as d:
            path,_=self.prepare(Path(d));mesh,c,v=self.mesh(path);ids=mesh['excavation_sets']['Pit']
            self.assertAlmostEqual(v[ids].sum(),32.,places=9)
            self.assertTrue(np.all(c[ids,2]>=6));self.assertEqual(set(np.array(mesh['material_ids'])[ids]),{2})
            provenance=json.loads(path.with_name('mesh.provenance.json').read_text())
            self.assertAlmostEqual(provenance['excavation_boolean_volume_m3']['Pit'],32.,places=9)

    def test_curved_tool_uses_boolean_membership_rather_than_its_bounding_box(self):
        with tempfile.TemporaryDirectory() as d:
            path,_=self.prepare(Path(d),'cylinder',(5000,5000,6000,0,0,3000,2000))
            mesh,c,v=self.mesh(path);ids=mesh['excavation_sets']['Pit'];r2=np.square(c[:,:2]-5).sum(axis=1)
            self.assertTrue(np.all(r2[ids]<4.));self.assertLess(v[ids].sum(),8*np.pi);self.assertGreater(v[ids].sum(),.90*8*np.pi)
            outside=(c[:,2]>6)&np.all(np.abs(c[:,:2]-5)<2,axis=1)&(r2>4)
            self.assertTrue(outside.any());self.assertFalse(set(np.flatnonzero(outside))&set(ids))

    def test_tool_crosses_material_interface_without_changing_materials(self):
        with tempfile.TemporaryDirectory() as d:
            path,_=self.prepare(Path(d),values=(3000,3000,3000,4000,4000,6000))
            mesh,c,v=self.mesh(path);ids=mesh['excavation_sets']['Pit']
            self.assertAlmostEqual(v[ids].sum(),80.,places=9);self.assertEqual(set(np.array(mesh['material_ids'])[ids]),{1,2})
            interface=np.asarray(mesh['face_sets']['interface:Clay0:Clay1']);points=np.asarray(mesh['nodes'])[interface]
            self.assertTrue(np.allclose(points[:,:,2],4.,rtol=0,atol=1e-10))
            self.assertAlmostEqual(np.linalg.norm(np.cross(points[:,1]-points[:,0],points[:,2]-points[:,0]),axis=1).sum()/2,100.,places=9)

    def test_no_intersection_and_tampered_source_are_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            path,manifest=self.prepare(Path(d),values=(20000,20000,6000,4000,4000,3000))
            with self.assertRaisesRegex(ValueError,'no positive-volume'):generate(path,path.parent/'empty.json',2.)
            manifest['excavation_tools'][0]['source_sha256']='0'*64;path.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError,'hash mismatch'):generate(path,path.parent/'tampered.json',2.)

    def test_curved_tet10_projects_edge_nodes_and_passes_rust_jacobian_certification(self):
        import fem_reference as ref
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);path,_=self.prepare(folder,'cylinder',(5000,5000,6000,0,0,3000,2000))
            output=folder/'quadratic.json'
            with contextlib.redirect_stdout(io.StringIO()):generate(path,output,2.,order=2)
            mesh=json.loads(output.read_text());xyz=np.array(mesh['nodes']);cells=np.array(mesh['tetrahedra10'])
            self.assertEqual(mesh['schema'],'soil-fem-mesh/3')
            self.assertTrue(any(np.linalg.norm(xyz[t[4]]-(xyz[t[0]]+xyz[t[1]])/2)>1e-5 for t in cells))
            total=sum(ref.volume(mesh,i,independent=True)['weight'].sum() for i in range(len(cells)))
            self.assertAlmostEqual(total,800.,places=8)
            solver=Path(__file__).resolve().parents[1]/'target/release/soil-fem-3d.exe'
            result=subprocess.run([str(solver),'validate',str(output)],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)


if __name__=='__main__':unittest.main()
