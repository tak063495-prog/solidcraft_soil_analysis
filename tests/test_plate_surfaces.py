"""CAD face-to-plate correspondence without geometry edits or lost midnodes."""
import copy
import json
from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
import fem_reference as ref
from plate_surfaces import attach


class PlateSurfaces(unittest.TestCase):
    def fixture(self):
        root=Path(__file__).resolve().parents[1]
        mesh=json.loads((root/'examples/hardening-column-tet10-mesh.json').read_text())
        faces=[list(f) for f,owners in ref.topology(mesh).items() if len(owners)==1 and all(mesh['nodes'][n][2]==1. for n in f)]
        mesh['face_sets']={'wall':faces}
        case=dict(beam_system=None,stages=[dict(kind='install_structures',members=['wall'])])
        section=dict(id='concrete',young_kpa=3e7,poisson=.2,thickness_m=.2,unit_weight_kn_m3=24.,drilling_factor=1e-4)
        spec=dict(schema='soil-fem-plate-surfaces/1',sections=[section],surfaces=[dict(id='wall',face_set='wall',owner_body='UpperClay',section_id='concrete',initially_active=False,boundary_traction_kpa=[0.,0.,-5.])])
        return mesh,case,spec

    def test_all_tri6_nodes_are_shared_and_source_geometry_and_inputs_are_unchanged(self):
        mesh,case,spec=self.fixture();before=copy.deepcopy((mesh,case,spec))
        output,stages,binding=attach(mesh,case,spec)
        self.assertEqual((mesh,case,spec),before)
        self.assertEqual({k:v for k,v in output.items() if k!='node_sets'},{k:v for k,v in mesh.items() if k!='node_sets'})
        cfg=stages['beam_system'];record=binding['surfaces'][0]
        self.assertEqual(record['plate_triangles'],4*len(mesh['face_sets']['wall']))
        self.assertAlmostEqual(record['area_m2'],1.,places=12)
        expected={int(n) for f in mesh['face_sets']['wall'] for n in ref.face_nodes(mesh,ref.topology(mesh)[tuple(sorted(f))][0],f)}
        self.assertEqual(set(map(int,binding['soil_node_to_structure_id'])),expected)
        self.assertEqual(set(stages['stages'][0]['members']),{p['id'] for p in cfg['plates']})
        for node in cfg['nodes']:
            self.assertEqual(len(output['node_sets'][node['soil_node_set']]),1)
            self.assertFalse(node.get('match_soil_node',False))
        resultant=np.array(spec['surfaces'][0]['boundary_traction_kpa'])*record['area_m2']
        np.testing.assert_array_equal(resultant,[0.,0.,-5.])

    def test_curved_shared_duplicate_and_colliding_faces_fail_without_partial_mutation(self):
        for failure in ('curved','internal','duplicate','collision'):
            mesh,case,spec=self.fixture();face=mesh['face_sets']['wall'][0]
            owner=ref.topology(mesh)[tuple(sorted(face))][0]
            if failure=='curved':mesh['nodes'][int(ref.face_nodes(mesh,owner,face)[3])][2]+=.01
            if failure=='internal':mesh['face_sets']['wall']=[list(next(f for f,owners in ref.topology(mesh).items() if len(owners)==2))]
            if failure=='duplicate':mesh['face_sets']['wall'].append(face)
            if failure=='collision':
                node=int(ref.face_nodes(mesh,owner,face)[0]);mesh['node_sets'][f'plate:mesh-node:{node}']=[node]
            before=copy.deepcopy((mesh,case,spec))
            with self.assertRaises(ValueError):attach(mesh,case,spec)
            self.assertEqual((mesh,case,spec),before)

    def test_existing_unique_structural_node_is_reused_and_ambiguous_rotation_is_rejected(self):
        mesh,case,spec=self.fixture();face=mesh['face_sets']['wall'][0]
        node=int(ref.face_nodes(mesh,ref.topology(mesh)[tuple(sorted(face))][0],face)[0]);mesh['node_sets']['BeamRoot']=[node]
        structure=dict(schema='soil-fem-beams/1',nodes=[dict(id='beam_root',position_m=mesh['nodes'][node],soil_node_set='BeamRoot')],sections=[],members=[],rotation_length_scale_m=1.,maximum_condensed_dofs=3000)
        case['beam_system']=structure
        _,result,binding=attach(mesh,case,spec)
        self.assertEqual(binding['soil_node_to_structure_id'][str(node)],'beam_root')
        self.assertEqual(sum(n['id']=='beam_root' for n in result['beam_system']['nodes']),1)
        self.assertTrue(any('beam_root' in p['nodes'] for p in result['beam_system']['plates']))
        structure['nodes'].append(dict(id='other_rotation',position_m=mesh['nodes'][node],soil_node_set='BeamRoot'))
        with self.assertRaisesRegex(ValueError,'multiple existing rotations'):attach(mesh,case,spec)

    def interface_fixture(self):
        mesh,case,spec=self.fixture();regions={r['id']:r['body'] for r in mesh['regions']}
        mesh['face_sets']['wall']=[list(f) for f,owners in ref.topology(mesh).items()
            if len(owners)==2 and {regions[mesh['region_ids'][e]] for e in owners}=={'LowerClay','UpperClay'}]
        spec['schema']='soil-fem-plate-surfaces/2'
        spec['surfaces'][0].update(attachment='shared_interface',owner_body='LowerClay',opposite_body='UpperClay')
        return mesh,case,spec

    def test_shared_interface_keeps_both_soil_traces_and_support_after_opposite_body_removal(self):
        mesh,case,spec=self.interface_fixture();before=copy.deepcopy((mesh,case,spec))
        output,result,binding=attach(mesh,case,spec)
        self.assertEqual((mesh,case,spec),before)
        self.assertEqual(binding['schema'],'soil-fem-plate-surface-binding/2')
        record=binding['surfaces'][0]
        self.assertAlmostEqual(record['area_m2'],1.,places=12)
        self.assertEqual(record['attachment'],'shared_interface')
        self.assertEqual(len(record['owner_pairs']),len(mesh['face_sets']['wall']))
        retained={n for pair in record['owner_pairs'] for n in mesh['tetrahedra10'][pair['owner_element']]}
        self.assertTrue(set(map(int,binding['soil_node_to_structure_id']))<=retained)
        xyz=np.asarray(output['nodes'])
        by_id={n['id']:n for n in result['beam_system']['nodes']}
        for plate in result['beam_system']['plates']:
            points=np.array([by_id[n]['position_m'] for n in plate['nodes']])
            self.assertGreater(np.cross(points[1]-points[0],points[2]-points[0])[2],0.)

    def test_shared_interface_reverses_orientation_with_explicit_owner(self):
        mesh,case,spec=self.interface_fixture()
        spec['surfaces'][0].update(owner_body='UpperClay',opposite_body='LowerClay')
        _,result,binding=attach(mesh,case,spec)
        by_id={n['id']:n for n in result['beam_system']['nodes']}
        for plate in result['beam_system']['plates']:
            points=np.array([by_id[n]['position_m'] for n in plate['nodes']])
            self.assertLess(np.cross(points[1]-points[0],points[2]-points[0])[2],0.)

    def test_shared_interface_requires_version_and_exact_two_body_declaration(self):
        for failure in ('legacy','missing_opposite','same_body','wrong_body','exterior_mode'):
            with self.subTest(failure=failure):
                mesh,case,spec=self.interface_fixture();surface=spec['surfaces'][0]
                if failure=='legacy':spec['schema']='soil-fem-plate-surfaces/1'
                if failure=='missing_opposite':surface.pop('opposite_body')
                if failure=='same_body':surface['opposite_body']=surface['owner_body']
                if failure=='wrong_body':surface['opposite_body']='Missing'
                if failure=='exterior_mode':surface['attachment']='exterior'
                before=copy.deepcopy((mesh,case,spec))
                with self.assertRaises(ValueError):attach(mesh,case,spec)
                self.assertEqual((mesh,case,spec),before)

    def test_shared_interface_rejects_unshared_quadratic_midnode_without_input_mutation(self):
        mesh,case,spec=self.interface_fixture();face=mesh['face_sets']['wall'][0]
        owners=ref.topology(mesh)[tuple(sorted(face))]
        owner=owners[0];mid=int(ref.face_nodes(mesh,owner,face)[3])
        new=len(mesh['nodes']);mesh['nodes'].append(copy.deepcopy(mesh['nodes'][mid]))
        cell=mesh['tetrahedra10'][owner];cell[cell.index(mid)]=new
        before=copy.deepcopy((mesh,case,spec))
        with self.assertRaisesRegex(ValueError,'nonconforming'):attach(mesh,case,spec)
        self.assertEqual((mesh,case,spec),before)


if __name__=='__main__':unittest.main()
