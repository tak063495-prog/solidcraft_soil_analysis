"""Full-interface CAD mesh splitting keeps geometry and undeclared continuity."""
import copy
import json
from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
import fem_reference as ref
from split_joint_mesh import split


class JointSplit(unittest.TestCase):
    def test_linear_and_quadratic_interfaces_preserve_material_geometry_and_volume(self):
        root=Path(__file__).resolve().parents[1]
        for quadratic in (False,True):
            mesh=json.loads((root/'examples/mc-column-mesh.json').read_text())
            if quadratic:
                edges={};cells=[]
                for t in mesh['tetrahedra']:
                    cell=t.copy()
                    for a,b in ref.EDGES:
                        key=tuple(sorted((t[a],t[b])))
                        if key not in edges:
                            edges[key]=len(mesh['nodes']);mesh['nodes'].append(((np.asarray(mesh['nodes'][key[0]])+mesh['nodes'][key[1]])/2).tolist())
                        cell.append(edges[key])
                    cells.append(cell)
                mesh.update(schema='soil-fem-mesh/3',tetrahedra10=cells,quadrature_order=3)
            spec=json.loads((root/'examples/joint-column-split.json').read_text())
            before=copy.deepcopy(mesh);result,summary=split(mesh,spec)
            self.assertEqual(mesh,before);self.assertEqual(summary['interfaces'],{'Seam':2})
            self.assertEqual(summary['added_displacement_nodes'],9 if quadratic else 4)
            self.assertEqual(result['material_ids'],mesh['material_ids']);self.assertEqual(result['region_ids'],mesh['region_ids'])
            np.testing.assert_array_equal(np.asarray(result['nodes'])[ref.cells(result)],np.asarray(mesh['nodes'])[ref.cells(mesh)])
            owners=ref.topology(result)
            for a,b in zip(result['face_sets']['Seam:a'],result['face_sets']['Seam:b']):
                self.assertEqual(len(owners[tuple(a)]),1);self.assertEqual(len(owners[tuple(b)]),1)
                na=ref.face_nodes(result,owners[tuple(a)][0],a);nb=ref.face_nodes(result,owners[tuple(b)][0],b)
                self.assertFalse(set(na)&set(nb))
                self.assertEqual({tuple(p) for p in np.asarray(result['nodes'])[na]},{tuple(p) for p in np.asarray(result['nodes'])[nb]})
            # Every original glued same-material face remains shared by two cells.
            self.assertEqual(sum(len(o)==2 for o in owners.values()),sum(len(o)==2 for o in ref.topology(mesh).values())-2)
            json.dumps(result,allow_nan=False)

    def test_missing_duplicate_and_partial_interface_declarations_are_rejected(self):
        root=Path(__file__).resolve().parents[1];mesh=json.loads((root/'examples/mc-column-mesh.json').read_text())
        spec=json.loads((root/'examples/joint-column-split.json').read_text())
        bad=copy.deepcopy(spec);bad['interfaces'][0]['body_b']='missing'
        with self.assertRaises(ValueError):split(mesh,bad)
        bad=copy.deepcopy(spec);bad['interfaces'].append(copy.deepcopy(bad['interfaces'][0]));bad['interfaces'][1]['id']='again'
        with self.assertRaises(ValueError):split(mesh,bad)
        # A single tetra on side b shares a crack front through undeclared glued
        # faces. Splitting only that patch must not tear the surrounding mesh.
        mesh['regions'].append({**mesh['regions'][1],'id':3,'body':'Patch'})
        mesh['region_ids'][6]=3;bad=copy.deepcopy(spec);bad['interfaces'][0]['body_b']='Patch'
        with self.assertRaisesRegex(ValueError,'tied crack fronts'):split(mesh,bad)


if __name__=='__main__':unittest.main()
