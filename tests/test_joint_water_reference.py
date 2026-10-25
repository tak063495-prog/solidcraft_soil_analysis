"""Closed-form gates for the independent hydraulic/pressure reference."""
import json
from pathlib import Path
import sys
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tools'))
from joint_water_reference import Water

def read(n):return json.loads((ROOT/'examples'/n).read_text(encoding='utf-8'))

class Reference(unittest.TestCase):
    def test_series_head_and_flux_with_separated_mechanical_nodes(self):
        for kind in ('mesh','tet10-mesh'):
            mesh=read(f'joint-water-column-{kind}.json');ground=read('ground-joint-water-mc.json');case=read('stages-joint-water-flow.json');w=Water(mesh,ground,case);active=[True]*len(mesh['tetrahedra']);r=w.steady(case['stages'][0]['flow'],active,mesh['material_ids'],.75)
            for e,g in r['solution']['quadrature'].items():
                h=r['solution']['relative_head'][g['nodes']];grad=np.stack([g['b'][:,a,a::3]@(h-h[0]) for a in range(3)],axis=1);assert np.max(np.abs(grad[:,:2]))<1e-11
                conductivity=1e-4 if mesh['material_ids'][e]==1 else 2e-4;self.assertLess(np.max(np.abs(-conductivity*grad[:,2]+1/7500)),1e-13)
            interface=w.from_heads(active,r['field'],.75);self.assertLess(np.max(np.abs(interface['pressure']-10*(2+2/3-.5))),1e-10)
            self.assertTrue(all(a!=b for p in w.contact.points for a,b in zip(p['nodes_a'],p['nodes_b'])))

    def test_prescribed_interface_is_a_reservoir_and_structure_is_not_aquifer(self):
        mesh=read('joint-water-column-tet10-mesh.json');ground=read('ground-joint-water-mc.json');active=[True]*len(mesh['tetrahedra'])
        for name,expected in (('soil',15.),('prescribed',2.5)):
            case=read(f'stages-joint-water-{name}.json');w=Water(mesh,ground,case);r=w.steady(case['stages'][0]['flow'],active,mesh['material_ids'],.75);interface=w.from_heads(active,r['field'],.75);self.assertLess(np.max(np.abs(interface['pressure']-expected)),1e-10)
            if name=='soil':self.assertEqual(sum(r['soil_active']),len(active)//2)

    def test_partial_wet_sloped_interface_matches_exact_projected_resultant(self):
        mesh=read('joint-water-column-mesh.json');ground=read('ground-joint-water-mc.json');case=read('stages-joint-water-column.json')
        for p in mesh['nodes']:
            if abs(p[2]-.5)<1e-12:p[2]=.2+.6*p[0]
        w=Water(mesh,ground,case);field=w.hydrostatic([True]*len(mesh['tetrahedra']),.5);nodes=np.unique(np.asarray(mesh['tetrahedra'])[np.asarray(mesh['region_ids'])==1]);self.assertTrue(np.allclose(field['load'][nodes].sum(axis=0),[.45,0.,-.75],rtol=0.,atol=1e-12));self.assertLess(np.linalg.norm(field['load'].sum(axis=0)),1e-12)

if __name__=='__main__':unittest.main()
