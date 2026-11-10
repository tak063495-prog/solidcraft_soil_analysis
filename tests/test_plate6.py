"""Quadratic plate trace, independent polynomial mechanics and native exports."""
import os
for name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[name]='1'
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
import meshio
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from plate6_reference import Element
from plate6_study import inputs,study
import plate_surfaces
import plate_io

class Plate6(unittest.TestCase):
    def element(self):
        p=np.array([[0.,0,0],[2.,0,0],[.4,1.3,0],[1.,0,0],[1.2,.65,0],[.2,.65,0]])
        s=dict(young_kpa=2e8,poisson=.2,thickness_m=.2,drilling_factor=1e-4)
        return p,Element(p,s)

    def test_explicit_legacy_selection_matches_rust_canonical_input(self):
        import fem_reference as ref
        root=Path(__file__).resolve().parents[1]
        original=json.loads((root/'examples/stages-plate-beam-column.json').read_text())
        explicit=copy.deepcopy(original)
        for plate in explicit['beam_system']['plates']:plate['formulation']='cst_dkt'
        self.assertEqual(ref.normalize_stage_case(explicit),ref.normalize_stage_case(original))

    def test_independent_rank_and_constant_shear_energy(self):
        p,e=self.element();v=np.linalg.eigvalsh(e.stiffness)
        self.assertEqual(np.count_nonzero(v>v.max()*1e-10),30)
        d=np.zeros((6,6));d[:,2]=.003*p[:,0]-.002*p[:,1]
        result=e.recover(d.ravel())
        for s in result['samples']:np.testing.assert_allclose(s['transverse_shear_strain'],[.003,-.002],atol=1e-14,rtol=0.)
        self.assertAlmostEqual(result['transverse_shear_energy_kn_m'],d.ravel()@e.stiffness@d.ravel()/2,places=9)

    def test_consistent_load_rejects_equal_six_virtual_work(self):
        _,e=self.element();f=e.load([0.,0.,1.]).reshape(6,6)[:,2]
        np.testing.assert_allclose(f[:3],0.,atol=1e-14,rtol=0.)
        np.testing.assert_allclose(f[3:],e.area/3,atol=1e-14,rtol=0.)
        variation=np.array([1.,0,0,.25,0,.25]);exact=e.area/6
        self.assertAlmostEqual(f@variation,exact,places=14)
        self.assertGreater(abs(np.full(6,e.area/6)@variation-exact),.01)

    def test_surface_default_retains_four_split_and_high_order_is_explicit(self):
        mesh,_,case,binding=inputs()
        self.assertEqual(len(case['beam_system']['plates']),2)
        self.assertTrue(all(len(p['nodes'])==6 and p['formulation']=='mindlin_tri6' for p in case['beam_system']['plates']))
        self.assertIn('same straight planar Tri6',binding['interpolation_limitation'])
        raw=copy.deepcopy(case);raw.pop('beam_system')
        mesh=copy.deepcopy(mesh);mesh['node_sets']={k:v for k,v in mesh['node_sets'].items() if not k.startswith('plate:mesh-node:')}
        body=mesh['regions'][-1]['body']
        spec=dict(schema='soil-fem-plate-surfaces/1',sections=[case['beam_system']['plate_sections'][0]],surfaces=[dict(id='p',face_set='top',owner_body=body,section_id='wall_section',initially_active=True)])
        _,legacy,_=plate_surfaces.attach(mesh,raw,spec)
        self.assertEqual(len(legacy['beam_system']['plates']),8)
        self.assertTrue(all(len(p['nodes'])==3 and 'formulation' not in p for p in legacy['beam_system']['plates']))
        for field,value in [('formulation','unknown'),('thickness_offset_m',.01),('load_distribution','equal_six')]:
            bad=copy.deepcopy(spec);bad['surfaces'][0][field]=value
            with self.subTest(field=field),self.assertRaises(ValueError):plate_surfaces.attach(mesh,raw,bad)

    def test_quadratic_vtu_roundtrip_and_mixed_legacy_contract(self):
        p,e=self.element();d=np.arange(36)*1e-6;r=e.recover(d);r['id']='high'
        # Mixed legacy/high-order output, with unequal sample counts.
        from plate_reference import Element as Legacy
        lo=Legacy(p[:3],dict(young_kpa=2e8,poisson=.2,thickness_m=.2,drilling_factor=1e-4)).recover(d[:18]);lo['id']='low'
        cfg=dict(nodes=[dict(id=str(i)) for i in range(6)],plates=[dict(id='high',nodes=list(map(str,range(6))),formulation='mindlin_tri6'),dict(id='low',nodes=['0','1','2'])])
        frame=dict(name='both',beam_results=dict(plates=[r,lo],node_positions_m=p),state=dict(beams=dict(active_plates=[True,True],node_displacement_m=d.reshape(6,6)[:,:3],node_rotation_rad=d.reshape(6,6)[:,3:])))
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp);plate_io.export(dict(case=dict(beam_system=cfg),frames=[frame]),path)
            mesh=meshio.read(path/'plates-000.vtu');self.assertEqual(set(mesh.cells_dict),{'triangle','triangle6'})
            self.assertEqual(mesh.cells_dict['triangle6'].tolist(),[list(range(6))])
            samples=meshio.read(path/'plate-samples-000.vtu');self.assertEqual(len(samples.points),9)
            self.assertIn('transverse_shear_force_kn_m',samples.point_data)
            self.assertEqual(json.loads((path/'plate-fields.json').read_text())['schema'],'soil-fem-plate-export/2')

    def test_small_coupled_study_direct_reference_restarts_and_falsification(self):
        with tempfile.TemporaryDirectory() as tmp:
            r=study(Path(tmp)/'study');self.assertTrue(r['passed']);self.assertFalse(r['p6_complete'])
            self.assertEqual(r['frames'],5);self.assertEqual(len(r['separate_process_restart_stops']),5)

if __name__=='__main__':unittest.main()
