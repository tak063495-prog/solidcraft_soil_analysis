"""New offset scope: independent work map, guards, VTU and full dry lifecycle."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
import meshio
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from plate6_offset_reference import Element,System
from plate6_offset_study import inputs,study,command,write,semantic
from plate6_reference import Element as Planar
from beam_reference import System as OldSystem
from plate6_scope_rejections import probe

ROOT=Path(__file__).resolve().parents[1]
BINARY=ROOT/'target/release'/('soil-fem-3d.exe' if os.name=='nt' else 'soil-fem-3d')

class Plate6Offset(unittest.TestCase):
    def element(self,e):
        p=np.array([[0.,0.,0.],[2.,0.,0.],[.4,1.3,0.],[1.,0.,0.],[1.2,.65,0.],[.2,.65,0.]])
        section=dict(young_kpa=2e8,poisson=.2,thickness_m=.2,drilling_factor=1e-4)
        return p,section,Element(p,section,e)

    def test_independent_offset_rank_rigid_modes_and_trace(self):
        for e in (-.1,.1):
            p,_,k=self.element(e);v=np.linalg.eigvalsh(k.stiffness)
            self.assertEqual(np.count_nonzero(v>1e-10*v.max()),30)
            for a in range(6):
                d=np.zeros((6,6))
                if a<3:d[:,a]=1.
                else:d[:,3+a-3]=1.;d[:,:3]=np.cross(np.eye(3)[a-3],p)
                for sample in k.recover(d.ravel())['samples']:
                    for key in ('membrane_strain','curvature_per_m','transverse_shear_strain'):
                        np.testing.assert_allclose(sample[key],0.,atol=1e-13,rtol=0.)
            d=np.arange(36).reshape(6,6)*1e-5
            for l in ([.2,.3,.5],[0.,.37,.63],[.23,0.,.77],[.61,.39,0.]):
                n,*_=k.neutral.shapes(*(np.array(l)@k.neutral.xy[:3]))
                np.testing.assert_allclose(k.attachment_translation(l,d.ravel()),n@d[:,:3],atol=1e-10,rtol=0.)

    def test_independent_work_map_and_zero_exactness(self):
        p,s,k=self.element(0.);old=Planar(p,s);np.testing.assert_array_equal(k.stiffness,old.stiffness)
        d=np.arange(36)*1e-6
        for e in (-.1,.1):
            _,_,k=self.element(e);r=k.recover(d)
            total=sum(r[key] for key in ('membrane_energy_kn_m','bending_energy_kn_m','transverse_shear_energy_kn_m','drilling_penalty_energy_kn_m'))
            np.testing.assert_allclose(total,d@k.stiffness@d/2,atol=1e-8,rtol=1e-7)
            f=k.load([1.,2.,3.]);self.assertAlmostEqual(float(f@d),float(k.neutral.load([1.,2.,3.])@(k.map@d)),places=14)

    def test_old_reference_does_not_silently_accept_offset(self):
        m,_,c,_=inputs()
        with self.assertRaises(ValueError):OldSystem(m,c['beam_system'])
        c=copy.deepcopy(c)
        for p in c['beam_system']['plates']:p['formulation']='mindlin_tri6'
        with self.assertRaisesRegex(ValueError,'separate'):OldSystem(m,c['beam_system'])
        with self.assertRaises(ValueError):System(m,c['beam_system'])

    def test_native_preparation_keeps_offset_scope_explicitly_refused(self):
        from prepare_native_plates import derive
        m,_,c,_=inputs()
        before=copy.deepcopy((m,c))
        with self.assertRaisesRegex(ValueError,'Native offset'):derive(m,c)
        self.assertEqual((m,c),before)

    def test_zero_offset_complete_physical_run_exactly_matches_planar(self):
        m,g,c,_=inputs(0.)
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            for name,value in [('mesh',m),('ground',g),('stages',c)]:write(p/(name+'.json'),value)
            command(BINARY,['solve-stages',p/'mesh.json',p/'ground.json',p/'stages.json',p/'offset.json'])
            for plate in c['beam_system']['plates']:plate['formulation']='mindlin_tri6';plate.pop('attachment_offset_m')
            write(p/'planar-stages.json',c)
            command(BINARY,['solve-stages',p/'mesh.json',p/'ground.json',p/'planar-stages.json',p/'planar.json'])
            a,b=(json.loads((p/n).read_text()) for n in ('offset.json','planar.json'))
            self.assertTrue(a['complete'] and b['complete'])
            for key in ('frames','last_committed','attempts','geostatic','transition_audits'):
                self.assertEqual(semantic(a[key]),semantic(b[key]),key)

    def test_curves_and_unverified_combinations_remain_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            r=probe(Path(tmp)/'refused',BINARY)
            self.assertTrue(r['refusal_contract_passed']);self.assertFalse(r['curved_tri6_accepted'])
            self.assertTrue(all(not b['passed'] for b in r['classical_benchmarks']))

    def test_positive_and_negative_offset_direct_system_restart_and_quadratic_vtu(self):
        with tempfile.TemporaryDirectory() as tmp:
            for e in (-.05,.05):
                output=Path(tmp)/str(e);r=study(output,e=e)
                self.assertTrue(r['passed']);self.assertTrue(r['accepted_run_exact_comparison'])
                self.assertLess(r['maximum_attachment_trace_error_m'],1e-10)
                self.assertEqual(len(r['separate_process_restart_stops']),5)
                self.assertFalse(r['p6_complete'])
                vtu=meshio.read(output/'plate-output/plates-002.vtu')
                self.assertIn('triangle6',vtu.cells_dict)
                np.testing.assert_allclose(vtu.cell_data['attachment_offset_m'][0],e,atol=0.,rtol=0.)
                np.testing.assert_allclose(vtu.points[:,2]-vtu.point_data['neutral_reference_position_m'][:,2],e,atol=1e-14,rtol=0.)
                self.assertEqual(json.loads((output/'plate-output/plate-fields.json').read_text())['schema'],'soil-fem-plate-export/3')

if __name__=='__main__':unittest.main()
