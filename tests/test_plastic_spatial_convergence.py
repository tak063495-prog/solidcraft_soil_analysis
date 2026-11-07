"""Physical h-series, exact overlap counterexamples, and unchanged-solver regressions."""
import copy
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
import plastic_spatial_convergence as convergence
import prepare_slope_h_series as slope
import assess_slope_h_series as slope_assessment

BINARY=ROOT/'target/release'/('soil-fem-3d.exe' if os.name=='nt' else 'soil-fem-3d')


class ReferenceFreeOrders(unittest.TestCase):
    def test_known_second_order_recovers_limit_without_reference(self):
        report=convergence.richardson([3.25,3.0625,3.015625])
        self.assertAlmostEqual(report['order'],2.)
        self.assertAlmostEqual(report['estimate'],3.)
        self.assertAlmostEqual(report['finest_discretization_component'],.015625)

    def test_oscillation_growth_floor_and_two_levels_do_not_invent_limits(self):
        for values in ([1.,2.,1.9],[1.,2.,4.],[1.,1.,1.],[1.,1.+1e-12,1.+1.25e-12],[1.,1.5]):
            with self.subTest(values=values):
                self.assertIsNone(convergence.richardson(values)['estimate'])
        with self.assertRaises(ValueError):
            convergence.richardson([1.,float('nan'),2.])


class ExactFields(unittest.TestCase):
    @staticmethod
    def manufactured(mesh, shift=0., shear=0.):
        u=np.asarray(mesh['nodes'])*[.01,.02,.03]+shift
        state=dict(displacement_m=u.tolist(),active_elements=[True]*len(mesh['tetrahedra']),
                   material_states=[dict(effective_stress_compression_positive_kpa=[1.,2.,3.,shear,0.,0.])
                                    for _ in mesh['tetrahedra']])
        return dict(frames=[dict(name='same',state=state)])

    def test_affine_displacement_and_constant_tensor_norms_on_nonmatching_tets(self):
        ma,_,_=convergence.inputs(2); mb,_,_=convergence.inputs(4)
        a=self.manufactured(ma); b=self.manufactured(mb,.001,2.)
        result=convergence.compare(ma,a,mb,b)['stages'][0]
        self.assertAlmostEqual(result['displacement']['rms_difference_m'],np.sqrt(3)*.001,places=12)
        self.assertAlmostEqual(result['stress']['rms_difference_kpa'],np.sqrt(8),places=12)
        self.assertLess(result['displacement']['maximum_element_coverage_relative_error'],1e-8)

    def test_unequal_active_domains_are_rejected(self):
        mesh,_,_=convergence.inputs(2)
        a=self.manufactured(mesh); b=copy.deepcopy(a)
        b['frames'][0]['state']['active_elements'][0]=False
        with self.assertRaisesRegex(ValueError,'same active physical domains'):
            convergence.compare(mesh,a,mesh,b)


class PlasticExcavation(unittest.TestCase):
    def test_ci_in_plane_and_finite_width_fields_and_half_increment(self):
        self.assertTrue(BINARY.is_file())
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp)
            for finite in (False,True):
                report=convergence.run_series(folder/str(finite),BINARY,finite_width=finite)
                self.assertTrue(report['complete'],report['failures'])
                self.assertFalse(report['p6_complete'])
                self.assertFalse(report['general_3d_spatial_convergence_demonstrated'])
                self.assertFalse(report['griffiths_slope_five_percent_closed'])
                self.assertEqual(report['fs']['status'],'not_attempted')
                self.assertEqual([v['elements'] for v in report['levels']],
                                 [48,192,768] if finite else [24,96,384])
                for level in report['levels']:
                    cut=level['frames'][-1]
                    self.assertGreater(cut['active_plastic_volume_m3'],0)
                    self.assertAlmostEqual(cut['active_volume_m3'],.875 if finite else .75,places=11)
                    self.assertLessEqual(cut['relative_force_residual'],1e-8)
                    self.assertLessEqual(cut['global_force_relative_error'],1e-8)
                    self.assertLessEqual(cut['maximum_normalized_yield'],1e-7)
                    self.assertLess(cut['maximum_engineering_strain'],.05)
                orders={v['field']:v for v in report['field_difference_orders'] if v['stage']=='cut'}
                for field in ('displacement','construction_displacement','stress'):
                    a,b=orders[field]['difference_norms']
                    self.assertLess(b,a)
                # Do not turn the plane-strain stress's insufficient reduction into a pass.
                if not finite:
                    self.assertIsNone(orders['stress']['order'])
                # Solve the finest mesh with a half load increment and the same material substep.
                mesh,ground,case=convergence.inputs(8,finite,1/32,5e-5)
                sub=folder/f'half-{finite}';sub.mkdir()
                for name,value in [('mesh',mesh),('ground',ground),('case',case)]:
                    convergence.write(sub/f'{name}.json',value)
                process=subprocess.run([str(BINARY),'solve-stages',str(sub/'mesh.json'),str(sub/'ground.json'),
                    str(sub/'case.json'),str(sub/'result.json')],capture_output=True,text=True,timeout=120)
                self.assertEqual(process.returncode,0,process.stderr)
                half=convergence.read(sub/'result.json');self.assertTrue(half['complete'])
                normal=convergence.read(folder/str(finite)/'n8/result.json')
                comparison=convergence.compare(mesh,normal,mesh,half)['stages'][-1]
                # Existing spatial difference must dominate this numerical-control change.
                for field,key in [('displacement','rms_difference_m'),('stress','rms_difference_kpa')]:
                    self.assertLess(comparison[field][key],orders[field]['difference_norms'][-1])

    def test_invalid_series_and_existing_destination_are_refused(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)
            with self.assertRaises(FileExistsError):
                convergence.run_series(path,BINARY)
            with self.assertRaises(ValueError):
                convergence.run_series(path/'bad',BINARY,divisions=[2,8,16])
            self.assertFalse((path/'bad').exists())


class SlopePreparation(unittest.TestCase):
    def test_assessment_of_known_constant_prism_stress_preserves_missing_reference(self):
        with tempfile.TemporaryDirectory() as temp:
            out=Path(temp)/'slope';report=slope.prepare(out)
            for level in report['levels']:
                sub=out/level['folder'];mesh=convergence.read(sub/'mesh.json')
                xyz=np.asarray(mesh['nodes']);tet=np.asarray(mesh['tetrahedra'])
                volumes=np.linalg.det(np.transpose(xyz[tet[:,1:]]-xyz[tet[:,:1]],(0,2,1)))/6
                result=dict(complete=True,ground=convergence.read(sub/'ground.json'),case=convergence.read(sub/'case.json'),
                    provenance=dict(solver_binary_sha256='manufactured assessment, not a solver run'),
                    frames=[dict(relative_force_residual=0.,maximum_normalized_yield=0.,
                        maximum_absolute_engineering_strain_from_seed=0.,state=dict(
                            displacement_m=np.zeros_like(xyz).tolist(),active_elements=[True]*len(tet),
                            material_states=[dict(effective_stress_compression_positive_kpa=[3.,3.,3.,0.,0.,0.])]*len(tet)))],
                    integration_points=[dict(element_index=i,weight_m3=float(v)) for i,v in enumerate(volumes)])
                convergence.write(sub/'result.json',result)
            assessed=slope_assessment.assess(out,Path(temp)/'assessment.json')
            self.assertFalse(assessed['p6_complete'])
            self.assertFalse(assessed['griffiths_slope_five_percent_closed'])
            self.assertIsNone(assessed['monitor']['estimate'])
            for pair in assessed['pairs']:
                self.assertLess(pair['source_prism_mean_stress_rms_difference_kpa'],1e-12)
            changed=out/report['levels'][0]['folder']/'case.json'
            changed.write_text('{}')
            with self.assertRaisesRegex(ValueError,'prepared input changed'):
                slope_assessment.assess(out,Path(temp)/'rejected.json')

    def test_pinned_geometry_fixed_depth_and_preserved_supports_without_displacements(self):
        geometry=convergence.read(slope.GEOMETRY)
        self.assertEqual(geometry['source_revision'],'440d246d9e32f5892c79539d4a7b797ad29442ee')
        self.assertNotIn('displacements',geometry)
        with tempfile.TemporaryDirectory() as temp:
            out=Path(temp)/'slope';report=slope.prepare(out)
            self.assertFalse(report['solved'])
            self.assertFalse(report['griffiths_slope_five_percent_closed'])
            self.assertEqual([v['elements'] for v in report['levels']],[1200,4800,19200])
            self.assertEqual([v['depth_divisions'] for v in report['levels']],[1,1,1])
            mesh=convergence.read(out/'xy1-ny1/mesh.json')
            xyz=np.asarray(mesh['nodes'])
            volumes=np.linalg.det(np.transpose(xyz[np.asarray(mesh['tetrahedra'])[:,1:]]-
                xyz[np.asarray(mesh['tetrahedra'])[:,:1]],(0,2,1)))/6
            self.assertGreater(volumes.min(),0)
            self.assertAlmostEqual(float(volumes.sum()),220.,places=10)
            self.assertEqual(len({tuple(p) for p in xyz}),len(xyz))
            result=subprocess.run([str(BINARY),'validate',str(out/'xy1-ny1/mesh.json')],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            cases=[convergence.read(out/v['folder']/'case.json') for v in report['levels']]
            self.assertTrue(all(case==cases[0] for case in cases))


if __name__=='__main__':
    unittest.main()
