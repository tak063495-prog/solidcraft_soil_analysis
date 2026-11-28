"""Separate depth geometry, field diagnostics, and frozen MCC evidence regressions."""
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
import plastic_spatial_convergence as convergence
import diagnose_plastic_spatial as diagnostic
import summarize_mcc_spatial_records as mcc

BINARY=ROOT/'target/release'/('soil-fem-3d.exe' if os.name=='nt' else 'soil-fem-3d')
EVIDENCE=ROOT/'interim/2026-11-07-depth-spatial-diagnostics'
OLD=ROOT/'interim/2026-11-07-plastic-spatial-convergence'


class DepthGeometry(unittest.TestCase):
    def test_default_inputs_still_match_all_original_fingerprints(self):
        with tempfile.TemporaryDirectory() as temp:
            for finite,name in [(False,'plane'),(True,'finite')]:
                archived=convergence.read(OLD/f'{name}.json')
                for level in archived['levels']:
                    mesh,ground,case=convergence.inputs(level['divisions'],finite)
                    for key,value in [('mesh',mesh),('ground',ground),('case',case)]:
                        path=Path(temp)/f'{name}-{level["divisions"]}-{key}.json'
                        convergence.write(path,value)
                        self.assertEqual(convergence.sha(path),level['files_sha256'][key+'.json'])

    def test_depth_refinement_preserves_cut_volume_supports_and_in_plane_grid(self):
        for finite in (False,True):
            for ny in (2,4,8):
                mesh,_,case=convergence.inputs(2,finite,depth_divisions=ny)
                xyz=np.asarray(mesh['nodes']);tets=np.asarray(mesh['tetrahedra'])
                volume=np.linalg.det(np.transpose(xyz[tets[:,1:]]-xyz[tets[:,:1]],(0,2,1)))/6
                cut=np.asarray(mesh['region_ids'])==2
                self.assertEqual(len(tets),24*ny)
                self.assertTrue((volume>0).all())
                self.assertAlmostEqual(float(volume.sum()),1.)
                self.assertAlmostEqual(float(volume[cut].sum()),.125 if finite else .25)
                self.assertEqual(sorted(set(xyz[:,0])),[0.,.5,1.])
                self.assertEqual(sorted(set(xyz[:,2])),[0.,.5,1.])
                self.assertEqual(len(set(xyz[:,1])),ny+1)
                centers=xyz[tets].mean(axis=1)
                self.assertTrue((centers[cut,0]<.5).all())
                self.assertTrue((centers[cut,2]>.5).all())
                if finite:
                    self.assertTrue((centers[cut,1]<.5).all())
                names=[c['node_set'] for c in case['boundary']['constraints']]
                self.assertEqual(names,['sides','ends' if finite else 'all','bottom'])

    def test_invalid_or_mixed_axes_do_not_create_output_or_solve(self):
        with tempfile.TemporaryDirectory() as temp, patch.object(convergence.subprocess,'run') as solve:
            for options in (dict(axis='depth',depth_divisions=2),dict(axis='depth',divisions=[2,4,16]),
                            dict(axis='depth',in_plane_divisions=3),dict(depth_divisions=3),
                            dict(finite_width=True,depth_divisions=1),dict(axis='depth',divisions=[32,64,128])):
                path=Path(temp)/'invalid'
                with self.subTest(options=options),self.assertRaises(ValueError):
                    convergence.run_series(path,BINARY,**options)
                self.assertFalse(path.exists())
            solve.assert_not_called()


class TensorDiagnostic(unittest.TestCase):
    def test_constant_tensor_energy_decompositions_on_depth_meshes(self):
        ma,_,_=convergence.inputs(2,depth_divisions=2)
        mb,_,_=convergence.inputs(2,depth_divisions=4)
        def manufactured(mesh,stress):
            return dict(frames=[dict(name='cut',state=dict(displacement_m=np.zeros_like(mesh['nodes']).tolist(),
                active_elements=[True]*len(mesh['tetrahedra']),material_states=[dict(
                    effective_stress_compression_positive_kpa=stress,equivalent_plastic_strain=.1)]*len(mesh['tetrahedra'])))])
        # diag(1,2,3), xy=4: norm² = 14+32=46, hydro=12, dev=34.
        pair=convergence.compare(ma,manufactured(ma,[0.]*6),mb,manufactured(mb,[1.,2.,3.,4.,0.,0.]),diagnostics=True)['stages'][0]
        d=pair['stress_diagnostic']
        self.assertAlmostEqual(pair['stress']['rms_difference_kpa']**2,46.)
        np.testing.assert_allclose(d['component_squared_integrals_kpa2_m3'],[1.,4.,9.,32.,0.,0.],atol=1e-11)
        self.assertAlmostEqual(d['hydrostatic_squared_integral_kpa2_m3'],12.)
        self.assertAlmostEqual(d['deviatoric_squared_integral_kpa2_m3'],34.)
        self.assertAlmostEqual(sum(v['squared_integral_kpa2_m3'] for v in d['plastic_partition']),46.)
        self.assertAlmostEqual(sum(v['volume_m3'] for v in d['plastic_partition']),1.)

    def test_archived_indeterminate_diagnosis_and_control_oscillation_persist(self):
        with tempfile.TemporaryDirectory() as temp:
            report=diagnostic.diagnose(EVIDENCE/'plane-inplane-ny1.json',OLD/'plane.json',OLD/'plane-half.json',Path(temp)/'diagnostic.json')
        self.assertAlmostEqual(report['stress']['recorded_ratio'],.7110377721615777)
        self.assertIsNone(report['stress']['estimate'])
        monitor=report['monitor']
        self.assertIsNone(monitor['original_richardson']['estimate'])
        self.assertIsNone(monitor['half_control_richardson']['estimate'])
        self.assertLess(monitor['construction_signed_changes_m'][0],0)
        self.assertGreater(monitor['construction_signed_changes_m'][1],0)
        for delta,control in zip(monitor['construction_signed_changes_m'],monitor['half_control_pair_change_bounds_m']):
            self.assertGreater(abs(delta),control)
        self.assertFalse(report['p6_complete'])
        self.assertFalse(report['fs']['nonconvergence_is_upper_bound'])

    def test_diagnostic_refuses_mixed_depth_and_changed_input(self):
        for mutation in ('depth','input'):
            report=convergence.read(EVIDENCE/'plane-inplane-ny1.json')
            if mutation=='depth':
                report['levels'][0]['depth_divisions']=2
            else:
                report['levels'][0]['files_sha256']['mesh.json']='changed'
            with tempfile.TemporaryDirectory() as temp:
                path=Path(temp)/'changed.json';convergence.write(path,report)
                with self.assertRaises(ValueError):
                    diagnostic.diagnose(path,OLD/'plane.json',OLD/'plane-half.json',Path(temp)/'never.json')
                self.assertFalse((Path(temp)/'never.json').exists())


class SmallDepthSolves(unittest.TestCase):
    def test_two_separate_depth_series_and_finest_half_control(self):
        with tempfile.TemporaryDirectory() as temp:
            for finite in (False,True):
                folder=Path(temp)/str(finite)
                report=convergence.run_series(folder,BINARY,axis='depth',finite_width=finite)
                self.assertTrue(report['complete'],report['failures'])
                self.assertEqual(report['refinement_axis'],'depth')
                self.assertEqual([v['elements'] for v in report['levels']],[48,96,192])
                self.assertEqual([v['h_xz_m'] for v in report['levels']],[.5]*3)
                self.assertEqual([v['depth_divisions'] for v in report['levels']],[2,4,8])
                self.assertFalse(report['general_3d_spatial_convergence_demonstrated'])
                self.assertFalse(report['p6_complete'])
                for v in report['levels']:
                    self.assertAlmostEqual(v['frames'][-1]['active_volume_m3'],.875 if finite else .75)
                    self.assertGreater(v['frames'][-1]['active_plastic_volume_m3'],0)
                with self.assertRaises(FileExistsError):
                    convergence.run_series(folder,BINARY,axis='depth',finite_width=finite)
                if finite:
                    order=next(v for v in report['field_difference_orders'] if v['stage']=='cut' and v['field']=='stress')
                    self.assertIsNone(order['order'])
                    mesh,ground,case=convergence.inputs(2,True,1/32,5e-5,depth_divisions=8)
                    half=Path(temp)/'half';half.mkdir()
                    for key,value in [('mesh',mesh),('ground',ground),('case',case)]:
                        convergence.write(half/f'{key}.json',value)
                    solve=subprocess.run([str(BINARY),'solve-stages',str(half/'mesh.json'),str(half/'ground.json'),
                        str(half/'case.json'),str(half/'result.json')],capture_output=True,text=True,timeout=120)
                    self.assertEqual(solve.returncode,0,solve.stderr)
                    comparison=convergence.compare(mesh,convergence.read(folder/'ny8/result.json'),
                        mesh,convergence.read(half/'result.json'))['stages'][-1]
                    self.assertLess(comparison['stress']['rms_difference_kpa'],order['difference_norms'][-1])

    def test_timeout_stays_incomplete_without_order_or_upper_factor(self):
        with tempfile.TemporaryDirectory() as temp, patch.object(convergence.subprocess,'run',side_effect=subprocess.TimeoutExpired('solve',1)):
            report=convergence.run_series(Path(temp)/'timed-out',BINARY,axis='depth',timeout=1)
        self.assertFalse(report['complete'])
        self.assertEqual(report['failures'][0]['reason'],'wall_time_limit')
        self.assertEqual(report['field_difference_orders'],[])
        self.assertEqual(report['fs']['status'],'not_attempted')
        self.assertFalse(report['p6_complete'])


class FrozenMccEvidence(unittest.TestCase):
    def test_existing_increment_records_are_separate_and_not_accepted(self):
        report=mcc.collect()
        self.assertEqual(report['new_solves'],0)
        self.assertEqual(report['new_field_integrations'],0)
        self.assertTrue(report['fixed_wall_5144_8061_14909_increment_difference_non_decreasing'])
        rows=report['rows']
        np.testing.assert_allclose([r['frames'][-1]['stage_increment']['relative_l2'] for r in rows[:4]],
            [.181664468373731,.1278766023138621,.0996904476826082,.06368386176664693],rtol=1e-14)
        self.assertTrue(all(r['wall_fixed'] for r in rows[3:]))
        self.assertTrue(all(not r['spatial_convergence_accepted'] for r in rows))
        self.assertIsNone(report['order'])
        self.assertFalse(report['p6_complete'])

    def test_missing_or_modified_frozen_report_never_falls_back_to_solver(self):
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaises(FileNotFoundError):
                mcc.collect(Path(temp))
            path=Path(temp)/mcc.RECORDS[0][1];path.parent.mkdir(parents=True)
            shutil.copyfile(ROOT/mcc.RECORDS[0][1],path)
            path.write_bytes(path.read_bytes()+b' ')
            with self.assertRaisesRegex(ValueError,'frozen MCC comparison changed'):
                mcc.collect(Path(temp))


if __name__=='__main__':
    unittest.main()
