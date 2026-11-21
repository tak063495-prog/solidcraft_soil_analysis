"""Analytic falsification and short moving-plane paths; no NC acceptance."""
import copy
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
import fold_stability as fold


def legacy():
    folder=ROOT/'interim/2026-11-09-nc-beyond-fold'
    sys.path.insert(0,str(folder.parent/'2026-11-09-nc-pseudo-arclength'))
    spec=importlib.util.spec_from_file_location('fold_existing_continuation',folder/'continue_path.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


from nc_fold_analytic import Analytic
from nc_fold_path import control_binding
from analyze_nc_fold import linear_audit

class FoldChecks(unittest.TestCase):
    def test_true_nonsymmetric_left_and_right_modes_at_exact_fold(self):
        k=np.array([[0.,0.],[.7,1.]])
        d=fold.point(k,[1.,0.],[1/np.sqrt(1.49),-.7/np.sqrt(1.49),0.],
                     contributions=[dict(kind='gp',matrix=k)])
        self.assertAlmostEqual(d['stiffness']['minimum'],0.,places=14)
        self.assertEqual(d['stiffness']['numerical_rank'],1)
        self.assertAlmostEqual(d['normalized_transversality'],1.)
        self.assertEqual(d['augmented']['numerical_rank'],3)
        self.assertGreater(d['stiffness_asymmetry'],0.)
        self.assertLess(d['right_mode_residual'],1e-14)
        self.assertLess(d['left_mode_residual'],1e-14)
        self.assertLess(d['contribution_sum_relative_error'],1e-14)
        np.testing.assert_allclose(k@d['right_weak_mode'],[0.,0.],atol=1e-14)
        np.testing.assert_allclose(k.T@d['left_weak_mode'],[0.,0.],atol=1e-14)

    def path(self,xs):
        previous=np.array([1.,-.7,-2*xs[0]]);previous/=np.linalg.norm(previous);points=[]
        for x in xs:
            d=fold.point([[2*x,0.],[.7,1.]],[1.,0.],previous)
            points.append(dict(accepted=True,**{'lambda':.95-x*x},diagnostic=d,
                               relative_force_residual=0.,force_block_linear_residual=0.,augmented_linear_residual=0.))
            previous=np.array(d['scaled_tangent'])
        return points

    def test_known_fold_detected_and_regular_segment_not_detected(self):
        self.assertEqual(fold.judge(self.path([-.02,-.01,.01,.02]))['conclusion'],'local_fold_bracket')
        self.assertEqual(fold.judge(self.path([-.04,-.03,-.02]))['conclusion'],'undetermined')

    def test_short_rejected_and_failed_force_or_linear_paths_do_not_certify(self):
        points=self.path([-.02,-.01,.01,.02])
        self.assertEqual(fold.judge(points[:2])['conclusion'],'undetermined')
        for key,value in [('relative_force_residual',1.01e-8),('force_block_linear_residual',1.01e-10),
                          ('augmented_linear_residual',1.01e-10),('accepted',False)]:
            changed=copy.deepcopy(points);changed[1][key]=value
            self.assertEqual(fold.judge(changed)['conclusion'],'undetermined')

    def test_every_load_component_including_restrained_component_is_checked(self):
        final=np.array([1.,2.,3.]);h=np.array([3.,-1.,0.]);actual=final+.2*h
        fold.audit_load(actual,final,h,.8)
        actual[2]+=.1
        with self.assertRaisesRegex(ValueError,'component load'):fold.audit_load(actual,final,h,.8)

    def test_low_stress_captures_removed_points_and_never_floors_pressure(self):
        points=[dict(point=i,active=i<2,stress_kpa=[p,p,p,0,0,0],void_ratio=.8,kappa=.04,
                     position_m=[i,0,0],weak_mode_strain_norm=float(3-i)) for i,p in enumerate([0.,1e-8,100.])]
        result=fold.low_stress(points)
        self.assertEqual(result['all']['pressure_kpa']['count'],3)
        self.assertEqual(result['active']['pressure_kpa']['count'],2)
        self.assertEqual(result['points'][0]['bulk_modulus_kpa'],0.)
        self.assertAlmostEqual(result['points'][1]['bulk_modulus_kpa'],4.5e-7)
        self.assertEqual(result['active']['near_zero_pressure'][-1]['point_ids'],[0,1])

    def test_invalid_arrays_scales_and_unoriented_exact_fold_are_rejected(self):
        for k,h in [([[float('nan')]], [1.]),([[1.,0.]],[1.]),([[1.]],[0.])]:
            with self.assertRaises(ValueError):fold.point(k,h)
        with self.assertRaises(ValueError):fold.point([[0.]],[1.])
        with self.assertRaises(ValueError):fold.point([[1.]],[1.],displacement_scale_m=0.)

    def test_existing_moving_plane_solver_crosses_fold_and_actual_restart_matches(self):
        module=legacy();binding=control_binding(dict(test='analytic fixed equations'))
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            for name in ('whole','first','second'): (root/name).mkdir()
            whole=module.follow(Analytic(),root/'whole',.005,12,30.,binding)
            module.follow(Analytic(),root/'first',.005,6,30.,binding)
            module.follow(Analytic(),root/'second',.005,6,30.,binding,root/'first/restart.json')
            rows=lambda name:[json.loads(line) for line in (root/name/'samples.jsonl').read_text().splitlines()]
            direct=rows('whole');self.assertEqual(direct,rows('first')+rows('second'))
            self.assertTrue(any(r['tangent'][-1]>0 for r in direct));self.assertTrue(any(r['tangent'][-1]<0 for r in direct))
            self.assertFalse(whole['passed']);self.assertEqual(whole['stop_reason'],'bounded step count; inconclusive')
            payload=json.loads((root/'first/restart.json').read_text());payload['payload']['progress']+=.1
            (root/'tampered.json').write_text(json.dumps(payload))
            with self.assertRaisesRegex(ValueError,'payload differs'):module.restore(root/'tampered.json',binding)
            with self.assertRaisesRegex(ValueError,'binding'):module.restore(root/'first/restart.json',binding|dict(arc_lambda_scale_m=2.))

    def test_existing_solver_rejection_keeps_short_origin_and_acceptance_false(self):
        module=legacy();m=Analytic();before=m.old.copy()
        original=m.evaluate
        def reject(du,old,label):
            if label!='seed-zero':raise ValueError('intentional rejected trial')
            return original(du,old,label)
        m.evaluate=reject
        with tempfile.TemporaryDirectory() as directory:
            report=module.follow(m,Path(directory),.005,3,30.,dict(test='rejection'))
            self.assertEqual(report['equilibrium_samples'],0);self.assertFalse(report['passed'])
            np.testing.assert_array_equal(m.old,before)
            checkpoint=json.loads((Path(directory)/'initial.json').read_text())
            self.assertEqual(checkpoint['payload']['states'],before.tolist())

    def test_force_block_failure_is_visible_when_the_augmented_norm_hides_it(self):
        a=np.eye(2);a[1,1]=1000.;b=np.array([1e-10,1.]);x=np.linalg.solve(a,b);x[0]*=1.000000005
        result=linear_audit([dict(matrix=a.tolist(),rhs=b.tolist(),solution=x.tolist())])[0]
        self.assertLess(result['augmented'],1e-10)
        self.assertGreater(result['force_block'],1e-10)

    def test_nontransverse_degenerate_problem_is_not_a_simple_fold(self):
        with self.assertRaisesRegex(ValueError,'ambiguous'):
            fold.point([[0.,0.],[0.,1.]],[0.,1.],[1.,0.,0.])

    def test_new_checkpoint_binding_covers_control_scales_and_branch_payload(self):
        module=legacy();binding=control_binding(dict(test='new binding',step_size_m=.005,geometry_sha256='fixed geometry',runtime='fixed runtime'))
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'checkpoint.json'
            payload=dict(progress=.9,tangent=[1.,0.,0.],states=[1.,2.],displacement_m=[1.,2.],contacts=[],arc_length_m=.1)
            module.checkpoint(path,binding,payload)
            for field in ('fold_control','arc_lambda_scale_m','step_size_m','geometry_sha256','runtime'):
                changed=copy.deepcopy(binding);changed[field]='changed'
                with self.assertRaisesRegex(ValueError,'binding'):module.restore(path,changed)
            # A real accepted small-case checkpoint plus its original evidence.
            real=Path(directory)/'real';real.mkdir()
            module.follow(Analytic(),real,.005,3,30.,binding)
            import hashlib
            hashes={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in real.iterdir() if f.is_file()}
            (real/'evidence-sha256.json').write_text(json.dumps(hashes))
            from nc_fold_path import validate_checkpoint
            original=validate_checkpoint(real/'restart.json',binding)
            for field in ('progress','tangent'):
                changed=json.loads((real/'restart.json').read_text())
                changed['payload'][field]=.8 if field=='progress' else [-float(v) for v in changed['payload'][field]]
                changed['payload_sha256']=module.digest(changed['payload'])
                fake=real/(field+'-resealed.json');fake.write_text(json.dumps(changed))
                with self.assertRaisesRegex(ValueError,'original emitted'):validate_checkpoint(fake,binding)
            self.assertEqual(original,module.restore(real/'restart.json',binding))
            for field in ('progress','tangent','states'):
                changed=json.loads(path.read_text());changed['payload'][field]=0.
                other=Path(directory)/(field+'.json');other.write_text(json.dumps(changed))
                with self.assertRaisesRegex(ValueError,'payload differs'):module.restore(other,binding)

    def test_native_load_reconstruction_keeps_prescribed_components_and_does_not_edit_origin(self):
        import zipfile
        import probe_low_preload_nc as low
        from nc_fold_path import native_loads
        from verify_high_order import internal,prepare
        archive=ROOT/'interim/2026-11-08-low-preload-nc-followup/native-evidence.zip'
        with zipfile.ZipFile(archive) as stream:
            run=json.loads(stream.read('dense-1-128-min16384-v2/job/recovery-checkpoint-1.json'))['run']
        before=copy.deepcopy(run);mesh,_,_=low.expected_inputs(1/128,minimum_increment=1/16384,dense_lu=True)
        origin=run['continuation']['origin'];active=run['last_committed']['active_elements']
        final,holding=native_loads(mesh,run['ground'],run['case'],origin,active)
        np.testing.assert_allclose(holding,run['excavation_audits'][0]['initial_holding_load_kn'],rtol=1e-9,atol=1e-7)
        state=run['last_committed'];fold.audit_load(state['applied_nodal_load_kn'],final,holding,state['stage_progress'])
        q,offset=prepare(mesh);retained,_=internal(mesh,run['ground'],origin['material_states'],np.asarray(active),q,offset)
        with self.assertRaisesRegex(ValueError,'component load'):
            fold.audit_load(state['applied_nodal_load_kn'],final,retained.ravel()-final,state['stage_progress'])
        self.assertEqual(run,before)

    def test_original_proportional_series_changes_only_the_declared_initial_increment(self):
        from prepare_nc_fold_series import run
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)/'series';run(output)
            read=lambda n,k:json.loads((output/str(n)/(k+'.json')).read_text())
            self.assertEqual(read(128,'mesh'),read(256,'mesh'));self.assertEqual(read(128,'ground'),read(256,'ground'))
            a=read(128,'stages');b=read(256,'stages');self.assertEqual(a['stages'][2]['initial_increment'],1/128)
            self.assertEqual(b['stages'][2]['initial_increment'],1/256)
            b['stages'][2]['initial_increment']=1/128;self.assertEqual(a,b)
            self.assertEqual(a['settings']['relative_force_tolerance'],1e-8)
            self.assertEqual(a['settings']['relative_linear_tolerance'],1e-10)
            with self.assertRaises(FileExistsError):run(output)

    def test_fresh_native_adapter_refuses_changed_physics_and_unlisted_initial_steps(self):
        from nc_fold_fresh import FreshNative
        from probe_nc_midpoint import inputs
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            for index,change in enumerate(('physical','step')):
                folder=root/str(index);folder.mkdir();mesh,ground,case=inputs('midpoint')
                if change=='physical':ground['materials'][0]['mcc']['kappa']*=1.01
                else:case['stages'][2]['initial_increment']=1/64
                (folder/'mesh.json').write_text(json.dumps(mesh));(folder/'run.json').write_text(json.dumps(dict(ground=ground,case=case)))
                with self.assertRaisesRegex(ValueError,'fixed inputs|only original'):
                    FreshNative(folder,Path('must-not-launch'),folder/'material')


if __name__=='__main__':unittest.main()
