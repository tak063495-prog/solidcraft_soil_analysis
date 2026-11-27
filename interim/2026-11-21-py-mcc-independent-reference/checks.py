"""Explicit tests; outside the unchanged tests/tools discovery and CI counts."""
import ast
import copy
import csv
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np
from mcc_reference import (Material, initial, elastic, yield_value, integrate,
                           isotropic, richardson, mc_column)
from compare_history import validate, compare, load, read_json
from run_reference import history
from assess_series import assess
from column_reference import load_point, point_history, integrated_stages

HERE=Path(__file__).resolve().parent


def fixture():
    rows,_=integrate('k0',64,end=.01,samples=2)
    return history(rows,'k0',1.,Material(),'a'*64)


def series(axis='spatial'):
    common=dict(physics_sha256='a'*64,observable='uz',unit='m',sample_ids=['a','b'],
                depth_m=1.,boundary_id='same',load_path_id='same',load_increment=1/32,
                material_increment=1e-5,mesh_sha256='b'*64)
    ls=[]
    varying={'spatial':'h_m','load_increment':'load_increment','material_increment':'material_increment'}[axis]
    for h in [.25,.125,.0625]:
        ls.append(dict(common,**{varying:h},values=[3+h*h,7+2*h*h]))
    return dict(axis=axis,parameter_unit='m' if axis=='spatial' else '1',run_status='completed',nested_h_refinement=True,levels=ls)


class AnalyticPhysics(unittest.TestCase):
    def test_documented_initial_nc_and_oc_void_ratios(self):
        m=Material()
        self.assertAlmostEqual(initial(m)[3]-1,1.5789659628023816)
        self.assertAlmostEqual(initial(m,ocr=4)[3]-1,1.3571588650231994)

    def test_isotropic_compression_lambda_and_elastic_unloading_kappa(self):
        a,b,c=isotropic([100.,130.,100.])
        self.assertAlmostEqual((a['e']-b['e'])/math.log(1.3),.2)
        self.assertAlmostEqual((c['e']-b['e'])/math.log(1.3),.04)
        self.assertEqual(b['pc'],c['pc'])
        self.assertGreater(c['eps_v'],0.)

    def test_isotropic_overconsolidated_yield_split_has_exact_kappa_lambda(self):
        r=isotropic([100.,200.,800.],pc0=400.)
        self.assertAlmostEqual(r[1]['e']-r[2]['e'],.04*math.log(2)+.2*math.log(2))
        self.assertEqual(r[2]['pc'],800.)

    def test_elastic_undrained_shear_modulus(self):
        m=Material();y=initial(m,ocr=4.)
        out=elastic(.001,y,'undrained',m)
        self.assertEqual(out[0],y[0]);self.assertEqual(out[3],y[3])
        self.assertAlmostEqual(out[1]/.001,3*m.g_over_k*y[0]*y[3]/m.kap)

    def test_elastic_k0_ratio_is_poisson_lateral_ratio(self):
        m=Material();y=initial(m,ocr=4.)
        out=elastic(.001,y,'k0',m)
        dp,dq=out[:2]-y[:2]
        self.assertAlmostEqual((dp-dq/3)/(dp+2*dq/3),m.nu/(1-m.nu))
        self.assertAlmostEqual(out[3],y[3]*math.exp(-.001))

    def test_elastic_drained_radial_stress_constant(self):
        m=Material();y=initial(m,ocr=4.)
        for ea in [.001,.003,.01]:
            a=elastic(ea,y,'drained',m)
            self.assertAlmostEqual(a[0]-a[1]/3,100.)

    def test_exact_elastic_plastic_event_on_surface(self):
        m=Material()
        for mode in ('k0','drained','undrained'):
            rows,event=integrate(mode,64,ocr=4.)
            self.assertGreater(event,0.)
            self.assertLess(abs(yield_value(elastic(event,initial(m,ocr=4.),mode,m),m)),1e-13)
            self.assertEqual(rows[0]['pc'],math.exp(math.log(400.)))

    def test_k0_zero_lateral_strain_and_exact_volume(self):
        rows,_=integrate('k0',128)
        for r in rows:
            self.assertAlmostEqual(r['eps_v'],r['axial_strain'],places=12)
            self.assertAlmostEqual(1+r['e'],(1+rows[0]['e'])*math.exp(-r['eps_v']),places=9)

    def test_undrained_total_volume_constant(self):
        rows,_=integrate('undrained',128)
        self.assertTrue(all(r['eps_v']==0 and r['e']==rows[0]['e'] for r in rows))
        self.assertGreater(rows[-1]['q'],0)

    def test_drained_radial_effective_stress_constant(self):
        rows,_=integrate('drained',128)
        self.assertLess(max(abs(r['p']-r['q']/3-100) for r in rows),1e-10)

    def test_midpoint_order_from_independent_step_differences_all_six_paths(self):
        for mode in ('k0','drained','undrained'):
            for ocr in (1.,4.):
                with self.subTest(mode=mode,ocr=ocr):
                    rs=[integrate(mode,n,ocr=ocr)[0] for n in (32,64,128)]
                    for field in ('p','q','pc'):
                        out=richardson([[r[field] for r in rows] for rows in rs])
                        self.assertGreater(out['observed_order'],1.8)
                        self.assertLess(out['observed_order'],2.2)
                        self.assertIsNone(out['accepted'])

    def test_consistency_and_ncl_drift_reduce_without_projection(self):
        m=Material()
        for mode in ('k0','drained','undrained'):
            errors=[]
            for n in (32,64,128):
                rows,_=integrate(mode,n)
                errors.append(max(abs((r['q']/r['pc'])**2+m.M**2*r['p']/r['pc']*(r['p']/r['pc']-1)) for r in rows))
            self.assertLess(errors[1],.3*errors[0]);self.assertLess(errors[2],.3*errors[1])

    def test_physical_oc_dilation_not_clipped_to_monotonic_pc(self):
        rows,_=integrate('undrained',128,ocr=4.)
        self.assertLess(rows[-1]['pc'],rows[0]['pc'])
        with self.assertRaisesRegex(ValueError,'pc regression'):
            validate(history(rows,'undrained',4.,Material(),'a'*64))

    def test_invalid_material_and_nonpositive_p_rejected(self):
        for args in (dict(kap=.3),dict(nu=.5),dict(M=0),dict(lam=float('nan'))):
            with self.assertRaises(ValueError):Material(**args)
        with self.assertRaises(ValueError):isotropic([0.])
        with self.assertRaises(ValueError):initial(p=0.)

    def test_original_mc_column_three_printed_displacements(self):
        for stage,expected in [('gravity',-.0003333333333333333),('load',-.002700167229974),('unload',-.000700167229974)]:
            self.assertAlmostEqual(mc_column([.2],stage)['uz'][0],expected,places=14)

    def test_original_mc_column_equilibrium_and_unloading_increment(self):
        z=np.linspace(0,.2,17)
        load=mc_column(z,'load');unload=mc_column(z,'unload')
        self.assertTrue(np.allclose(np.diff(load['sigma_v'])/np.diff(z),-20))
        self.assertTrue(np.allclose(np.array(unload['uz'])-load['uz'],12*z/1200))
        self.assertAlmostEqual(load['uz'][0],0.)

    def test_spatial_p1_quadratic_interpolation_error_has_order_two(self):
        z=np.linspace(0,.2,4097);exact=np.array(mc_column(z,'gravity')['uz'])
        errors=[]
        for n in (4,8,16):
            nodes=np.linspace(0,.2,n+1)
            values=np.interp(z,nodes,mc_column(nodes,'gravity')['uz'])
            errors.append(np.linalg.norm(values-exact))
        for a,b in zip(errors,errors[1:]):self.assertAlmostEqual(math.log(a/b,2),2.,places=7)

    def test_reference_imports_never_call_existing_solvers(self):
        for name in ('mcc_reference.py','run_reference.py','compare_history.py','assess_series.py','column_reference.py'):
            tree=ast.parse((HERE/name).read_text())
            imports=[]
            for node in ast.walk(tree):
                if isinstance(node,ast.Import):imports.extend(x.name.split('.')[0] for x in node.names)
                elif isinstance(node,ast.ImportFrom):imports.append(node.module.split('.')[0])
            self.assertTrue(set(imports)<={'dataclasses','math','numpy','scipy','argparse','csv','hashlib','json','pathlib','platform','time','mcc_reference','compare_history','column_reference'},imports)


class DocumentedMccColumn(unittest.TestCase):
    def test_vertical_equilibrium_load_unload_cut_and_removed_point_freeze(self):
        for z in (.125,.375,.625,.875):
            hs=point_history(z,128)
            expected=[100+20*(1-z),130+20*(1-z),100+20*(1-z),100+20*(.5-z) if z<=.5 else 100+20*(1-z)]
            for y,sv in zip(hs,expected):
                self.assertAlmostEqual(y[0]+2*y[1]/3,sv,places=9)
            self.assertEqual(hs[1][2],hs[2][2]);self.assertEqual(hs[2][2],hs[3][2])
            if z>.5:self.assertTrue(np.array_equal(hs[2],hs[3]))

    def test_column_material_step_order_is_two_at_fixed_physical_point(self):
        out=richardson([load_point(.25,n) for n in (64,128,256)])
        self.assertGreater(out['observed_order'],1.9)
        self.assertLess(out['observed_order'],2.1)

    def test_column_spatial_quadrature_order_at_fixed_material_steps(self):
        out=richardson([integrated_stages(.5,n=128,cells=c) for c in (4,8,16)])
        self.assertGreater(out['observed_order'],1.9)
        self.assertLess(out['observed_order'],2.1)

    def test_unloading_has_rebound_but_retains_plastic_settlement(self):
        load,unload,cut=integrated_stages(.5,n=128,quadrature=8)
        self.assertLess(load,unload);self.assertLess(unload,0.)
        self.assertGreater(cut,unload)


class HistoryContract(unittest.TestCase):
    def setUp(self):self.data=fixture()
    def reject(self,change,pattern=None):
        a=copy.deepcopy(self.data);change(a)
        with self.assertRaisesRegex((ValueError,KeyError,TypeError),pattern or '.*'):validate(a)

    def test_valid_full_history_comparison_is_diagnostic_only(self):
        r=compare(self.data,self.data)
        self.assertIsNone(r['accepted'])
        self.assertTrue(all(v['max_abs']==0 and v['accepted'] is None for v in r['metrics'].values()))

    def test_units_and_stress_sign_and_strain_measure_are_not_guessed(self):
        for key,value in [('stress_unit','Pa'),('length_unit','mm'),('stress_sign','tension_positive'),('strain_sign','tension_positive'),('strain_measure','engineering'),('stress_kind','total')]:
            self.reject(lambda a:a['metadata'].update({key:value}),'metadata rejected')

    def test_nonpositive_effective_pressure_rejected(self):
        for value in (0,-1):self.reject(lambda a:a['rows'][1].update(p=value),'p<=0')

    def test_pc_regression_rejected(self):
        self.reject(lambda a:a['rows'][1].update(pc=99.),'pc regression')

    def test_nonfinite_missing_duplicate_and_invented_history_rejected(self):
        self.reject(lambda a:a['rows'][1].update(p=float('nan')),'nonfinite')
        self.reject(lambda a:a['rows'].pop(),'missing history')
        self.reject(lambda a:a['rows'].append(a['rows'][0]),'duplicate')
        self.reject(lambda a:a['rows'][1].pop('e'),'history fields')
        self.reject(lambda a:a['metadata'].update(history_completion='forward_fill'),'metadata rejected')

    def test_e_and_volumetric_strain_inconsistent_rejected(self):
        self.reject(lambda a:a['rows'][1].update(eps_v=.1),'inconsistent')

    def test_timeout_cancelled_partial_not_completed(self):
        for value in ('timeout','cancelled','partial'):
            self.reject(lambda a:a['metadata'].update(run_status=value),'metadata rejected')

    def test_missing_initial_state_load_grid_and_wrong_step_rejected(self):
        self.reject(lambda a:a['rows'].pop(0),'missing history')
        self.reject(lambda a:a['rows'][1].update(load_coordinate=999),'coordinate mismatch')
        self.reject(lambda a:a['rows'][1].update(step=.5),'invalid step')
        self.reject(lambda a:a['rows'][1].update(p=True),'boolean')

    def test_different_physical_problem_or_sampling_refused(self):
        for key,value in [('physics_sha256','b'*64),('case_id','different'),('ids',['other'])]:
            b=copy.deepcopy(self.data);b['metadata'][key]=value
            with self.assertRaises(ValueError):compare(self.data,b)

    def test_node_history_metrics_and_zero_reference_relative_error_null(self):
        a=copy.deepcopy(self.data);a['metadata']['kind']='node'
        a['rows']=[dict(step=r['step'],id=r['id'],load_coordinate=r['load_coordinate'],ux=0.,uy=0.,uz=-r['step']*.001) for r in a['rows']]
        b=copy.deepcopy(a);b['rows'][-1]['uz']-=.0001
        r=compare(a,b)
        self.assertAlmostEqual(r['metrics']['uz']['max_abs'],.0001)
        self.assertIsNone(r['metrics']['ux']['relative_l2'])

    def test_json_csv_equivalence_and_no_implicit_csv_metadata(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'a.json').write_text(json.dumps(self.data));(p/'m.json').write_text(json.dumps(self.data['metadata']))
            with (p/'a.csv').open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=list(self.data['rows'][0]));w.writeheader();w.writerows(self.data['rows'])
            self.assertEqual(compare(load(p/'a.json'),load(p/'a.csv',p/'m.json'))['metrics']['p']['max_abs'],0.)
            with self.assertRaises(ValueError):load(p/'a.csv')

    def test_duplicate_json_keys_and_csv_headers_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'bad.json').write_text('{"a":1,"a":2}')
            with self.assertRaisesRegex(ValueError,'duplicate JSON'):read_json(p/'bad.json')
            (p/'bad.csv').write_text('step,step\n0,0\n');(p/'meta.json').write_text(json.dumps(self.data['metadata']))
            with self.assertRaisesRegex(ValueError,'duplicate CSV'):load(p/'bad.csv',p/'meta.json')

    def test_cli_rejects_before_output_and_never_overwrites(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);source=p/'a.json';source.write_text(json.dumps(self.data));out=p/'result.json'
            cmd=[sys.executable,str(HERE/'compare_history.py'),'--reference',str(source),'--candidate',str(source),'--output',str(out)]
            r=subprocess.run(cmd,capture_output=True,text=True)
            self.assertEqual(r.returncode,0,r.stderr)
            self.assertIsNone(read_json(out)['accepted'])
            before=out.read_bytes();r=subprocess.run(cmd,capture_output=True)
            self.assertNotEqual(r.returncode,0);self.assertEqual(out.read_bytes(),before)
            source.write_text('{}');other=p/'reject.json';cmd[-1]=str(other)
            r=subprocess.run(cmd,capture_output=True);self.assertEqual(r.returncode,2);self.assertFalse(other.exists())


class RefinementContract(unittest.TestCase):
    def test_known_limit_is_recovered_from_mesh_differences(self):
        out=assess(series())
        self.assertAlmostEqual(out['observed_order'],2.)
        self.assertTrue(np.allclose(out['extrapolated'],[3,7]))
        self.assertIsNone(out['accepted'])

    def test_material_and_global_load_series_are_separate(self):
        for axis in ('material_increment','load_increment'):
            out=assess(series(axis));self.assertAlmostEqual(out['observed_order'],2.)
            a=series(axis);a['levels'][1]['mesh_sha256']='other'
            with self.assertRaisesRegex(ValueError,'mesh changed'):assess(a)

    def test_depth_material_or_load_change_confounds_spatial_series(self):
        for key,value in [('depth_m',2.),('material_increment',1e-6),('load_increment',1/64),('physics_sha256','other'),('boundary_id','other')]:
            a=series();a['levels'][1][key]=value
            with self.assertRaisesRegex(ValueError,'confounded'):assess(a)

    def test_missing_shared_samples_non_h_series_and_incomplete_refused(self):
        a=series();a['levels'][1]['values'].pop()
        with self.assertRaisesRegex(ValueError,'missing common'):assess(a)
        a=series();a['nested_h_refinement']=False
        with self.assertRaisesRegex(ValueError,'nested'):assess(a)
        a=series();a['run_status']='timeout'
        with self.assertRaisesRegex(ValueError,'incomplete'):assess(a)

    def test_nonconvergent_roundoff_and_direction_reversal_do_not_get_orders(self):
        for values in ([1,2,1.9],[1,2,4],[1,1,1],[1,1+1e-15,1+1.25e-15]):
            out=richardson(values)
            self.assertIsNone(out['observed_order']);self.assertIsNone(out['extrapolated'])
        with self.assertRaises(ValueError):richardson([1,float('nan'),2])


if __name__=='__main__':
    unittest.main(verbosity=2)
