"""Explicit standalone unittest entry; intentionally outside test_*.py discovery."""
import copy
import csv
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from reference import Model,material,contact,EVENTS,write
from compare import packet,validate,compare,read

class Constitutive(unittest.TestCase):
    def test_vg_limits_and_saturation_are_distinct(self):
        h=np.array([-100.,-.5,0.,1.]);w,c,k,dk=material(h,np.zeros(4))
        self.assertTrue(np.all(np.diff(w)>0));self.assertEqual(k[-1],.02)
        self.assertEqual(w[-1],.4117647058823529+.02)
        self.assertTrue(np.all(c>0));self.assertTrue(np.all(dk>=0))
    def test_analytic_storage_and_conductivity_derivatives(self):
        h=np.array([-2.,-.5,-.02,.2]);delta=1e-7
        w,c,k,dk=material(h,np.zeros(4));plus=material(h+delta,np.zeros(4));minus=material(h-delta,np.zeros(4))
        np.testing.assert_allclose(c,(plus[0]-minus[0])/(2*delta),rtol=2e-7,atol=1e-10)
        np.testing.assert_allclose(dk,(plus[2]-minus[2])/(2*delta),rtol=2e-7,atol=1e-10)
    def test_linear_saturated_limit_matches_vg(self):
        for a,b in zip(material(np.array([2.,3.]),np.array([0.,1.]),'linear'),material(np.array([2.,3.]),np.array([0.,1.]),'vg'),strict=True):
            np.testing.assert_array_equal(a,b)
    def test_normal_complementarity_and_external_spring_equilibrium(self):
        for dry in (True,False):
            c=contact(np.array([.2,.4,.8]),np.full(3,.25),.25,np.zeros(3),dry)
            np.testing.assert_allclose(c['effective']-2000*c['gap'],c['load']-c['pressure'],atol=1e-14)
            np.testing.assert_array_equal(c['effective']*c['gap'],np.zeros(3))
    def test_coulomb_bound_and_open_traction(self):
        old=np.zeros(3)
        for t in np.linspace(0,1,81):
            c=contact(np.array([.1,.4,.8]),np.full(3,.25),t,old)
            self.assertTrue(np.all(np.abs(c['shear'])<=.3*c['effective']))
            self.assertTrue(np.all(c['shear'][c['gap']>0]==0));old=c['plastic']

class FV(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r=Model(6,1).run(.05)
    def test_vg_one_dimensional_path_contains_unsaturated_cells(self):
        self.assertLess(min(h-.32 for frame in self.r['hydraulic_frames'] for h in frame['head']),0.)
    def test_small_initial_saturated_flux_preserves_existing_ledger_bound(self):
        r=Model(16,1,'linear').run(.003125,.00625)
        validate(packet(r))
    def test_closed_interface_has_common_water_head(self):
        for row in self.r['rows']:
            if not row['hydraulic_open']:self.assertAlmostEqual(row['head_a_m'],row['head_b_m'],places=11)
    def test_gap_has_real_storage_and_distinct_in_out_flows(self):
        self.assertGreater(max(abs(r['q_a_to_joint_m3_s']-r['q_joint_to_b_m3_s']) for r in self.r['rows']),1e-7)
    def test_cumulative_water_conservation_is_not_only_endpoint(self):
        self.assertLess(self.r['maximum_water_balance_defect_m3'],1e-10)
        self.assertLess(abs(self.r['final_storage_m3']-self.r['initial_storage_m3']-self.r['cumulative_net_inflow_m3']),1e-10)
    def test_reservoir_exchange_is_included_in_budget(self):
        r=Model(6,1,policy='PRESCRIBED').run(.05)
        self.assertLess(r['maximum_water_balance_defect_m3'],1e-10)
        for row in r['rows']:
            if row['hydraulic_open']:self.assertAlmostEqual(row['joint_head_m'],.3+.07*np.sin(np.pi*row['time_s'])**2,places=10)
    def test_one_way_hydraulics_are_independent_of_mechanics(self):
        dry=Model(6,1,dry=True).run(.05)
        np.testing.assert_array_equal(self.r['final_head'],dry['final_head'])
        self.assertNotEqual([r['gap_m'] for r in self.r['rows']],[r['gap_m'] for r in dry['rows']])
    def test_two_way_opening_changes_hydraulic_solution(self):
        wet=Model(6,1,mode='TWO_WAY').run(.05);dry=Model(6,1,mode='TWO_WAY',dry=True).run(.05)
        self.assertGreater(np.linalg.norm(np.array(wet['final_head'])-dry['final_head']),1e-5)
        self.assertLess(wet['maximum_water_balance_defect_m3'],1e-10)
    def test_two_dimensional_local_flux_equations_and_jacobian(self):
        m=Model(6,3,mode='TWO_WAY');old=np.full(m.size,.3);x=np.linspace(.325,.38,m.size)
        r,J,_=m.residual(x,old,0.,.2);direction=np.sin(np.arange(m.size)+1)
        eps=1e-7
        numerical=(m.residual(x+eps*direction,old,0.,.2)[0]-m.residual(x-eps*direction,old,0.,.2)[0])/(2*eps)
        np.testing.assert_allclose(J@direction,numerical,rtol=1e-6,atol=2e-9)
    def test_shifted_event_endpoints_have_exactly_zero_aperture(self):
        for shift in (-.0125,0.,.0125):
            m=Model(switch_shift=shift)
            for event in EVENTS:
                self.assertEqual(m.aperture(np.array([.35]),round(event+shift,12))[0][0],0.)
    def test_physical_events_and_fixed_generation_are_retained(self):
        ts=sorted(set(r['time_s'] for r in self.r['rows']))
        for event in EVENTS:self.assertIn(event,ts)
        self.assertEqual({r['generation'] for r in self.r['rows']},{0})
    def test_no_false_generation_or_unsupported_soil_connection(self):
        with self.assertRaisesRegex(ValueError,'unconfirmed'):Model(policy='SOIL')
        with self.assertRaisesRegex(ValueError,'not defined'):Model(mode='TWO_WAY',policy='PRESCRIBED')
    def test_saturated_interface_eigenmode_converges(self):
        errors=[]
        for n in (4,8,16):
            r=Model(n,1,'linear',schedule=False,eigen=True).run(.0005,.02)
            exact=2+.1*np.sin(np.pi*(np.arange(n)+.5)/n)*np.exp(-np.pi**2*.02)
            errors.append(np.linalg.norm(np.array(r['final_head'][:n])-exact)/np.sqrt(n))
        self.assertGreater(errors[0]/errors[1],2);self.assertGreater(errors[1]/errors[2],2)

class Contract(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.p=packet(Model(4,1).run(.1))
    def test_roundtrip_never_accepts_a04_c04(self):
        r=compare(self.p,self.p);self.assertIsNone(r['accepted']);self.assertEqual(r['tolerance'],'未定義')
        self.assertFalse(r['actual_body_comparison'])
    def test_two_way_label_without_aperture_feedback_is_rejected(self):
        p=packet(Model(4,1,mode='TWO_WAY').run(.1))
        p['rows'][0]['hydraulic_aperture_m']+=.001
        p['rows'][0]['hydraulic_open']=True
        with self.assertRaisesRegex(ValueError,'TWO_WAY aperture'):validate(p)
    def test_numeric_difference_remains_visible_with_undefined_tolerance(self):
        p=copy.deepcopy(self.p);p['rows'][0]['head_a_m']+=.1
        r=compare(self.p,p);self.assertIsNone(r['accepted']);self.assertAlmostEqual(r['metrics']['head_a_m']['maximum_absolute'],.1)
    def test_csv_and_json_carry_identical_explicit_contract(self):
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'candidate.csv';write(path.with_suffix('.metadata.json'),self.p['metadata'])
            with path.open('w',newline='') as f:
                writer=csv.DictWriter(f,fieldnames=list(self.p['rows'][0]));writer.writeheader()
                for row in self.p['rows']:writer.writerow(dict(row,hydraulic_open=str(row['hydraulic_open']).lower()))
            self.assertEqual(read(path),self.p)


def mutation_test(change):
    def test(self):
        p=copy.deepcopy(self.p);change(p)
        with self.assertRaises((ValueError,KeyError,TypeError)):compare(self.p,p)
    return test

MUTATIONS={
 'reversed_numeric_flow':lambda p:p['rows'][0].update(q_a_to_joint_m3_s=-p['rows'][0]['q_a_to_joint_m3_s']),
 'missing_reservoir_source':lambda p:p['rows'][0].pop('reservoir_inflow_m3_s'),
 'wrong_units':lambda p:p['metadata']['units'].update(head='mm'),
 'wrong_flow_sign':lambda p:p['metadata']['signs'].update(flow_b='outward_from_b'),
 'pressure_head_instead_of_total':lambda p:p['metadata']['signs'].update(head='pressure_head'),
 'load_progress_instead_of_seconds':lambda p:p['metadata'].update(time_basis='load_progress'),
 'implicit_interpolation':lambda p:p['metadata'].update(interpolation='linear'),
 'missing_point':lambda p:p['rows'].pop(),
 'duplicate_point':lambda p:p['rows'].append(p['rows'][-1]),
 'reordered_rows':lambda p:p['rows'].reverse(),
 'changed_generation':lambda p:p['rows'][0].update(generation=1),
 'missing_generation':lambda p:p['rows'][0].pop('generation'),
 'wrong_one_way_mode':lambda p:p['metadata']['configuration'].update(mode='TWO_WAY'),
 'unfinished_output':lambda p:p['metadata'].update(complete=False),
 'nonfinite_head':lambda p:p['rows'][0].update(head_a_m=float('nan')),
 'negative_effective_pressure':lambda p:p['rows'][0].update(effective_pressure_kpa=-1),
 'open_flag_as_integer':lambda p:p['rows'][0].update(hydraulic_open=0),
 'missing_joint_id':lambda p:p['rows'][0].pop('joint_id'),
 'wrong_contact_state':lambda p:p['rows'][0].update(contact_state='unknown'),
 'changed_time_sample':lambda p:p['rows'][0].update(time_s=.10001),
 'unknown_soil_mapping':lambda p:p['metadata']['configuration'].update(policy='SOIL'),
 'averaged_input_parameters':lambda p:p['metadata']['parameters'].update(ks_m_s=.01),
}
for name,change in MUTATIONS.items():setattr(Contract,'test_reject_'+name,mutation_test(change))

if __name__=='__main__':unittest.main(verbosity=2)
