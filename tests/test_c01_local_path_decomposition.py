"""Small independent analytic checks; no archived global re-solves in CI."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
HERE=ROOT/'interim/2026-12-15-c01-local-path-decomposition'
sys.path.insert(0,str(HERE))
from c01_local_oracle import advance_exact,initial,refine
from c01_path_replay import integrate
from c01_local_study import make_request,compare
from c01_decomposition import dot,metric,pair_budget,STAGES,FIELDS


class LocalPathDecompositionTests(unittest.TestCase):
    def test_constant_pressure_closed_form(self):
        s=advance_exact(initial(),[0.,-.1],100.,100.,.5)
        s=advance_exact(s,[.2,-.1],100.,100.,.5)
        self.assertAlmostEqual(s['tau'],5.)
        self.assertAlmostEqual(s['slip'],.15)
        self.assertAlmostEqual(s['dissipation'],.75)

    def test_linear_pressure_exact_yield_root(self):
        s=advance_exact(initial(),[0.,-.1],100.,100.,.5)
        s=advance_exact(s,[.2,-.2],100.,100.,.5)
        self.assertAlmostEqual(s['tau'],10.)
        self.assertAlmostEqual(s['slip'],.1)
        self.assertAlmostEqual(s['dissipation'],5/6)
        q,dp,d,_,_=integrate(np.array([[0.,0.,-.1]]),np.array([[0.,0.,-10.]]),np.array([[.2,0.,-.2]]),np.array([[0.,0.,1.]]),np.array([100.]),np.array([100.]),np.array([.5]),1)
        self.assertAlmostEqual(float(d[0]),1.)
        self.assertGreater(float(d[0]),s['dissipation'])

    def test_pressure_release_at_fixed_displacement(self):
        s=dict(jump=[.1,-.1],tau=5.,slip=.05,accumulated=.05,dissipation=.25)
        end=advance_exact(s,[.1,-.02],100.,100.,.5)
        self.assertAlmostEqual(end['slip']-s['slip'],.04)
        self.assertAlmostEqual(end['dissipation']-s['dissipation'],.12)

    def test_reversal_opening_and_subdivision_invariant_exact_integral(self):
        vertices=[[0.,-.1],[.2,-.15],[-.25,-.05],[.3,.1],[-.1,.05],[-.4,-.1],[.4,-.2]]
        values=[]
        for count in (1,2,7,32):
            s=initial();previous=0.
            for _,_,target in refine(vertices,count):
                s=advance_exact(s,target,100.,100.,.5)
                self.assertGreaterEqual(s['dissipation'],previous-1e-14)
                previous=s['dissipation']
            values.append(s)
        for v in values[1:]:
            for key in ('tau','slip','accumulated','dissipation'):
                self.assertAlmostEqual(v[key],values[0][key],places=12)

    def test_signed_components_require_cross_terms(self):
        total=np.array([1.,2.]);a=np.array([3.,-1.]);b=total-a
        self.assertEqual(dot(total,total),dot(a,a)+dot(b,b)+2*dot(a,b))
        self.assertNotEqual(dot(total,total),dot(a,a)+dot(b,b))
        self.assertEqual(metric(total,np.ones(2)),metric(total.copy(order='F'),np.ones(2).copy(order='F')))

    def test_signed_path_budget_with_known_switch_increment(self):
        coords=[(stage,1.) for stage in STAGES]
        normal=np.array([[0.,0.,1.],[0.,0.,1.]])
        series=[]
        for fine in (False,True):
            fields={};refined={}
            for key in FIELDS:
                base=np.array([[i+1.,2*(i+1.)] for i in range(7)]) if fine else np.zeros((7,2))
                local=np.array([3.,-1.]) if fine else np.zeros(2)
                if key=='plastic_slip':
                    base=np.stack((base,np.zeros_like(base),np.zeros_like(base)),axis=-1)
                    local=np.stack((local,np.zeros_like(local),np.zeros_like(local)),axis=-1)
                refined[key]=base;fields[key]=base+local
            fields['u']=np.zeros((7,2,3));fields['jump']=np.zeros((7,2,3))
            contacts={(stage,gp):[] for stage in STAGES for gp in range(2)}
            if fine:contacts[('open',0)]=[(.25,.5)]
            series.append(dict(coords=coords,index={c:i for i,c in enumerate(coords)},fields=fields,refined=refined,contacts=contacts))
        summary,budget=pair_budget(*series,1,normal)
        np.testing.assert_array_equal(budget['components']['plastic_dissipation']['local_integration'][2],[3.,-1.])
        np.testing.assert_array_equal(budget['components']['plastic_dissipation']['displacement_path_switch'][2],[1.,0.])
        np.testing.assert_array_equal(budget['components']['plastic_dissipation']['displacement_path_same_contact'][2],[2.,6.])
        self.assertEqual(summary['fields']['plastic_dissipation']['maximum_identity_residual'],0.)

    def test_native_local_law_against_exact_paths_and_independent_endpoint(self):
        # Four analytical paths, a bounded representative outer series. The
        # complete 168-run study is archived, not repeated as a CI experiment.
        subprocess.run(['cargo','build','--locked','--release','--example','wet_joint_local_path_probe'],cwd=ROOT,check=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        exe=ROOT/'target/release/examples'/('wet_joint_local_path_probe.exe' if sys.platform=='win32' else 'wet_joint_local_path_probe')
        req,_=make_request()
        req['runs']=[r for r in req['runs'] if '/n1/' in r['id'] or '/n32/fixed-1' in r['id']]
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);(p/'request.json').write_text(json.dumps(req),encoding='utf-8')
            subprocess.run([str(exe),str(p/'request.json'),str(p/'response.json')],cwd=ROOT,check=True)
            summary=compare(json.loads((p/'response.json').read_text(encoding='utf-8')))
        for row in summary['runs'].values():
            self.assertLess(row['maximum_slip_error_m'],1e-14)
            self.assertLess(row['maximum_traction_error_kpa'],1e-10)
            self.assertGreaterEqual(row['minimum_increment_dissipation_kpa_m'],0.)
            self.assertLess(row['maximum_independent_endpoint_difference'],1e-10)
        for path in ('rising_pressure','reverse_variable_pressure','opening_recontact_reverse','pressure_release_at_fixed_slip'):
            coarse=summary['runs'][path+'/n1/fixed-1'];fine=summary['runs'][path+'/n32/fixed-1']
            self.assertLess(fine['maximum_cumulative_dissipation_error_kpa_m'],coarse['maximum_cumulative_dissipation_error_kpa_m'])
        # A positive, nonnegative dissipative algorithm need not converge monotonically.
        a=summary['runs']['opening_recontact_reverse/n1/fixed-2']['final_dissipation_error_kpa_m']
        b=summary['runs']['opening_recontact_reverse/n1/fixed-4']['final_dissipation_error_kpa_m']
        self.assertGreater(abs(b),abs(a))

    def test_new_sealed_evidence_integrity_and_original_failure(self):
        import c01_local_evidence
        c01_local_evidence.verify()
        acceptance=json.loads((HERE/'acceptance.json').read_text(encoding='utf-8'))
        self.assertIs(acceptance['passed'],False)
        self.assertEqual(acceptance['c01_decision'],'未達')

    def test_documents_survive_real_crlf(self):
        for name in ('README.ja.md','REPRODUCE.ja.md'):
            text=(HERE/name).read_text(encoding='utf-8')
            self.assertTrue(text.endswith('C-01は未達、passed=false。\n'))
            with tempfile.TemporaryDirectory() as tmp:
                path=Path(tmp)/name;path.write_bytes(text.replace('\n','\r\n').encode('utf-8'))
                self.assertEqual(path.read_text(encoding='utf-8'),text)

if __name__=='__main__':unittest.main()
