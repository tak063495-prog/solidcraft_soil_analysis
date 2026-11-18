"""Small input, order, preservation and seal contracts; no new FEM solve."""
import copy
import hashlib
import json
from pathlib import Path, PureWindowsPath
import subprocess
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
HERE=ROOT/'interim/2027-01-07-a03-explicit-shear-initial-stress'
sys.path.insert(0,str(HERE))
from a03shear_fixture import fixture, audit, plan, STRESSES
from a03shear_binding import dot, norm, common, fields
from a03shear_evidence import verify


def read(name):return json.loads((HERE/name).read_text(encoding='utf-8'))


class A03ExplicitShear(unittest.TestCase):
    def test_fixed_matrix_and_counts(self):
        p=plan();self.assertEqual(p,read('plan.json'))
        self.assertEqual((len(STRESSES),len(p['specs'])),(9,36))
        self.assertEqual(sum(len(s['engines']) for s in p['specs']),63)
        self.assertEqual(sum(4*s['n']*len(s['engines']) for s in p['specs']),2304)
        self.assertEqual(len({s['name'] for s in p['specs']}),36)
        self.assertEqual({s['height'] for s in p['specs']},{.5})

    def test_explicit_fresh_stress_and_unchanged_bodies(self):
        base=fixture(4,.5,'zero')
        for key,(index,tau) in STRESSES.items():
            m,g,c=fixture(4,.5,key)
            self.assertEqual(m,base[0]);self.assertEqual(g,base[1])
            stress=[2.5,2.5,5.,0.,0.,0.]
            if index is not None:stress[index]=tau
            self.assertEqual(c['stages'][0]['states']['UpperClay']['effective_stress_compression_positive_kpa'],stress)
            self.assertTrue(audit(m,g,c)['allowed'])

    def test_outside_matrix_refused_without_solve(self):
        for args in ((32,.5,'zero'),(8,.4,'xz5p'),(8,.5,'generated_k0')):
            with self.assertRaises(ValueError):fixture(*args)
        for mutation in ('both_shear','partial','darcy','mcc','old_history','support'):
            m,g,c=fixture(4,.5,'xz5p')
            if mutation=='both_shear':c['stages'][0]['states']['UpperClay']['effective_stress_compression_positive_kpa'][4]=.25
            if mutation=='partial':g['initialization']['water_level_m']=.75
            if mutation=='darcy':c['stages'].insert(0,dict(kind='steady_darcy',name='flow'))
            if mutation=='mcc':g['materials'][1]['modified_cam_clay']={}
            if mutation=='old_history':c['stages'][0]['states']['UpperClay']['plastic_work_kpa']=1.
            if mutation=='support':c['settings']['kinematics']='affine'
            with self.subTest(mutation=mutation),self.assertRaises((ValueError,AssertionError)):audit(m,g,c)

    def test_fixed_fsum_c_order_and_common_coordinates(self):
        self.assertEqual(dot([1e16,1,-1e16],[1,1,1]),1.)
        self.assertEqual(norm([[3,0],[0,4]]),5.)
        rows=[(fields.row(s,i/8),{'u':np.array([i/8])}) for s in ('geostatic','fill','contact','reload') for i in range(1,9)]
        self.assertEqual(len(common(rows,4)),16)
        self.assertEqual([r['progress'] for r,_ in common(rows,4)][:4],[.25,.5,.75,1.])

    def test_universal_newline_terminal_status(self):
        for name in ('README.ja.md','REPRODUCE.ja.md'):
            text=(HERE/name).read_text(encoding='utf-8')
            self.assertTrue(text.endswith('A-03は未達、passed=false。\n'))
            with tempfile.TemporaryDirectory() as d:
                p=Path(d)/name;p.write_bytes(text.replace('\n','\r\n').encode())
                self.assertEqual(p.read_text(encoding='utf-8'),text)

    def test_scope_null_and_failed_acceptance(self):
        a=read('acceptance.json');m=read('measurements.json')
        for k in ('passed','a03_complete','new_limited_acceptance','production_rust_changed','nc_code_or_configuration_changed','acceptance_tolerance_changed','old_values_rewritten','other_items_complete'):
            self.assertIs(a[k],False,k)
        self.assertTrue(a['refusals_retained'])
        self.assertTrue(all(v is None for v in m['water_accounting'].values()))
        self.assertEqual((m['formal_paths'],m['formal_states']),(63,2304))
        self.assertEqual(m['failures'],[])
        self.assertTrue(all(s['all_mc_elastic'] and s['all_contact_closed_stick'] for s in m['scope'].values()))
        self.assertEqual(sum(len(r) for r in m['raw'].values()),2304)

    def test_preserved_prior_source_sha(self):
        for path,sha in read('prior-record-binding.json').items():
            self.assertEqual(hashlib.sha256((ROOT/path).read_bytes()).hexdigest(),sha,path)

    def test_sealed_evidence_sha(self):
        self.assertTrue(verify(HERE)['all_sha256_verified'])

    def test_ci_sparse_posix_nul_paths_include_archive(self):
        parts=[HERE/p['name'] for p in read('archive-index.json')['parts']]
        prior=json.loads((ROOT/'interim/2026-12-30-a03-wet-soil-birth-cloud-followup/validation/sparse-checkout-contract.json').read_text())
        for group in ('python-core-0','python-core-1'):
            with tempfile.TemporaryDirectory() as d:
                rules=Path(d)/'rules';rules.write_text('\n'.join(prior['actual_core_sparse_rules'][group]['rules'])+'\n')
                paths=[p.relative_to(ROOT).as_posix() for p in parts]
                inputs=b''.join(p.encode('utf-8')+b'\0' for p in paths)
                r=subprocess.run(['git','sparse-checkout','check-rules','-z','--no-cone','--rules-file',str(rules)],cwd=ROOT,input=inputs,capture_output=True)
                self.assertEqual(r.returncode,0,r.stderr);self.assertEqual(r.stdout,inputs)
        windows=[PureWindowsPath('interim')/HERE.name/p['name'] for p in read('archive-index.json')['parts']]
        self.assertEqual([p.as_posix() for p in windows],paths)
