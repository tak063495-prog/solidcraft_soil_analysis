"""Small input/order/archive contracts; no FEM solve or new Rust integration work."""
import hashlib, json
from pathlib import Path, PureWindowsPath
import subprocess, sys, tempfile, unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
HERE=ROOT/'interim/2027-01-06-a03-unequal-thickness-anisotropic-k'
sys.path.insert(0,str(HERE))
from a03aniso_fixture import fixture, audit, plan, HEIGHTS, STRESSES
from a03aniso_binding import dot, norm, common, fields
from a03aniso_evidence import verify

def read(name):return json.loads((HERE/name).read_text(encoding='utf-8'))

class A03UnequalAnisotropic(unittest.TestCase):
    def test_fixed_matrix(self):
        p=plan();self.assertEqual(p,read('plan.json'))
        self.assertEqual(sum(len(s['engines']) for s in p['specs']),63)
        self.assertEqual(sum(4*s['n']*len(s['engines']) for s in p['specs']),2304)
        self.assertEqual(len({s['name'] for s in p['specs']}),36)
        self.assertTrue(all(s['minimum']==1 and not any(s[k] for k in ('element','face','reverse')) for s in p['specs']))

    def test_unequal_layer_bounds_volume_and_midpoints(self):
        for h in HEIGHTS:
            m,g,c=fixture(4,h,'k58');a=audit(m,g,c)
            self.assertEqual(a['bounds_m'],[[0.,0.,0.],[1.,1.,1.]])
            self.assertEqual(a['body_volumes_m3'],{'LowerClay':h,'UpperClay':1-h})
            self.assertEqual(a['total_volume_m3'],1.)
            self.assertLess(a['midpoint_maximum_defect_m'],1e-15)
            lower=next(r for r in m['regions'] if r['body']=='LowerClay')
            self.assertEqual(lower['bounds_m'][1][2],h)

    def test_explicit_directional_constant_fresh_stress(self):
        for s,(kx,ky) in STRESSES.items():
            m,g,c=fixture(8,.4,s)
            for name,z in (('LowerClay',90.),('UpperClay',70.)):
                self.assertEqual(g['initialization']['geostatic']['seed']['states'][name]['effective_stress_compression_positive_kpa'],[kx*z,ky*z,z,0.,0.,0.])
            f=c['stages'][0]['states']['UpperClay']
            self.assertEqual(f['effective_stress_compression_positive_kpa'],[kx*5,ky*5,5.,0.,0.,0.])
            self.assertEqual(f['plastic_strain_compression_positive'],[0.]*6)
            self.assertEqual((f['equivalent_plastic_strain'],f['plastic_work_kpa'],f['last_return']),(0.,0.,'elastic'))

    def test_outside_matrix_is_rejected_before_solve(self):
        for args in ((32,.4,'k58'),(8,.3,'k58'),(8,.4,'generated_k0')):
            with self.assertRaises(ValueError):fixture(*args)
        m,g,c=fixture(4,.4,'k58');g['initialization']['water_level_m']=.75
        with self.assertRaises(AssertionError):audit(m,g,c)

    def test_fixed_fsum_and_common_coordinates(self):
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

    def test_scope_failure_and_null_are_not_acceptance(self):
        a=read('acceptance.json');m=read('measurements.json')
        for k in ('passed','a03_complete','new_limited_acceptance','production_rust_changed','nc_code_or_configuration_changed','acceptance_tolerance_changed','old_values_rewritten','other_items_complete'):
            self.assertIs(a[k],False,k)
        self.assertTrue(a['refusals_retained'])
        self.assertTrue(all(v is None for v in m['water_accounting'].values()))
        self.assertEqual((m['formal_paths'],m['formal_states']),(63,2304))
        self.assertEqual(m['failures'],[])
        self.assertTrue(all(x['u_stress_traction_below_original_1e7'] and x['force_absolute_below_original_1e10'] and x['categorical_mismatches']==0 for x in m['implementation_gates'].values()))

    def test_sealed_raw_evidence_sha(self):
        self.assertTrue(verify(HERE)['all_sha256_verified'])

    def test_prior_records_sha_unchanged(self):
        for p,sha in read('prior-record-binding.json').items():
            self.assertEqual(hashlib.sha256((ROOT/p).read_bytes()).hexdigest(),sha,p)

    def test_ci_sparse_selects_new_archive(self):
        parts=[HERE/p['name'] for p in read('archive-index.json')['parts']]
        prior=json.loads((ROOT/'interim/2026-12-30-a03-wet-soil-birth-cloud-followup/validation/sparse-checkout-contract.json').read_text())
        for group in ('python-core-0','python-core-1'):
            with tempfile.TemporaryDirectory() as d:
                rules=Path(d)/'rules'
                rules.write_text('\n'.join(prior['actual_core_sparse_rules'][group]['rules'])+'\n')
                paths=[p.relative_to(ROOT).as_posix() for p in parts]
                inputs=b''.join(p.encode('utf-8')+b'\0' for p in paths)
                r=subprocess.run(['git','sparse-checkout','check-rules','-z','--no-cone','--rules-file',str(rules)],cwd=ROOT,input=inputs,capture_output=True)
                self.assertEqual(r.returncode,0,r.stderr)
                self.assertEqual(r.stdout,inputs)

    def test_sparse_nul_protocol_preserves_special_paths_and_exclusions(self):
        # These are logical Git paths, not files created on either host OS.
        paths=['kept/space name.zip', 'kept/quote"name.zip', 'kept/back\\slash.zip',
               'kept/tab\tname.zip', 'kept/line\nname.zip', 'kept/日本語.zip',
               'kept/duplicate.zip', 'kept/duplicate.zip', 'excluded/absent.zip']
        inputs=b''.join(p.encode('utf-8')+b'\0' for p in paths)
        expected=b''.join(p.encode('utf-8')+b'\0' for p in paths[:-1])
        with tempfile.TemporaryDirectory() as d:
            rules=Path(d)/'rules';rules.write_bytes(b'/kept/\n')
            for quote in ('true','false'):
                with self.subTest(core_quotePath=quote):
                    r=subprocess.run(['git','-c','core.quotePath='+quote,'sparse-checkout','check-rules','-z','--no-cone','--rules-file',str(rules)],cwd=ROOT,input=inputs,capture_output=True)
                    self.assertEqual(r.returncode,0,r.stderr)
                    self.assertEqual(r.stdout,expected)

    def test_sparse_windows_repo_paths_are_posix_before_git(self):
        paths=[PureWindowsPath('interim')/HERE.name/p['name'] for p in read('archive-index.json')['parts']]
        expected=[(HERE/p['name']).relative_to(ROOT).as_posix() for p in read('archive-index.json')['parts']]
        self.assertEqual([p.as_posix() for p in paths],expected)
        self.assertTrue(all('\\' in str(p) and '\\' not in p.as_posix() for p in paths))
