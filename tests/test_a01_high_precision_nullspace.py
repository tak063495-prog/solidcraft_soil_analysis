"""Light sealed-evidence contracts; no mpmath or solver invocation in CI."""
from decimal import Decimal
import hashlib
from pathlib import Path
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1]
HERE=ROOT/'interim/2027-01-04-a01-high-precision-nullspace'
sys.path.insert(0,str(HERE))
from a01hp_common import read
from a01hp_evidence import verify
from a01hp_report import collect
INDEX_SHA='5d091e84503335c7bb795f7f5011cb55c3b88a242c75437c5523df02126ab1c0'


class A01HighPrecisionNullspace(unittest.TestCase):
    def test_scope_and_failed_global_acceptance_retained(self):
        a=read(HERE/'acceptance.json');m=read(HERE/'measurements.json')
        for value in (a,m):
            self.assertIs(value['passed'],False)
            self.assertIs(value['a01_complete'],False)
            self.assertIs(value['diagnostic_only'],True)
            self.assertIs(value['diagnostic_completed'],True)
        for key in ('production_rust_changed','nc_code_or_configuration_changed',
                    'acceptance_tolerance_changed','old_values_rewritten',
                    'new_support_or_regularization','physical_instability_proved'):
            self.assertIs(a[key],False)
        self.assertIs(a['refusals_retained'],True)
        self.assertEqual(m['new_global_solves'],0)
        self.assertEqual((m['element_count'],m['reduced_dofs']),(18,48))
        self.assertEqual(m['physical_instability_classification_d'],'判定不能')

    def test_seal_and_report_reproduce_without_arbitrary_precision(self):
        self.assertEqual(hashlib.sha256((HERE/'archive-index.json').read_bytes()).hexdigest(),INDEX_SHA)
        with tempfile.TemporaryDirectory() as temp:
            output=Path(temp)/'evidence'
            proof=verify(HERE,output)
            self.assertTrue(proof['all_sha_verified'])
            result=Path(temp)/'measurements.json'
            collect(output/'raw',result)
            self.assertEqual(result.read_bytes(),(HERE/'measurements.json').read_bytes())

    def test_saved_inputs_remain_bound_to_old_members(self):
        binding=read(HERE/'input-binding.json')
        self.assertIs(binding['saved_state_and_displacement_unchanged'],True)
        self.assertEqual(binding['new_global_solves'],0)
        index=read(HERE/'archive-index.json')
        for name,entry in binding['selected_members'].items():
            self.assertEqual(index['entries']['source/'+name],entry)
        old=ROOT/binding['old_archive_index_path']
        self.assertEqual(hashlib.sha256(old.read_bytes()).hexdigest(),binding['old_archive_index_sha256'])

    def test_zero_rank_observations_are_separate_from_acceptance(self):
        m=read(HERE/'measurements.json')
        self.assertIs(m['cutoff_is_acceptance_tolerance'],False)
        self.assertEqual(m['exact_analytic_null_directions'],2)
        self.assertEqual(m['exact_fd_null_directions'],1)
        self.assertEqual(m['fd_probe_branch_changes'],[])
        for name,value in m['matrices'].items():
            self.assertEqual(len(value['spectrum_50']),48)
            self.assertEqual(len(value['spectrum_80']),48)
            self.assertEqual(value['rank_50'],46 if name=='analytic' else 47)
            self.assertEqual(value['rank_80'],value['rank_50'])
            self.assertEqual(len(value['signed_projections']),2)

    def test_rhs_sign_active_holding_and_left_range_are_explicit(self):
        m=read(HERE/'measurements.json')
        self.assertEqual(m['residual_sign'],'external_minus_internal')
        self.assertIs(m['birth_holding_active_at_stop'],False)
        self.assertIs(m['frozen_linearized_rhs_is_in_range'],False)
        self.assertIs(m['roundoff_is_sole_stopping_mechanism'],False)
        for value in m['matrices'].values():
            for projected in value['signed_projections']:
                self.assertEqual(Decimal(projected['birth_holding_at_stop']['left_component']),0)
            self.assertGreater(Decimal(value['true_null_projection_norms']['saved_residual']),1)
        self.assertIs(m['physical_percentages_summing_to_100'],False)

    def test_documents_accept_actual_crlf_without_touching_sealed_bytes(self):
        for name in ('README.ja.md','REPRODUCE.ja.md'):
            text=(HERE/name).read_text(encoding='utf-8')
            self.assertTrue(text.endswith('A-01は未達、passed=false。\n'))
            with tempfile.TemporaryDirectory() as temp:
                p=Path(temp)/name;p.write_bytes(text.replace('\n','\r\n').encode('utf-8'))
                self.assertEqual(p.read_text(encoding='utf-8'),text)

    def test_no_production_diagnostic_import_or_switch(self):
        for path in (ROOT/'src').rglob('*.rs'):
            text=path.read_text(encoding='utf-8')
            self.assertNotIn('a01hp_',text)
            self.assertNotIn('A01HP_',text)
        for path in (ROOT/'requirements.lock.txt',ROOT/'requirements-reference.lock.txt'):
            self.assertNotIn('mpmath',path.read_text(encoding='utf-8'))

    def test_source_access_and_prior_records_are_not_relabelled(self):
        s=read(HERE/'primary-source-review.json')
        self.assertIs(s['book_body_claims_invented'],False)
        self.assertIs(s['external_reference_values_introduced'],False)
        prior=read(HERE/'prior-record-binding.json')
        for path,digest in prior['files'].items():
            self.assertEqual(hashlib.sha256((ROOT/path).read_bytes()).hexdigest(),digest)
        books=[v for v in s['sources'] if any(x in v['title'] for x in ('Wriggers','Laursen','Nocedal','Dennis'))]
        self.assertEqual(len(books),4)
        self.assertTrue(all('本文未取得' in v['access'] for v in books))

if __name__=='__main__':
    unittest.main()
