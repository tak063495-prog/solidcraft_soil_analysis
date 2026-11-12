"""Isolated mutation tests; no solver, skip, xfail or existing collection changes."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('p6_status_audit', HERE/'audit_status.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class AuditMutationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.payload = b'{"passed": false, "accepted": true, "a/b": {"~x": [false]}}\n'
        self.root.joinpath('proof.json').write_bytes(self.payload)
        for path in audit.DOCS:
            p = self.root/path
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text('test\n', encoding='utf-8', newline='')
        self.ledger = {'base_commit': audit.BASE, 'sources': {}, 'evidence': {}, 'rows': [], 'claims': [], 'ci_records': []}
        for i, path in enumerate(['proof.json']+audit.DOCS):
            data = (self.root/path).read_bytes()
            self.ledger['sources'][f's{i}'] = dict(path=path, commit=audit.BASE, sha256=audit.sha(data),
                git_blob_sha1=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest(),
                kind='json' if i == 0 else 'markdown', status_keys=[])
        self.ledger['evidence'] = {
            'no': dict(source='s0', pointer='/passed', value=False, supports='未達', scope='all'),
            'yes': dict(source='s0', pointer='/accepted', value=True, supports='限定受入', scope='small'),
            'quote': dict(source='s1', lines=[1,1], value='test')}
        self.ledger['rows'] = [dict(id=t, target=t, state='未達', scope='all', reason='not done',
                                  evidence=['no'], completed=False) for t in sorted(audit.REQUIRED)]

    def codes(self):
        return {f['code'] for f in audit.audit(self.root, self.ledger, verify_git=False)['findings']}

    def test_valid_evidence(self):
        self.assertEqual(self.codes(), set())

    def test_missing_file(self):
        (self.root/'proof.json').unlink()
        self.assertIn('SOURCE_MISSING', self.codes())

    def test_sha_mutation(self):
        (self.root/'proof.json').write_bytes(self.payload+b' ')
        self.assertIn('SOURCE_SHA256', self.codes())

    def test_blob_mutation(self):
        self.ledger['sources']['s0']['git_blob_sha1'] = '0'*40
        self.assertIn('SOURCE_BLOB', self.codes())

    def test_wrong_base(self):
        self.ledger['base_commit'] = '0'*40
        self.assertIn('BASE_COMMIT', self.codes())

    def test_wrong_source_commit(self):
        self.ledger['sources']['s0']['commit'] = '0'*40
        self.assertIn('SOURCE_COMMIT', self.codes())

    def test_key_value_mutation(self):
        self.ledger['evidence']['no']['value'] = True
        self.assertIn('EVIDENCE_VALUE', self.codes())

    def test_boolean_is_not_integer(self):
        self.ledger['evidence']['no']['value'] = 0
        self.assertIn('EVIDENCE_VALUE', self.codes())

    def test_missing_pointer(self):
        self.ledger['evidence']['no']['pointer'] = '/missing'
        self.assertIn('EVIDENCE_SELECTOR', self.codes())

    def test_rfc6901_arrays(self):
        self.ledger['evidence']['no']['pointer'] = '/a~1b/~0x/0'
        self.assertEqual(self.codes(), set())

    def test_bad_pointer_escape(self):
        self.ledger['evidence']['no']['pointer'] = '/bad~2key'
        self.assertIn('EVIDENCE_SELECTOR', self.codes())

    def test_quote_mutation(self):
        self.ledger['evidence']['quote']['value'] = 'done'
        self.assertIn('EVIDENCE_VALUE', self.codes())

    def test_quote_range(self):
        self.ledger['evidence']['quote']['lines'] = [0,1]
        self.assertIn('EVIDENCE_SELECTOR', self.codes())

    def test_unsupported_completion(self):
        self.ledger['rows'][0]['completed'] = True
        self.assertIn('UNSUPPORTED_COMPLETION', self.codes())

    def test_acceptance_without_positive_evidence(self):
        self.ledger['rows'][0]['state'] = '限定受入'
        self.assertIn('UNSUPPORTED_ACCEPTANCE', self.codes())

    def test_scope_broadening(self):
        self.ledger['rows'][0].update(state='限定受入', scope='general 3D', evidence=['yes'])
        self.assertIn('SCOPE_EXPANSION', self.codes())

    def test_limited_acceptance(self):
        self.ledger['rows'][0].update(state='限定受入', scope='small', evidence=['yes'])
        self.assertEqual(self.codes(), set())

    def test_same_scope_document_conflict(self):
        self.ledger['claims'] = [dict(target='X', scope='small', state=s, evidence=e)
                                 for s, e in [('限定受入','yes'), ('未達','no')]]
        self.assertIn('DOCUMENT_CONTRADICTION', self.codes())

    def test_different_scope_is_not_conflict(self):
        self.ledger['claims'] = [dict(target='X', scope='small', state='限定受入', evidence='yes'),
                                 dict(target='X', scope='all', state='未達', evidence='no')]
        self.assertNotIn('DOCUMENT_CONTRADICTION', self.codes())

    def test_historical_conflict_reported(self):
        self.ledger['claims'] = [dict(target='X', scope='small', state='限定受入', evidence='yes'),
                                 dict(target='X', scope='small', state='未達', evidence='no', historical=True)]
        self.assertIn('HISTORICAL_SCOPE_CONFLICT', self.codes())

    def test_cancelled_not_success(self):
        self.ledger['ci_records'] = [dict(id=1, status='completed', conclusion='cancelled', claimed_success=True)]
        self.assertIn('CI_NOT_SUCCESS', self.codes())

    def test_timeout_not_success(self):
        self.ledger['ci_records'] = [dict(id=1, status='completed', conclusion='timed_out', claimed_success=True)]
        self.assertIn('CI_NOT_SUCCESS', self.codes())

    def test_queued_not_success(self):
        self.ledger['ci_records'] = [dict(id=1, status='queued', conclusion=None, claimed_success=True)]
        self.assertIn('CI_NOT_SUCCESS', self.codes())

    def test_wrong_ci_head(self):
        self.ledger['ci_records'] = [dict(id=1, status='completed', conclusion='success', claimed_success=True,
                                        head_sha='old', claimed_head_sha='new')]
        self.assertIn('CI_WRONG_HEAD', self.codes())

    def test_coverage(self):
        self.ledger['rows'].pop()
        self.assertIn('COVERAGE_MISSING', self.codes())

    def test_new_source_not_inventoried(self):
        p = self.root/'interim/new/acceptance.json'
        p.parent.mkdir(parents=True)
        p.write_text('{}')
        self.assertIn('SOURCE_NOT_INVENTORIED', self.codes())

    def test_path_traversal(self):
        self.ledger['sources']['s0']['path'] = '../proof.json'
        self.assertIn('SOURCE_MISSING', self.codes())

    def test_inventory_value(self):
        self.ledger['sources']['s0']['status_keys'] = [dict(pointer='/passed', value=True)]
        self.assertIn('INVENTORY_KEY_VALUE', self.codes())

    def test_duplicate_id(self):
        self.ledger['rows'].append(copy.deepcopy(self.ledger['rows'][0]))
        self.assertIn('DUPLICATE_ROW', self.codes())


    def test_row_same_scope_conflict(self):
        other = copy.deepcopy(self.ledger['rows'][0])
        other.update(id='other', state='限定受入', scope='all', evidence=['yes'])
        self.ledger['rows'].append(other)
        self.assertIn('ROW_STATE_CONTRADICTION', self.codes())

    def test_updated_hashes_cannot_reseal_changed_base_bytes(self):
        changed = self.payload.replace(b'false', b'true', 1)
        (self.root/'proof.json').write_bytes(changed)
        source = self.ledger['sources']['s0']
        source['sha256'] = audit.sha(changed)
        source['git_blob_sha1'] = hashlib.sha1(b'blob '+str(len(changed)).encode()+b'\0'+changed).hexdigest()
        self.ledger['evidence']['no']['value'] = True
        def original_bytes(root, *args):
            if args[0] in ('diff', 'ls-files'):
                return b''
            original = self.payload if args[-1].endswith(':proof.json') else b'test\n'
            return hashlib.sha1(b'blob '+str(len(original)).encode()+b'\0'+original).hexdigest().encode()+b'\n'
        with patch.object(audit, 'git', side_effect=original_bytes):
            result = audit.audit(self.root, self.ledger)
        self.assertIn('BASE_BYTES', {f['code'] for f in result['findings']})

    def test_document_status_column_completion(self):
        (self.root/audit.DOCS[0]).write_text('|P6-01|all|完了（全体）|no proof|\n', encoding='utf-8', newline='')
        self.assertIn('UNSUPPORTED_COMPLETION', self.codes())

    def test_json_p6_completion(self):
        (self.root/'proof.json').write_text('{"passed":false,"accepted":true,"P6_complete":true}')
        self.assertIn('UNSUPPORTED_COMPLETION', self.codes())

    def test_negated_completion_does_not_trigger(self):
        (self.root/audit.DOCS[0]).write_text('|P6-01|all|未完了|no proof|\n', encoding='utf-8', newline='')
        self.assertNotIn('UNSUPPORTED_COMPLETION', self.codes())

    def test_infra_checkout_failure_cannot_be_success(self):
        self.ledger['ci_records'] = [dict(id=43191750591, status='completed', conclusion='success',
            infrastructure_failure=True, claimed_success=True, head_sha='a'*40, claimed_head_sha='a'*40,
            jobs=[dict(name='python-spatial-evidence', conclusion='failure')])]
        self.assertIn('CI_INFRA_OR_JOB_FAILURE', self.codes())

    def test_false_cannot_be_acceptance(self):
        self.ledger['evidence']['no'].update(supports='限定受入', scope='all')
        self.ledger['rows'][0].update(state='限定受入')
        self.assertIn('FALSE_OR_NULL_ACCEPTANCE', self.codes())

    def test_overall_cannot_inherit_small_acceptance(self):
        row = next(r for r in self.ledger['rows'] if r['target'] == '一般3D空間収束')
        row.update(state='限定受入', scope_kind='overall', scope='small', evidence=['yes'])
        self.assertIn('SCOPE_EXPANSION', self.codes())

    def test_exact_document_appendix_allowed(self):
        appendix = '\n監査追補: 未達維持\n'
        self.ledger['sources']['s1']['appendix'] = appendix
        (self.root/audit.DOCS[0]).write_text('test\n'+appendix, encoding='utf-8', newline='')
        self.assertEqual(self.codes(), set())

    def test_document_appendix_mutation(self):
        self.ledger['sources']['s1']['appendix'] = '\n未達維持\n'
        (self.root/audit.DOCS[0]).write_text('test\n完了\n', encoding='utf-8', newline='')
        self.assertIn('DOCUMENT_APPENDIX', self.codes())

    def test_failed_job_even_with_success_run(self):
        self.ledger['ci_records'] = [dict(id=1, status='completed', conclusion='success',
            claimed_success=True, head_sha='a'*40, claimed_head_sha='a'*40,
            jobs=[dict(name='checks', conclusion='cancelled')])]
        self.assertIn('CI_INFRA_OR_JOB_FAILURE', self.codes())

    def test_unsupported_prose_completion(self):
        (self.root/audit.DOCS[0]).write_text('P6全体は完了しました。\n', encoding='utf-8', newline='')
        self.assertIn('UNSUPPORTED_COMPLETION', self.codes())

    def test_nc_diagnostic_is_not_nc_completion(self):
        (self.root/'proof.json').write_text('{"passed":false,"accepted":true,"nc01_complete":true}')
        self.assertIn('UNSUPPORTED_COMPLETION', self.codes())

    def test_approved_appendix_still_cannot_claim_completion(self):
        appendix = '\nP6全体は完了しました。\n'
        self.ledger['sources']['s1']['appendix'] = appendix
        (self.root/audit.DOCS[0]).write_text('test\n'+appendix, encoding='utf-8', newline='')
        self.assertIn('UNSUPPORTED_COMPLETION', self.codes())

    def test_git_source_cannot_be_relabelled_as_external(self):
        self.ledger['sources']['s0']['origin'] = 'github_snapshot'
        self.assertIn('EXTERNAL_SOURCE_SCOPE', self.codes())


class FixedRepositoryTests(unittest.TestCase):
    def test_ledger_and_table(self):
        ledger = json.loads((HERE/'ledger.json').read_text(encoding='utf-8'))
        self.assertEqual((HERE/'ledger.ja.md').read_text(encoding='utf-8'), audit.table(ledger))
        self.assertEqual({r['state'] for r in ledger['rows']}, set(audit.STATES))
        self.assertTrue(all(r['completed'] is False for r in ledger['rows']))

    def test_all_fixed_sources_and_classifications(self):
        ledger = json.loads((HERE/'ledger.json').read_text(encoding='utf-8'))
        result = audit.audit(HERE.parents[1], ledger)
        self.assertTrue(result['integrity_ok'], result['findings'])
        self.assertFalse(result['numerical_acceptance_performed'])
        self.assertIn('HISTORICAL_SCOPE_CONFLICT', result['finding_counts'])


if __name__ == '__main__':
    unittest.main(verbosity=2)
