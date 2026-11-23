"""Bounded C-03 diagnostics: deterministic arithmetic and sealed observations."""
import hashlib
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
STUDY=ROOT/'interim/2026-12-18-c03-plastic-bulk-cause-diagnosis'
sys.path.insert(0,str(ROOT/'tools'))
sys.path.insert(0,str(STUDY))
import c03diag_reference
from c03diag_metrics import norm,dot,metric,compare
from c03diag_mcc_point import step


class C03PlasticBulkCauseDiagnosis(unittest.TestCase):
    def test_fixed_reduction_order(self):
        self.assertEqual(dot([1e16,1.,-1e16],[1.,1.,1.]),1.)
        self.assertEqual(norm([[3.,4.],[0.,12.]]),13.)
        self.assertEqual(metric([0.,0.],[3.,4.])['relative_l2'],1.)

    def test_invalid_correspondence_and_nonfinite(self):
        with self.assertRaises(ValueError):metric(np.zeros(2),np.zeros(3))
        with self.assertRaises(ValueError):dot([1.],[1.,2.])
        for bad in (float('nan'),float('inf')):
            with self.assertRaises(ValueError):metric([bad],[0.])

    def test_endpoint_gate_cannot_hide_history(self):
        rows=[dict(stage='cycle',progress=t) for t in (.5,1.)]
        zeros=np.zeros((2,1,3));changed=zeros.copy();changed[0,0,0]=1.
        a=(rows,dict(u=changed,traction=changed));b=(rows,dict(u=zeros,traction=zeros))
        result=compare(a,b)
        self.assertEqual(result['all_fields']['u']['absolute_l2'],1.)
        self.assertEqual(result['endpoints']['u']['absolute_l2'],0.)
        self.assertFalse(result['numerical_acceptance'])
        other=[dict(stage='cycle',progress=t) for t in (.25,1.)]
        with self.assertRaises(ValueError):compare(a,(other,b[1]))

    def test_mcc_independent_hydrostatic_return_and_unloading(self):
        mat=dict(critical_state_m=1.2,**{'lambda':.2},kappa=.04,poisson=.3)
        old=dict(effective_stress_compression_positive_kpa=[100.,100.,100.,0.,0.,0.],
                 void_ratio=.7,preconsolidation_pressure_kpa=100.,plastic_strain_compression_positive=[0.]*6)
        nxt=step(mat,old,np.array([.001,.001,.001,0.,0.,0.]))
        p=nxt['effective_stress_compression_positive_kpa'][0];pc=nxt['preconsolidation_pressure_kpa']
        self.assertAlmostEqual(p,pc,places=9)
        dp=sum(nxt['plastic_strain_compression_positive'][:3])
        self.assertAlmostEqual(math.log(pc/100),1.7/(.2-.04)*dp,places=11)
        unloaded=step(mat,nxt,np.array([-.0001,-.0001,-.0001,0.,0.,0.]))
        self.assertEqual(unloaded['preconsolidation_pressure_kpa'],pc)
        self.assertEqual(unloaded['plastic_strain_compression_positive'],nxt['plastic_strain_compression_positive'])
        self.assertEqual(old['preconsolidation_pressure_kpa'],100.)

    def test_document_crlf_and_lf_have_same_terminal_assertion(self):
        text=(STUDY/'README.ja.md').read_text(encoding='utf-8')
        self.assertTrue(text.endswith('C-03は未達、passed=false。\n'))
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'readme.md'
            for newline in ('\n','\r\n'):
                p.write_bytes(text.replace('\n',newline).encode('utf-8'))
                self.assertEqual(p.read_text(encoding='utf-8'),text)

    def test_sealed_metric_fixture_exact_without_platform_tolerance(self):
        fixture=json.loads((STUDY/'diagnostic-fixture.json').read_text(encoding='utf-8'))
        for record in fixture['comparisons']:
            self.assertEqual(metric(record['a'],record['b']),record['expected'])
        self.assertGreater(fixture['native_contact_adjacent_l2'][1],fixture['native_contact_adjacent_l2'][0])
        self.assertGreater(fixture['independent_contact_adjacent_l2'][1],fixture['independent_contact_adjacent_l2'][0])

    def test_sealed_parts_and_index(self):
        acceptance=json.loads((STUDY/'acceptance.json').read_text(encoding='utf-8'))
        raw=(STUDY/'archive-index.json').read_bytes()
        self.assertEqual(hashlib.sha256(raw).hexdigest(),acceptance['archive_index_sha256'])
        self.assertEqual(hashlib.sha256((STUDY/'diagnostic-fixture.json').read_bytes()).hexdigest(),acceptance['diagnostic_fixture_sha256'])
        self.assertEqual(hashlib.sha256((STUDY/'return-trace-summary.json').read_bytes()).hexdigest(),acceptance['return_trace_summary_sha256'])
        with (STUDY/'return-trace-evidence.zip').open('rb') as stream:
            self.assertEqual(hashlib.file_digest(stream,'sha256').hexdigest(),acceptance['return_trace_archive_sha256'])
        index=json.loads(raw)
        for part in index['parts']:
            p=STUDY/part['name'];self.assertEqual(p.stat().st_size,part['size'])
            with p.open('rb') as stream:
                self.assertEqual(hashlib.file_digest(stream,'sha256').hexdigest(),part['sha256'])

    def test_acceptance_and_scope_stay_unachieved(self):
        a=json.loads((STUDY/'acceptance.json').read_text(encoding='utf-8'))
        for key in ('passed','c03_complete','c01_complete','c02_complete','c04_complete','nc01_complete',
                    'p6_complete','ab_items_completed','pc_executed','original_acceptance_all_met'):
            self.assertIs(a[key],False,key)
        self.assertTrue(a['production_defaults_unchanged'])
        self.assertTrue(a['old_evidence_unchanged'])


class PrincipalOrderCutback(unittest.TestCase):
    def run_shear(self, failure, *, minimum=.4, cutback=True, attempts=None):
        attempts = [] if attempts is None else attempts
        system = SimpleNamespace(points=[], c=np.zeros((0, 3, 1)),
                                 weight=np.zeros(0))
        seed = {"effective_stress_compression_positive_kpa": np.zeros(6)}
        strain = np.zeros((6, 1))
        strain[0, 0] = 1.
        model = SimpleNamespace(
            system=system, nd=1, seeds=[seed],
            case=dict(stages=[dict(name="shear", target_factor=1.)],
                      settings=dict(minimum_increment=minimum)),
            t=np.ones((1, 1)), base=np.zeros(1), seed_force=np.zeros(1),
            boundary=np.ones(1), mats=[None], b=[strain], ids=[[0]],
            weights=[1.], cr=np.zeros((0, 3, 1)),
            strains=lambda u: (strain @ u)[None, :],
            internal=lambda states, traction: -states[0]["effective_stress_compression_positive_kpa"][:1],
        )
        load = 0.

        def offset(factor):
            nonlocal load
            load = factor
            return np.zeros(1)

        model.offset = offset

        def evaluate(material, old, delta, count, *, require_order):
            original = old["effective_stress_compression_positive_kpa"]
            displacement = -delta[0]
            if load == .5 and original[0] == 0 and displacement > 1e-7:
                reason = failure(displacement, require_order)
                if reason is not None:
                    raise ArithmeticError(reason)
            return {"effective_stress_compression_positive_kpa": original + delta}

        with patch.object(c03diag_reference, "integrate", side_effect=evaluate):
            result = c03diag_reference.run(
                model, steps=2, fixed_material=1, fixed_contact=1,
                cutback=cutback, attempt_sink=attempts.append)
        self.assertTrue(result["complete"])
        return [a for a in attempts if a["stage"] == "shear"]

    def test_cutback_retries_same_target_when_every_candidate_has_wrong_order(self):
        trials = []

        def failure(displacement, ordered):
            trials.append(ordered)
            return None if ordered else ("inadmissible principal order",)

        actual = self.run_shear(failure)
        # The accepted candidate is evaluated again at the next Newton iteration.
        self.assertEqual(trials, [False] * 18 + [True] * 2)
        self.assertEqual([(a["target"], a["accepted"]) for a in actual],
                         [(.5, False), (.5, True), (1., True)])
        self.assertIn("inadmissible principal order", actual[0]["reason"])
        self.assertTrue(actual[1]["ordered_return"])

    def test_cutback_tries_smaller_damping_before_retrying_the_target(self):
        trials = []

        def failure(displacement, ordered):
            trials.append((displacement, ordered))
            return ("inadmissible principal order",) if len(trials) == 1 else None

        actual = self.run_shear(failure)
        self.assertEqual(trials[:2], [(.5, False), (.25, False)])
        self.assertEqual([(a["target"], a["ordered_return"]) for a in actual],
                         [(.5, False), (1., False)])

    def test_cutback_bisects_only_after_ordered_trial_fails(self):
        calls = []

        def failure(displacement, ordered):
            calls.append(ordered)
            return ("inadmissible principal order",)

        actual = self.run_shear(failure, minimum=.2)
        self.assertEqual(calls, [False] * 18 + [True] * 18)
        self.assertEqual([(a["target"], a["accepted"]) for a in actual],
                         [(.5, False), (.5, False), (.25, True), (.5, True), (1., True)])
        self.assertTrue(actual[2]["ordered_return"])

    def test_default_line_search_keeps_its_existing_failure(self):
        with self.assertRaisesRegex(ArithmeticError, "independent line search"):
            self.run_shear(lambda *_: ("inadmissible principal order",), cutback=False)

    def test_cutback_does_not_disguise_unexpected_arithmetic_failure(self):
        with self.assertRaisesRegex(ArithmeticError, "true linear residual"):
            self.run_shear(lambda *_: ("true linear residual", 1e-4))
