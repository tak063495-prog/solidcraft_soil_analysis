"""Opt-in global control: independent rollback, exact coordinates, unchanged legacy reference."""
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
HERE=ROOT/'interim/2026-12-11-c01-global-step-control'
spec=importlib.util.spec_from_file_location('c01_control_study_tests',HERE/'study.py')
st=importlib.util.module_from_spec(spec);spec.loader.exec_module(st)
spec=importlib.util.spec_from_file_location('c01_old_reference_tests',ROOT/'interim/2026-11-18-wet-joint-step-convergence/reference_steps.py')
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)

def inputs(n):
    mesh=st.read(ROOT/'examples/joint-water-column-tet10-mesh.json')
    ground=st.fem.normalize_ground(st.read(ROOT/'examples/ground-joint-water-mc.json'))
    case=st.prior.configuration(st.fem.normalize_stage_case(st.read(ROOT/'examples/stages-joint-water-cycle.json')),n,5e-5)
    case['settings']['minimum_increment']=min(case['settings']['minimum_increment'],1/(n*32))
    return mesh,ground,case

class GlobalStepControl(unittest.TestCase):
    def test_absent_control_keeps_every_independent_numerical_field_exact(self):
        mesh,ground,case=inputs(4)
        a=old._reference(mesh,ground,case,4,1,7,False,True,maximum_m=5e-5)
        b=st.reference._reference(mesh,ground,case,4,1,7,False,True,maximum_m=5e-5)
        for x,y in zip(a['frames'],b['frames']):
            for key in st.FIELDS: np.testing.assert_array_equal(x[key],y[key])
        self.assertEqual(a['maximum_relative_force_residual'],b['maximum_relative_force_residual'])

    def test_reference_owns_its_rollback_and_prescribed_coordinates(self):
        mesh,ground,case=inputs(8);before=copy.deepcopy((mesh,ground,case));fields=[]
        result=st.reference._reference(mesh,ground,case,8,1,7,False,True,lambda r,f:fields.append((r,f)),maximum_m=5e-5,switch_bisections=4)
        self.assertEqual((mesh,ground,case),before)
        actual={(r['stage'],r['progress']) for r in result['path']}
        self.assertTrue({(s,i/8) for s in st.STAGES for i in range(1,9)}<=actual)
        self.assertTrue(any(not a['accepted'] for a in result['attempts']))
        for a in result['attempts']:
            if not a['accepted']:
                self.assertLess(a['relative_force_residual'],1e-10)
                self.assertNotEqual(a['old_status'],a['trial_status'])
        for row,f in fields:
            self.assertTrue(all(np.isfinite(v).all() for key,v in f.items() if key!='name'))
        self.assertFalse(result['reference_uses_native_results'])

    def test_document_pc_command_and_decision_survive_crlf(self):
        # Text comparisons normalize line endings; binary evidence hashes remain strict.
        for name in ('README.ja.md','REPRODUCE.ja.md'):
            text=(HERE/name).read_text(encoding='utf-8')
            with tempfile.TemporaryDirectory() as tmp:
                path=Path(tmp)/name;path.write_bytes(text.replace('\r\n','\n').replace('\n','\r\n').encode('utf-8'))
                normalized=path.read_text(encoding='utf-8')
                self.assertEqual(normalized,text)
                self.assertIn('C-01',normalized)
                self.assertTrue(normalized.endswith('C-01は未達、passed=false。\n'))
        self.assertIn('--maximum 0.000025',(HERE/'REPRODUCE.ja.md').read_text(encoding='utf-8'))


class GlobalStepSavedEvidence(unittest.TestCase):
    def test_all_history_metrics_are_recomputed_and_forged_acceptance_is_refused(self):
        spec=importlib.util.spec_from_file_location('c01_global_saved_evidence_tests',HERE/'evidence.py')
        evidence=importlib.util.module_from_spec(spec);spec.loader.exec_module(evidence)
        with tempfile.TemporaryDirectory() as tmp:
            result=evidence.replay(evidence.extract(Path(tmp)/'replay'))
            self.assertTrue(result['measurements_exact']);self.assertFalse(result['numerical_acceptance'])
        spec=importlib.util.spec_from_file_location('c01_global_measure_tests',HERE/'measure.py')
        measure=importlib.util.module_from_spec(spec);spec.loader.exec_module(measure)
        report=json.loads((HERE/'measurements.json').read_text(encoding='utf-8'))
        for key in ('passed','c02_c03_c04_accepted','nc01_accepted','p6_complete','ab_items_completed'):
            changed=copy.deepcopy(report);changed[key]=True
            with self.assertRaisesRegex(ValueError,'unsupported'): measure.validate_report(changed)

if __name__=='__main__': unittest.main()
