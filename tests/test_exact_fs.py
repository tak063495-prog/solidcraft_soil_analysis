"""Real-section preparation never solves; only actual evidence may pass."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('real_srm', ROOT / 'tools/prepare_real_section_srm.py')
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)
BINARY = ROOT / 'target/release' / ('soil-fem-3d.exe' if os.name == 'nt' else 'soil-fem-3d')
FIXTURE = ROOT / 'examples/exact-fs-uniaxial'


class RealSectionPreparation(unittest.TestCase):
    def prepare(self, directory):
        source = Path(directory) / 'source.txt'
        source.write_text('試験用の一様応力立方体。実断面ではない。', encoding='utf-8')
        bundle = Path(directory) / 'bundle'
        result = MODULE.prepare(BINARY, *(FIXTURE / (n + '.json') for n in
                                ('mesh', 'ground', 'case', 'fresh-spec')),
                                source, 'TEST ONLY: not a real section', bundle)
        return bundle, result

    def test_prepare_inspects_inputs_and_creates_job_without_solving_or_passing(self):
        with tempfile.TemporaryDirectory() as d:
            bundle, result = self.prepare(d)
            for flag in ('search_executed_by_preparer', 'exact_fs_for_this_case', 'real_section_srm_complete'):
                self.assertIs(result[flag], False)
            self.assertFalse((bundle / 'job/report.json').exists())
            self.assertFalse((bundle / 'job/run.claim').exists())
            with self.assertRaises((ValueError, subprocess.CalledProcessError)):
                MODULE.verify(bundle, Path(d) / 'assessment.json')
            with self.assertRaises(FileExistsError):
                self.prepare(d)

    def test_missing_source_empty_identity_or_invalid_input_cannot_prepare(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(ValueError):
                MODULE.prepare(BINARY, *(FIXTURE / (n + '.json') for n in
                               ('mesh', 'ground', 'case', 'fresh-spec')), BINARY, '', Path(d) / 'empty')
            with self.assertRaises(FileNotFoundError):
                MODULE.prepare(BINARY, *(FIXTURE / (n + '.json') for n in
                               ('mesh', 'ground', 'case', 'fresh-spec')), Path(d) / 'absent', 'id', Path(d) / 'missing')
            bad = Path(d) / 'bad.json'
            bad.write_text('{}', encoding='utf-8')
            with self.assertRaises(subprocess.CalledProcessError):
                MODULE.prepare(BINARY, FIXTURE / 'mesh.json', FIXTURE / 'ground.json',
                               FIXTURE / 'case.json', bad, BINARY, 'id', Path(d) / 'invalid')

    def test_tampered_preparation_source_and_job_are_refused(self):
        for name in ('source', 'job', 'claim'):
            with self.subTest(name=name), tempfile.TemporaryDirectory() as d:
                bundle, _ = self.prepare(d)
                if name == 'source':
                    (bundle / 'section-source').write_text('changed', encoding='utf-8')
                elif name == 'job':
                    (bundle / 'job/job.json').write_text('{}', encoding='utf-8')
                else:
                    m = json.loads((bundle / 'preparation.json').read_text(encoding='utf-8'))
                    m['exact_fs_for_this_case'] = True
                    (bundle / 'preparation.json').write_text(json.dumps(m), encoding='utf-8')
                with self.assertRaises(ValueError):
                    MODULE.verify(bundle, Path(d) / 'assessment.json')

    def test_real_search_and_cli_verification_with_missing_or_modified_evidence(self):
        with tempfile.TemporaryDirectory() as d:
            bundle, _ = self.prepare(d)
            job = bundle / 'job'
            solver = job / 'solver/soil-fem-3d.exe'
            subprocess.run([str(solver), 'run-srm-search-job', str(job)], check=True, capture_output=True)
            result = MODULE.verify(bundle, Path(d) / 'assessment.json')
            self.assertTrue(result['exact_fs_for_this_case'])
            self.assertLess(result['attained_equilibrium_factor'], 4 / 3)
            self.assertGreater(result['independently_excluded_factor'], 4 / 3)
            status = json.loads((job / 'status.json').read_text(encoding='utf-8'))
            r = json.loads((job / 'report.json').read_text(encoding='utf-8'))
            r['trials'][0]['endpoint_audit'] = None
            (job / 'report.json').write_text(json.dumps(r), encoding='utf-8')
            with self.assertRaises((ValueError, subprocess.CalledProcessError)):
                MODULE.verify(bundle, Path(d) / 'bad-assessment.json')
            checked = subprocess.run([str(solver), 'verify-exact-fs',
                *(str(job / 'inputs' / (name + '.json')) for name in ('mesh', 'ground', 'case', 'specification')),
                str(job / 'report.json'), status['report_sha256'], str(Path(d) / 'raw-invalid.json')],
                capture_output=True)
            self.assertEqual(checked.returncode, 1)
            self.assertFalse((Path(d) / 'raw-invalid.json').exists())
