"""End-to-end mixed primary fields and structural history audit falsification."""
import os
for name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[name]='1'
import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
import mixed_mcc_io
import fem_reference as ref
from prepare_mcc_structures import prepare
from verify_projected_stage_package import audit


class MixedStructuralAudit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root=Path(__file__).resolve().parents[1]
        binary=root/'target/release'/('soil-fem-3d.exe' if os.name=='nt' else 'soil-fem-3d')
        if not binary.is_file():raise unittest.SkipTest('build the release CLI before integration audits')
        cls.temporary=tempfile.TemporaryDirectory();cls.addClassCleanup(cls.temporary.cleanup);cls.runs=[]
        for projected in (True,False):
            folder=Path(cls.temporary.name)/str(projected);prepare(root,folder,projected)
            paths=[folder/'ground/mesh.json',folder/'ground/ground.json',folder/'stages.json']
            result=folder/'stage-results.json'
            subprocess.run([str(binary),'solve-stages',*map(str,paths),str(result)],check=True,capture_output=True,text=True)
            mesh,ground,case=[json.loads(p.read_text()) for p in paths]
            cls.runs.append((mesh,ground,case,json.loads(result.read_text())))

    def test_both_primary_volume_forms_audit_member_loads_excavation_and_installation(self):
        for mesh,ground,case,run in self.runs:
            with self.subTest(schema=run['schema']):
                self.assertEqual(len(audit(mesh,ground,case,run)),9)
                self.assertTrue(any(not a['accepted'] for a in run['attempts']))
                self.assertTrue(any(any(abs(x)>1e-6 for x in p['plastic_strain_compression_positive']) for p in run['frames'][4]['state']['material_states']))

    def test_downgraded_schema_and_missing_structural_report_are_rejected(self):
        for mesh,_,_,original in self.runs:
            for failure in ('schema','report'):
                run=copy.deepcopy(original)
                if failure=='schema':run['schema']=run['schema'].replace('/2','/1')
                else:run['frames'][2].pop('beam_results')
                with self.subTest(schema=original['schema'],failure=failure),self.assertRaises(ValueError):mixed_mcc_io.validate(mesh,run)

    def test_modified_member_moment_is_rejected(self):
        for mesh,_,_,original in self.runs:
            run=copy.deepcopy(original);run['frames'][3]['beam_results']['external_moment_origin_kn_m'][0]+=.01
            with self.assertRaises((ValueError,AssertionError)):mixed_mcc_io.validate(mesh,run)

    def test_modified_installation_reference_and_holding_vector_are_rejected(self):
        for mesh,ground,case,original in self.runs:
            for failure in ('reference','holding'):
                run=copy.deepcopy(original);transition=run['transition_audits'][0]
                if failure=='reference':transition['prepared_state']['beams']['plate_reference_global_kinematics'][0][0]+=.001
                else:transition['initial_holding_load_kn'][0]+=.1
                with self.subTest(schema=original['schema'],failure=failure),self.assertRaises(ValueError):audit(mesh,ground,case,run)

    def test_modified_excavation_release_is_rejected(self):
        for mesh,ground,case,original in self.runs:
            run=copy.deepcopy(original);run['excavation_audits'][0]['initial_holding_load_kn'][0]+=.1
            with self.assertRaises(ValueError):audit(mesh,ground,case,run)

    def test_geostatic_base_load_is_counted_once_with_both_structural_forms(self):
        root=Path(__file__).resolve().parents[1]
        binary=root/'target/release'/('soil-fem-3d.exe' if os.name=='nt' else 'soil-fem-3d')
        for original_mesh,_,original_case,_ in self.runs:
            mesh=copy.deepcopy(original_mesh);case=copy.deepcopy(original_case)
            mesh['face_sets']={'top':[list(face) for face,owners in ref.topology(mesh).items()
                if len(owners)==1 and all(mesh['nodes'][n][2]==1. for n in face)]}
            ground=json.loads((root/'examples/ground-geostatic-mcc-column.json').read_text())
            case['stages']=case['stages'][1:3];case['verification_reject_once_stages']=['install']
            with tempfile.TemporaryDirectory() as tmp:
                paths=[Path(tmp)/name for name in ('mesh.json','ground.json','stages.json')]
                for path,value in zip(paths,(mesh,ground,case),strict=True):path.write_text(json.dumps(value))
                result=Path(tmp)/'result.json'
                subprocess.run([str(binary),'solve-stages',*map(str,paths),str(result)],check=True,capture_output=True,text=True)
                self.assertEqual(len(audit(mesh,ground,case,json.loads(result.read_text()))),3)


if __name__=='__main__':unittest.main()
