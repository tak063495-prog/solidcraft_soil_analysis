"""Dedicated view transport contracts; mechanical verification is mocked here.

Actual Native CAD/solver/reader evidence is a separate integration study.
"""
import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
import export_mc_epoch_view as module
from test_mc_epoch_export import fixture


def setup(folder):
    root=Path(folder);job=root/'job';native=root/'origin'
    job.mkdir();native.mkdir()
    mesh,origin,checkpoint=fixture()
    metrics=dict(relative_force_residual=1e-12,maximum_normalized_yield=1e-13)
    origin['frames']=[dict(state=origin['last_committed'],**metrics)]
    checkpoint['payload']['continuation']=dict(accepted=[metrics],progress=.5)
    def write(base,name,value):
        path=base/name;path.parent.mkdir(parents=True,exist_ok=True)
        data=json.dumps(value,allow_nan=False).encode();path.write_bytes(data)
        return hashlib.sha256(data).hexdigest()
    job_hashes={}
    for n,v in (('inputs/mesh.json',mesh),('inputs/origin.json',origin),('inputs/specification.json',{})):
        job_hashes[n]=write(job,n,v)
    (job/'solver').mkdir();(job/'solver/soil-fem-3d.exe').write_bytes(b'not executable; mocked audits only')
    job_hashes['solver/soil-fem-3d.exe']=module.sha(job/'solver/soil-fem-3d.exe')
    job_manifest=dict(job_id='a'*64,binding=dict(files_sha256=job_hashes));write(job,'job.json',job_manifest)
    for source,target in (('inputs/mesh.json','analysis/ground/mesh.json'),('inputs/origin.json','analysis/stage-results.json'),('solver/soil-fem-3d.exe','analysis/solver/soil-fem-3d.exe')):
        path=native/target;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes((job/source).read_bytes())
    geometry=dict(schema='soil-fem-view-geometry/1',length_unit='m',axis='Z-up',nodes_m=mesh['nodes'],cells=mesh['tetrahedra'],cell_kind='tetra4')
    geometry_hash=write(native,'view/geometry.json',geometry)
    sources={n:module.sha(native/n) for n in ('analysis/ground/mesh.json','analysis/stage-results.json','analysis/solver/soil-fem-3d.exe')}
    manifest=dict(schema='soil-fem-solidcraft-view/1',complete=True,input_sha256='b'*64,cad_sha256='c'*64,exporter_binary_sha256='d'*64,
        solver_binary_sha256=job_hashes['solver/soil-fem-3d.exe'],geometry='geometry.json',files_sha256={'geometry.json':geometry_hash},job_sources_sha256=sources)
    digest=write(native,'view/view-manifest.json',manifest)
    write(native,'view-verification.json',dict(passed=True,input_sha256=manifest['input_sha256'],view_manifest_sha256=digest))
    verified=dict(progress=.5,complete=False,stop_reason=None)
    snapshot=json.dumps(dict(checkpoint=checkpoint)).encode()
    accepted=(mesh,origin,checkpoint,verified,snapshot,job_manifest)
    return root,job,native,accepted


class EpochView(unittest.TestCase):
    def test_path_traversal_and_windows_paths_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder).resolve()
            for name in ('../secret','x/../secret','x//y','/root','C:/root',r'\\server\share','x\\y','./x'):
                with self.assertRaises(ValueError,msg=name):module.source(root,name)
            self.assertEqual(module.source(root,'view/geometry.json'),root/'view/geometry.json')

    def test_new_separate_output_required_before_verification(self):
        with tempfile.TemporaryDirectory() as folder:
            root,job,native,accepted=setup(folder)
            with patch.object(module,'accepted') as audit:
                for output in (job,job/'child',native/'child',root):
                    with self.assertRaises(ValueError):module.export(job,native,output)
                audit.assert_not_called()

    def test_partial_snapshot_keeps_dedicated_schema_and_true_increment(self):
        with tempfile.TemporaryDirectory() as folder:
            root,job,native,accepted=setup(folder);out=root/'view-packet'
            before={str(q):module.sha(q) for base in (job,native) for q in base.rglob('*') if q.is_file()}
            with patch.object(module,'accepted',return_value=accepted):manifest=module.export(job,native,out)
            self.assertEqual(manifest['schema'],'soil-fem-solidcraft-mc-epoch-view/1')
            self.assertTrue(manifest['complete']);self.assertFalse(manifest['epoch']['path_complete'])
            self.assertEqual(manifest['epoch']['parameter_progress'],.5)
            frame=json.loads((out/'view/frame-001.json').read_text('utf-8'))
            self.assertEqual(frame['schema'],'soil-fem-view-frame/4')
            self.assertIn('epoch_displacement_increment_m',frame['node_fields'])
            self.assertNotIn('stage_displacement_increment_m',frame['node_fields'])
            self.assertAlmostEqual(frame['node_fields']['epoch_displacement_increment_m'][0][2],-.003)
            self.assertEqual(frame['cell_fields']['material_id'],[7])
            self.assertEqual(frame['cell_fields']['effective_stress_compression_positive_kpa'][0][0],5.)
            self.assertEqual((out/'epoch/accepted-snapshot.json').read_bytes(),accepted[4])
            self.assertEqual(before,{n:module.sha(Path(n)) for n in before})

    def test_mismatched_origin_run_rejected_before_output(self):
        with tempfile.TemporaryDirectory() as folder:
            root,job,native,accepted=setup(folder)
            manifest_path=native/'view/view-manifest.json';manifest=json.loads(manifest_path.read_text('utf-8'))
            manifest['job_sources_sha256']['analysis/stage-results.json']='f'*64
            manifest_path.write_text(json.dumps(manifest),'utf-8')
            audit_path=native/'view-verification.json';audit=json.loads(audit_path.read_text('utf-8'))
            audit['view_manifest_sha256']=module.sha(manifest_path);audit_path.write_text(json.dumps(audit),'utf-8')
            with patch.object(module,'accepted',return_value=accepted):
                with self.assertRaisesRegex(ValueError,'origin mesh/Run/executable'):module.export(job,native,root/'output')
            self.assertFalse((root/'output').exists())

    def test_unverified_native_origin_rejected_before_output(self):
        with tempfile.TemporaryDirectory() as folder:
            root,job,native,accepted=setup(folder)
            audit=native/'view-verification.json';value=json.loads(audit.read_text('utf-8'));value['passed']=False;audit.write_text(json.dumps(value),'utf-8')
            with patch.object(module,'accepted',return_value=accepted):
                with self.assertRaisesRegex(ValueError,'origin view audit'):module.export(job,native,root/'output')
            self.assertFalse((root/'output').exists())

    def test_altered_native_source_checksum_rejected_before_output(self):
        with tempfile.TemporaryDirectory() as folder:
            root,job,native,accepted=setup(folder)
            with (native/'analysis/ground/mesh.json').open('ab') as stream:stream.write(b' ')
            with patch.object(module,'accepted',return_value=accepted):
                with self.assertRaisesRegex(ValueError,'source/packet checksum'):module.export(job,native,root/'output')
            self.assertFalse((root/'output').exists())


if __name__=='__main__':unittest.main()
