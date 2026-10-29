"""Display-contract tests with synthetic records, not mechanical validation."""
import copy
import json
import shutil
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tools'))
from export_solidcraft_view import EDGES,boundary,export,scene_frame
import export_solidcraft_view as view
from verify_acceptance_packet import sha


def write(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,allow_nan=False),encoding='utf-8')


def fixture():
    mesh=dict(nodes=[[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],[0.,0.,-1.]],
        tetrahedra=[[0,1,2,3],[0,2,1,4]],material_ids=[1,1],region_ids=[1,1],regions=[dict(id=1,body='Clay')])
    states=[dict(effective_stress_compression_positive_kpa=[p,p,p,0.,0.,0.],
        plastic_strain_compression_positive=[0.]*6,void_ratio=e,preconsolidation_pressure_kpa=300.) for p,e in [(10.,1.),(30.,3.),(50.,5.)]]
    a=dict(displacement_m=[[0.,0.,-.1]]*5,material_states=states,active_elements=[True,True])
    b=copy.deepcopy(a);b.update(displacement_m=[[0.,0.,-.5]]*5,active_elements=[True,False],display_reference_displacement_m=[[0.,0.,-.4]]*5)
    def frame(name,state):return dict(name=name,state=state,monitors={},relative_force_residual=0.,maximum_normalized_yield=-.1)
    layout=[dict(element_index=e,position_m=[.1,.1,z],weight_m3=w) for e,z,w in [(0,.1,1.),(0,.2,3.),(1,-.1,2.)]]
    run=dict(schema='soil-fem-stage-result/6',complete=True,failure_reason=None,frames=[frame('seed',a),frame('cut',b)],
        integration_points=layout,reference_displacement_m=[[0.,0.,-.05]]*5,case=dict(stages=[dict(kind='excavate')]),last_committed=b)
    return mesh,run


def job_fixture(root, specimen=None):
    """Minimal transport fixture. No solver or CAD binary is executed."""
    mesh,run=fixture() if specimen is None else specimen;analysis=root/'analysis';conformal=root/'conformal';conformal.mkdir(parents=True)
    ground=dict(schema='soil-fem-ground/1',units='m-kN-kPa',initialization=dict(method='equilibrate'),materials=[dict(id=1)],bodies=[dict(body='Clay',material_id=1)])
    stages=dict(schema='soil-fem-stages/1',settings=dict(relative_force_tolerance=1e-8),stages=[dict(kind='excavate')])
    if specimen is not None:
        stages=dict(schema='soil-fem-stages/1',settings=dict(relative_force_tolerance=1e-8),stages=[],beam_system=run['case']['beam_system'])
        run['case']=stages
    payload=dict(schema='solidcraft-soil-snapshot/1',cad_sha256='c'*64,cad_length_unit='mm',analysis_units='m-kN-kPa',
        ground=ground,stages=stages,bindings=[dict(analysis_body='Clay',cad_body='CAD layer')])
    (conformal/'body-0001.step').write_text('Transport-test placeholder only; never imported as geometry',encoding='utf-8')
    write(conformal/'ground.json',ground);write(conformal/'stages.json',stages);write(conformal/'snapshot-canonical.json',payload)
    digest=sha(conformal/'snapshot-canonical.json')
    write(conformal/'input-snapshot.json',dict(input_sha256=digest,snapshot=payload));write(conformal/'source.solidcraft.json',dict(test_only=True))
    record=dict(analysis_body='Clay',cad_body='CAD layer',step='body-0001.step',source_sha256=sha(conformal/'body-0001.step'),cad_length_unit='mm')
    capture=dict(schema='solidcraft-soil-capture/1',complete=True,input_sha256=digest,exporter_binary_sha256='e'*64,bodies=[record],
        files_sha256={p.name:sha(p) for p in conformal.iterdir()})
    write(conformal/'native-capture.json',capture)
    write(conformal/'native-capture-import.json',dict(passed=True,input_sha256=digest,source_capture_sha256=sha(conformal/'native-capture.json')))
    solver=analysis/'solver/soil-fem-3d.exe';solver.parent.mkdir(parents=True);solver.write_bytes(b'Not an executable: display-contract unit test')
    run['provenance']=dict(solver_binary_sha256=sha(solver))
    inputs={analysis/'ground/mesh.json':mesh,analysis/'ground/ground.json':ground,analysis/'stages.json':stages,analysis/'stage-results.json':run}
    for path,value in inputs.items():write(path,value)
    prepared=analysis/'ground'
    for path in conformal.iterdir():shutil.copyfile(path,prepared/('prepared-ground.json' if path.name=='ground.json' else path.name))
    write(prepared/'bodies.json',dict(regions=[dict(body='Clay',step=record['step'],source_sha256=record['source_sha256'])]))
    write(prepared/'mesh.provenance.json',dict(cad_manifest_sha256=sha(prepared/'bodies.json')))
    inventory={p.name:sha(p) for p in prepared.iterdir() if p.name!='prepared-ground.json'}
    write(prepared/'model-input-package.json',dict(schema='soil-fem-prepared-model-import/1',ground_sha256=sha(prepared/'ground.json'),
        original_inventory=dict(schema='soil-fem-prepared-model/1',files_sha256=inventory)))
    write(analysis/'accepted-stage-run.json',dict(schema='soil-fem-accepted-staged-run/1',inputs_and_results_sha256={str(p):sha(p) for p in inputs},solver_binary_sha256=sha(solver)))
    write(analysis/'stage-verification.json',dict(vtu_roundtrip='passed',raw_integration_point_roundtrip='passed',frames=[dict(independent_force_residual=0.)]*len(run['frames'])))
    write(analysis/'results/export-verification.json',dict(result_sha256=sha(analysis/'stage-results.json'),mesh_sha256=sha(analysis/'ground/mesh.json')))


class SolidCraftView(unittest.TestCase):
    def test_projected_export_audits_expanded_case_after_verifying_compact_cad_linkage(self):
        # Transport-contract regression. The numerical audit is a separate gate;
        # stop exactly at its boundary after checking which inputs it receives.
        with tempfile.TemporaryDirectory() as temp:
            # Hosted Windows TEMP may use an8.3 user alias. Match the exporter's
            # canonical path while preserving the directory-identity assertion.
            root=Path(temp).resolve()/'job';job_fixture(root)
            analysis=root/'analysis';conformal=root/'conformal'
            expanded=json.loads((analysis/'stages.json').read_text())
            declared=copy.deepcopy(expanded)
            declared['plate_surfaces']={'schema':'soil-fem-plate-surfaces/3','surfaces':[]}
            run=json.loads((analysis/'stage-results.json').read_text())
            run.update(schema='soil-fem-projected-mcc-stage-result/1',case=expanded)
            write(analysis/'stage-results.json',run)
            accepted=json.loads((analysis/'accepted-stage-run.json').read_text())
            accepted['inputs_and_results_sha256'][str(analysis/'stage-results.json')]=sha(analysis/'stage-results.json')
            write(analysis/'accepted-stage-run.json',accepted)
            capture=json.loads((conformal/'native-capture.json').read_text())
            ground=json.loads((conformal/'ground.json').read_text())
            def linkage(folder,prepared,mesh,case,captured):
                self.assertEqual(folder,conformal)
                self.assertEqual(case,expanded)
                self.assertEqual(captured,declared)
                return {'linked':True}
            def audited(mesh,actual_ground,case,result):
                self.assertEqual(actual_ground,ground)
                self.assertEqual(case,expanded)
                self.assertEqual(result['case'],expanded)
                raise RuntimeError('numerical audit boundary reached')
            with mock.patch('export_solidcraft_view.validate',return_value=(capture,ground,declared,{})), \
                 mock.patch('prepare_native_plates.verify',side_effect=linkage), \
                 mock.patch('verify_projected_stage_package.audit',side_effect=audited):
                with self.assertRaisesRegex(RuntimeError,'numerical audit boundary reached'):
                    export(root,root/'view')
            self.assertFalse((root/'view').exists())

    def test_weighted_cell_values_remain_separate_from_raw_points_and_removed_cells(self):
        mesh,run=fixture();view,points=scene_frame(mesh,run,1)
        self.assertEqual(view['active_element_ids'],[0]);self.assertEqual(len(view['boundary_faces']),4)
        self.assertEqual(view['cell_fields']['void_ratio'],[2.5])
        self.assertEqual(view['cell_fields']['effective_stress_compression_positive_kpa'][0],[25.,25.,25.,0.,0.,0.])
        self.assertEqual(points['fields']['void_ratio'],[1.,3.]);self.assertEqual(points['integration_point_ids'],[0,1])
        before,_=scene_frame(mesh,run,0);self.assertEqual(len(before['boundary_faces']),6)

    def test_reset_reference_seed_total_and_birth_increment_are_distinct(self):
        mesh,run=fixture();view,_=scene_frame(mesh,run,1)
        self.assertAlmostEqual(view['node_fields']['displacement_m'][0][2],-.1)
        self.assertAlmostEqual(view['node_fields']['stage_displacement_increment_m'][0][2],-.4)
        self.assertAlmostEqual(view['node_fields']['total_displacement_from_seed_m'][0][2],-.5)
        run['case']['stages'][0]['kind']='birth'
        run['transition_audits']=[dict(stage='cut',prepared_state=dict(displacement_m=[[0.,0.,-.45]]*5))]
        view,_=scene_frame(mesh,run,1)
        self.assertAlmostEqual(view['node_fields']['stage_displacement_increment_m'][0][2],-.05)

    def test_legacy_tet4_single_point_keeps_its_stored_value_and_true_centroid_volume(self):
        mesh,run=fixture();run['schema']='soil-fem-stage-result/1';run.pop('integration_points')
        for frame in run['frames']:
            frame['state']['material_states']=[frame['state']['material_states'][0],frame['state']['material_states'][2]]
        view,points=scene_frame(mesh,run,0)
        self.assertEqual(view['cell_fields']['void_ratio'],[1.,5.])
        self.assertEqual(points['fields']['void_ratio'],[1.,5.])
        np.testing.assert_allclose(points['weights_m3'],[1/6,1/6])
        np.testing.assert_allclose(points['positions_m'],[[.25,.25,.25],[.25,.25,-.25]])

    def test_quadratic_boundary_retains_midnodes_and_rejects_nonconforming_shared_faces(self):
        mesh,_=fixture();xyz=mesh['nodes'];edge_nodes={};cells=[]
        for tet in mesh['tetrahedra']:
            cell=tet.copy()
            for a,b in EDGES:
                key=tuple(sorted((tet[a],tet[b])))
                if key not in edge_nodes:
                    edge_nodes[key]=len(xyz);xyz.append(((np.array(xyz[key[0]])+np.array(xyz[key[1]]))/2).tolist())
                cell.append(edge_nodes[key])
            cells.append(cell)
        xyz=np.array(xyz);cells=np.array(cells)
        faces=boundary(cells,np.array([True,True]),xyz)
        self.assertEqual(len(faces),6);self.assertTrue(all(len(f['nodes'])==6 for f in faces))
        self.assertNotIn({0,1,2},[set(f['nodes'][:3]) for f in faces])
        after=boundary(cells,np.array([True,False]),xyz)
        self.assertIn({0,1,2},[set(f['nodes'][:3]) for f in after])
        for face in faces:
            corners=xyz[face['nodes'][:3]];normal=np.cross(corners[1]-corners[0],corners[2]-corners[0])
            self.assertGreater(np.dot(normal,corners.mean(axis=0)-xyz[cells[face['element_index'],:4]].mean(axis=0)),0.)
        cells[1,4]=len(xyz)
        with self.assertRaisesRegex(ValueError,'nonconforming'):boundary(cells,np.array([True,True]),np.vstack([xyz,xyz[edge_nodes[(0,2)]]]))

    def test_model_specific_missing_value_is_null_and_bad_required_values_are_rejected(self):
        mesh,run=fixture();run['frames'][0]['state']['material_states'][2].pop('preconsolidation_pressure_kpa')
        view,points=scene_frame(mesh,run,0)
        self.assertIsNone(view['cell_fields']['preconsolidation_pressure_kpa'][1]);self.assertIsNone(points['fields']['preconsolidation_pressure_kpa'][2])
        run['frames'][0]['state']['material_states'][0]['void_ratio']=float('nan')
        with self.assertRaisesRegex(ValueError,'non-finite .*material field'):scene_frame(mesh,run,0)

    def test_complete_packet_is_version_bound_relocatable_and_never_overwrites_or_accepts_tampering(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);job=root/'job';job_fixture(job)
            view=export(job,root/'view')
            self.assertEqual(view,json.loads((root/'view/view-manifest.json').read_text(encoding='utf-8')))
            self.assertTrue(view['complete']);self.assertEqual(view['input_sha256'],sha(job/'conformal/snapshot-canonical.json'))
            self.assertIn('not smoothed',view['semantics']['cell_fields'])
            for name,digest in view['files_sha256'].items():self.assertEqual(sha(root/'view'/name),digest)
            with self.assertRaises(FileExistsError):export(job,root/'view')
            # Receipt paths identify hashes, never files to read outside this job.
            accepted_path=job/'analysis/accepted-stage-run.json';a=json.loads(accepted_path.read_text())
            a['inputs_and_results_sha256']={'Z:/old-unavailable-job/'+Path(p).name:h for p,h in a['inputs_and_results_sha256'].items()};write(accepted_path,a)
            export(job,root/'relocated-view')
            p=job/'analysis/ground/mesh.json';p.write_bytes(p.read_bytes()+b' ')
            with self.assertRaisesRegex(ValueError,'hash mismatch'):export(job,root/'rejected')
            self.assertFalse((root/'rejected').exists())

    def test_resealed_solver_receipt_cannot_detach_result_from_its_prepared_cad_geometry(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);job=root/'job';job_fixture(job)
            mesh=job/'analysis/ground/mesh.json';m=json.loads(mesh.read_text());m['nodes'][3][2]=2.;write(mesh,m)
            path=job/'analysis/accepted-stage-run.json';a=json.loads(path.read_text());a['inputs_and_results_sha256'][str(mesh)]=sha(mesh);write(path,a)
            with self.assertRaisesRegex(ValueError,'prepared geometry file hash mismatch'):export(job,root/'rejected')
            self.assertFalse((root/'rejected').exists())


class SourceBudget(unittest.TestCase):
 def test_only_explicit_result_reader_accepts_a_larger_source(self):
  with tempfile.TemporaryDirectory() as folder:
   path=Path(folder)/'result.json';value=dict(material_history=list(range(100)));path.write_text(json.dumps(value),encoding='utf-8')
   # Bind small test budgets without constructing a giant test artifact.
   with self.assertRaises(ValueError):view.read(path,maximum=32)
   self.assertEqual(view.read_result(path),value)
 def test_source_limit_is_independent_of_packet_limits(self):
  old=view.MAX_SOURCE_RESULT
  try:
   view.MAX_SOURCE_RESULT=16
   with tempfile.TemporaryDirectory() as folder:
    path=Path(folder)/'source.json';path.write_text(json.dumps(dict(history=[1]*100)),encoding='utf-8')
    with self.assertRaisesRegex(ValueError,'source limit'):view.read_result(path)
    self.assertEqual(view.read(path),dict(history=[1]*100))
  finally:view.MAX_SOURCE_RESULT=old
  self.assertEqual(view.MAX_FILE,256<<20);self.assertEqual(view.MAX_TOTAL,512<<20)
 def test_nonfinite_source_numbers_remain_rejected(self):
  with tempfile.TemporaryDirectory() as folder:
   path=Path(folder)/'bad.json'
   for raw in ['{"a":NaN}','{"a":1e999}']:
    path.write_text(raw,encoding='utf-8')
    with self.assertRaises(ValueError):view.read_result(path)


if __name__=='__main__':unittest.main()
