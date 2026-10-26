"""Exact raw structural display fields, including inactive archived histories."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from test_mcc_factor_reference import specimen
from test_mcc_structural_reference import configuration
from beam_reference import System
from solidcraft_structure_fields import fields


def clean(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,dict):return {k:clean(v) for k,v in value.items()}
    if isinstance(value,list):return [clean(v) for v in value]
    return value


def fixture(active=True):
    mesh,_,_=specimen();cfg=configuration(mesh);model=System(mesh,cfg);history=model.initial_state()
    if not active:
        history['active_members']=[False];history['active_plates']=[False]
        history['node_displacement_m'][0]=[.02,.03,.04]
        history['node_rotation_rad'][0]=[.002,.003,.004]
    u=np.zeros((len(mesh['nodes']),3));response=model.solve_interface(history,u,.4)
    report=dict(node_positions_m=model.positions,node_displacement_m=response['state']['node_displacement_m'],
        node_rotation_rad=response['state']['node_rotation_rad'],node_active=response['node_active'],
        members=[dict(id=k,**v) for k,v in response['members'].items()],
        plates=[dict(id=k,**v) for k,v in response.get('plates',{}).items()])
    for name,key in (('external_force_kn','external_force'),('external_moment_origin_kn_m','external_moment'),
            ('support_force_kn','support_force'),('support_moment_origin_kn_m','support_moment')):report[name]=response[key]
    run=dict(case=dict(beam_system=cfg),frames=[dict(state=dict(beams=response['state'],displacement_m=u.tolist(),
        display_reference_displacement_m=[[.1,.2,.3]]*len(u),boundary_factor=.4),beam_results=clean(report))])
    return mesh,run


def scene_fixture():
    mesh,run=fixture();mesh.update(region_ids=[1],regions=[dict(id=1,body='Clay')])
    run.update(schema='soil-fem-stage-result/14',complete=True,failure_reason=None,
        integration_points=[dict(element_index=0,position_m=[.25,.25,.25],weight_m3=1/6)])
    frame=run['frames'][0];frame.update(name='transport',monitors={},relative_force_residual=0.,maximum_normalized_yield=-.1)
    frame['state'].update(active_elements=[True],material_states=[dict(effective_stress_compression_positive_kpa=[100.]*3+[0.]*3,
        preconsolidation_pressure_kpa=200.,void_ratio=1.,plastic_strain_compression_positive=[0.]*6)])
    run['last_committed']=copy.deepcopy(frame['state'])
    return mesh,run


class StructuralFields(unittest.TestCase):
    def test_legacy_beam_configuration_retains_the_same_versioned_transport(self):
        mesh,run=fixture();run['case']['beam_system']['schema']='soil-fem-beams/1'
        self.assertEqual(fields(mesh,run,0)['schema'],'soil-fem-view-structures/1')

    def test_local_axes_section_samples_and_installation_references_are_exact_and_unsmoothed(self):
        mesh,run=fixture();before=copy.deepcopy((mesh,run));packet=fields(mesh,run,0)
        self.assertEqual((mesh,run),before)
        self.assertEqual(packet['members'][0]['result'],run['frames'][0]['beam_results']['members'][0])
        self.assertEqual(packet['plates'][0]['result'],run['frames'][0]['beam_results']['plates'][0])
        self.assertEqual(packet['nodes']['displacement_m'],run['frames'][0]['state']['beams']['node_displacement_m'])
        self.assertEqual(packet['nodes']['soil_node_ids'],[0,1,2,None])
        json.dumps(packet,allow_nan=False)

    def test_soil_display_reset_does_not_change_member_coordinates_or_displacements(self):
        mesh,run=fixture();expected=fields(mesh,run,0)
        run['frames'][0]['state']['display_reference_displacement_m']=[[12.,23.,34.]]*len(mesh['nodes'])
        self.assertEqual(fields(mesh,run,0),expected)

    def test_inactive_results_are_null_and_archived_kinematics_remain_available(self):
        mesh,run=fixture(False);packet=fields(mesh,run,0)
        self.assertEqual(packet['nodes']['active'],[False]*4)
        self.assertEqual(packet['nodes']['displacement_m'][0],[.02,.03,.04])
        self.assertEqual(packet['nodes']['rotation_rad'][0],[.002,.003,.004])
        self.assertFalse(packet['members'][0]['active']);self.assertIsNone(packet['members'][0]['result'])
        self.assertFalse(packet['plates'][0]['active']);self.assertIsNone(packet['plates'][0]['result'])

    def test_missing_duplicate_and_modified_result_fields_are_rejected(self):
        for mutation in ('missing','duplicate','moment','rotation','position'):
            mesh,run=fixture();report=run['frames'][0]['beam_results']
            if mutation=='missing':report['plates']=[]
            elif mutation=='duplicate':report['plates']*=2
            elif mutation=='moment':report['plates'][0]['samples'][0]['bending_moment_kn'][0]+=.01
            elif mutation=='rotation':report['node_rotation_rad'][0][0]+=.01
            else:report['node_positions_m'][0][0]+=.01
            with self.subTest(mutation=mutation),self.assertRaises(ValueError):fields(mesh,run,0)

    def test_nonfinite_inactive_reference_is_rejected_and_absent_configuration_is_explicit(self):
        mesh,run=fixture(False);run['frames'][0]['state']['beams']['plate_reference_global_kinematics'][0][0]=float('nan')
        with self.assertRaises(ValueError):fields(mesh,run,0)
        run=dict(case={},frames=[dict(state={})]);self.assertIsNone(fields(mesh,run,0))
        run['frames'][0]['beam_results']={}
        with self.assertRaises(ValueError):fields(mesh,run,0)

    def test_separately_read_vtu_preserves_active_and_inactive_geometry_and_all_section_samples(self):
        from beam_io import export as beam_export
        from plate_io import export as plate_export
        from verify_structural_view_fields import verify
        for active in (True,False):
            mesh,run=fixture(active);run['frames'][0]['name']='transport'
            with tempfile.TemporaryDirectory() as tmp:
                folder=Path(tmp);beam_export(run,folder);plate_export(run,folder)
                packet=fields(mesh,run,0)
                self.assertEqual(verify(packet,folder,0),dict(members=int(active),plates=int(active)))
                if active:
                    packet['plates'][0]['result']['samples'][0]['bending_moment_kn'][0]+=.1
                    with self.assertRaisesRegex(ValueError,'raw plate'):verify(packet,folder,0)

    def test_structural_scene_uses_version_three_and_soil_raw_points_remain_separate(self):
        from export_solidcraft_view import scene_frame
        mesh,run=scene_fixture()
        view,points=scene_frame(mesh,run,0)
        self.assertEqual(view['schema'],'soil-fem-view-frame/3');self.assertEqual(points['schema'],'soil-fem-view-integration-points/3')
        self.assertEqual(view['structures'],fields(mesh,run,0))
        self.assertNotIn('bending_moment_kn',points['fields'])

    def test_complete_structural_transport_packet_has_explicit_version_and_sources(self):
        from test_solidcraft_view import job_fixture
        from export_solidcraft_view import export
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);job=root/'job';job_fixture(job,scene_fixture())
            manifest=export(job,root/'view')
            self.assertEqual(manifest['schema'],'soil-fem-solidcraft-view/3')
            self.assertEqual(manifest['semantics']['source_result_schema'],'soil-fem-stage-result/14')
            self.assertEqual(json.loads((root/'view/frame-000.json').read_text())['structures']['schema'],'soil-fem-view-structures/1')
            self.assertNotIn('beam/anchor/contact result layers',manifest['omissions'])

    def test_resealed_structural_moment_cannot_publish_a_complete_packet(self):
        from test_solidcraft_view import job_fixture,write
        from export_solidcraft_view import export
        from verify_acceptance_packet import sha
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);job=root/'job';job_fixture(job,scene_fixture());analysis=job/'analysis'
            path=analysis/'stage-results.json';run=json.loads(path.read_text());run['frames'][0]['beam_results']['plates'][0]['samples'][0]['bending_moment_kn'][0]+=.1;write(path,run)
            receipt=analysis/'accepted-stage-run.json';data=json.loads(receipt.read_text());data['inputs_and_results_sha256'][str(path)]=sha(path);write(receipt,data)
            exported=analysis/'results/export-verification.json';data=json.loads(exported.read_text());data['result_sha256']=sha(path);write(exported,data)
            with self.assertRaisesRegex(ValueError,'independent member equilibrium'):export(job,root/'rejected')
            self.assertFalse((root/'rejected/view-manifest.json').exists())


if __name__=='__main__':unittest.main()
