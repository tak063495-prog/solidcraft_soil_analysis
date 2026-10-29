"""Material IDs, RESET meaning, raw GP fields and refused unsafe export inputs."""
import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

import meshio
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
import export_mc_epoch as module


def fixture():
    mesh = dict(nodes=[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,0,-1]],
        tetrahedra=[[0,1,2,3],[0,2,1,4]], material_ids=[99,7], region_ids=[3,5])
    def parameter(c):
        return dict(young_kpa=1000.,poisson=.2,cohesion_kpa=c,friction_deg=20.,
                    dilation_deg=0.,tension_cutoff_kpa=None,flow='non_associated')
    def state(stress):
        return dict(effective_stress_compression_positive_kpa=[stress]*3+[0.]*3,
            plastic_strain_compression_positive=[.1,0.,0.,0.,0.,0.],void_ratio=.8,
            equivalent_plastic_strain=.2,plastic_work_kpa=.3)
    original = dict(active_elements=[True,False], element_material_ids=[7,99],
        displacement_m=[[0.,0.,-.01]]*5, display_reference_displacement_m=[[0.,0.,-.008]]*5,
        material_states=[state(2.),state(6.),state(800.)])
    # Unequal quadrature weights make the exact average 5, not the arithmetic 4.
    layout = [dict(element_index=0,position_m=[.25]*3,weight_m3=w,
                   shape_functions=[.25]*4) for w in (1.,3.)]
    layout.append(dict(element_index=1,position_m=[.25,.25,-.25],weight_m3=2.,shape_functions=[.25]*4))
    origin = dict(last_committed=original,integration_points=layout,
        ground=dict(materials=[dict(id=99,mohr_coulomb=parameter(19.)),dict(id=7,mohr_coulomb=parameter(11.))]))
    committed = copy.deepcopy(original)
    committed['displacement_m'] = [[0.,0.,-.013]]*5
    checkpoint = dict(payload=dict(material_checkpoint=dict(payload=dict(
        committed=committed,parameters=[parameter(9.5),parameter(5.5)]))))
    return mesh,origin,checkpoint


class EpochFields(unittest.TestCase):
    def test_material_ids_and_removed_points(self):
        mesh,origin,checkpoint = fixture()
        xyz,kind,cells,nodal,cell,gp_xyz,gp = module.fields(mesh,origin,checkpoint,True)
        self.assertEqual(kind,'tetra')
        self.assertEqual(cells.tolist(),[[0,1,2,3]])
        self.assertEqual(cell['material_id'][0].tolist(),[7])
        self.assertEqual(gp['integration_point_id'].tolist(),[0,1])
        self.assertEqual(gp['cohesion_kpa'].tolist(),[5.5,5.5])
        self.assertEqual(cell['volume_average_cohesion_kpa'][0].tolist(),[5.5])
        self.assertTrue(np.isnan(gp['tension_cutoff_kpa']).all())
        self.assertEqual(gp['has_tension_cutoff'].tolist(),[0,0])

    def test_display_reset_and_increment_are_distinct(self):
        values = module.fields(*fixture(),True)
        nodal = values[3]
        np.testing.assert_allclose(nodal['total_displacement_from_seed_m'][:,2],-.013)
        np.testing.assert_allclose(nodal['displacement_m'][:,2],-.005)
        np.testing.assert_allclose(nodal['epoch_displacement_increment_m'][:,2],-.003)
        np.testing.assert_allclose(values[6]['epoch_displacement_increment_m'][:,2],-.003)

    def test_raw_stress_is_not_smoothed_and_average_uses_volume(self):
        values = module.fields(*fixture(),True)
        self.assertEqual(values[6]['effective_stress_compression_positive_kpa'][:,0].tolist(),[2.,6.])
        self.assertEqual(values[4]['volume_average_effective_stress_compression_positive_kpa'][0][0,0],5.)
        self.assertNotIn('effective_stress_compression_positive_kpa',values[3])

    def test_origin_uses_origin_parameters_and_zero_increment(self):
        values = module.fields(*fixture(),False)
        self.assertEqual(values[6]['cohesion_kpa'].tolist(),[11.,11.])
        self.assertTrue((values[3]['epoch_displacement_increment_m']==0.).all())

    def test_full_packet_roundtrip_and_partial_progress(self):
        mesh,origin,checkpoint = fixture()
        info = dict(complete=False,progress=.375,stop_reason=None,generation=4)
        manifest = dict(job_id='test',binding=dict(files_sha256={'x':'y'}))
        source = b'{"dedicated_snapshot":true}'
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);job=root/'job';job.mkdir();out=root/'out'
            with patch.object(module,'accepted',return_value=(mesh,origin,checkpoint,info,source,manifest)):
                packet=module.export(job,out)
                self.assertFalse(packet['complete'])
                self.assertEqual(packet['stability_classification'],'not_evaluated')
                self.assertEqual(packet['parameter_progress'],.375)
                self.assertFalse(packet['native_solidcraft_reader_supported'])
                self.assertEqual((out/'accepted-snapshot.json').read_bytes(),source)
                self.assertEqual(packet['files_sha256']['accepted-snapshot.json'],hashlib.sha256(source).hexdigest())
                self.assertEqual(meshio.read(out/'current.vtu').cell_data['material_id'][0].tolist(),[7])
                self.assertEqual(json.loads((out/'export-verification.json').read_text('utf-8'))['schema'],
                                 'soil-fem-mc-epoch-export/1')
                with self.assertRaises(ValueError):module.export(job,out)

    def test_refuses_output_in_source_without_reading_or_writing(self):
        with tempfile.TemporaryDirectory() as folder:
            job=Path(folder)/'job';job.mkdir()
            with patch.object(module,'accepted') as reader:
                for output in (job,job/'child',job.parent):
                    with self.assertRaises(ValueError):module.export(job,output)
                reader.assert_not_called()
                self.assertEqual(list(job.iterdir()),[])

    def test_decode_rejects_nonfinite(self):
        for value in ('NaN','Infinity','-Infinity','1e999'):
            with self.assertRaises(ValueError):module.decode('{"x":'+value+'}')

    def test_running_writer_has_bounded_retries(self):
        with tempfile.TemporaryDirectory() as folder:
            job=Path(folder);(job/'inputs').mkdir();(job/'solver').mkdir()
            data=b'{}';files={n:hashlib.sha256(data).hexdigest() for n in module.FILES}
            for n in files:(job/n).write_bytes(data)
            (job/'job.json').write_text(json.dumps(dict(schema='soil-fem-mc-epoch-job/1',job_id='test',
                binding=dict(parent=None,files_sha256=files))),'utf-8')
            def pointer(generation):
                name=f'recovery-checkpoint-{generation%2}.json';(job/name).write_bytes(data)
                (job/'recovery.json').write_text(json.dumps(dict(schema='soil-fem-mc-epoch-recovery/1',
                    generation=generation,job_id='test',checkpoint=name,
                    checkpoint_sha256=hashlib.sha256(data).hexdigest())),'utf-8')
            pointer(1)
            def advance(*args,**kwargs):
                generation=json.loads((job/'recovery.json').read_text('utf-8'))['generation']
                pointer(generation+1)
                return SimpleNamespace(returncode=0,stdout=b'{}',stderr=b'')
            with patch.object(module.subprocess,'run',side_effect=advance) as process:
                with self.assertRaisesRegex(ValueError,'all three reads'):module.accepted(job)
                self.assertEqual(process.call_count,3)

    def test_tet10_keeps_midside_nodes_and_shape_interpolation(self):
        mesh,origin,checkpoint=fixture()
        corners=np.asarray(mesh['nodes'][:4],dtype=float)
        nodes=list(corners)+[(corners[a]+corners[b])/2 for a,b in ((0,1),(1,2),(2,0),(0,3),(1,3),(2,3))]
        mesh.update(nodes=np.asarray(nodes).tolist(),tetrahedra10=[list(range(10))],material_ids=[7],region_ids=[3])
        for state in (origin['last_committed'],checkpoint['payload']['material_checkpoint']['payload']['committed']):
            state.update(active_elements=[True],element_material_ids=[7],material_states=state['material_states'][:2])
            state['displacement_m']=[state['displacement_m'][0]]*10
            state['display_reference_displacement_m']=[state['display_reference_displacement_m'][0]]*10
        origin['integration_points']=origin['integration_points'][:2]
        for row in origin['integration_points']:row['shape_functions']=[-.125]*4+[.25]*6
        values=module.fields(mesh,origin,checkpoint,True)
        self.assertEqual(values[1],'tetra10');self.assertEqual(len(values[0]),10)
        self.assertEqual(values[2].tolist(),[list(range(10))])
        np.testing.assert_allclose(values[6]['epoch_displacement_increment_m'][:,2],-.003)


if __name__=='__main__':unittest.main()
