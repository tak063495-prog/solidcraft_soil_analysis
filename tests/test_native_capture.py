"""Importer contract gates; actual native-SolidCraft integration is recorded separately."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tools'))
from import_solidcraft_capture import import_capture,validate
from verify_acceptance_packet import sha


def write(path,value):path.write_text(json.dumps(value,ensure_ascii=False,sort_keys=True,separators=(',',':')),encoding='utf-8')


class NativeCapture(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='soil-fem-native-capture-');self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.capture=self.root/'capture';self.capture.mkdir()
        import gmsh
        gmsh.initialize()
        try:
            gmsh.option.setNumber('General.Terminal',0)
            gmsh.model.occ.addBox(0,0,0,10000,10000,4000);gmsh.model.occ.synchronize();gmsh.write(str(self.capture/'body-0001.step'))
        finally:gmsh.finalize()
        ground=json.loads((ROOT/'examples/ground-local.json').read_text());ground['layers']=ground['layers'][:1];ground['layers'][0]['body']='Clay';ground['excavation_bodies']=[]
        stages=dict(schema='soil-fem-stages/1',stages=[])
        payload=dict(schema='solidcraft-soil-snapshot/1',ground=ground,stages=stages,cad_sha256='c'*64,cad_length_unit='mm',analysis_units='m-kN-kPa',
            bindings=[dict(analysis_body='Clay',cad_body='CAD renamed earth',creator_feature=1,component=0)])
        write(self.capture/'ground.json',ground);write(self.capture/'stages.json',stages);write(self.capture/'snapshot-canonical.json',payload)
        fingerprint=sha(self.capture/'snapshot-canonical.json')
        write(self.capture/'input-snapshot.json',dict(input_sha256=fingerprint,snapshot=payload))
        write(self.capture/'source.solidcraft.json',{'format':'future/source-artifact','commands':'Never execute source designs in the importer'})
        record=dict(analysis_body='Clay',cad_body='CAD renamed earth',step='body-0001.step',source_sha256=sha(self.capture/'body-0001.step'),cad_length_unit='mm')
        self.receipt=dict(schema='solidcraft-soil-capture/1',complete=True,input_sha256=fingerprint,exporter_binary_sha256='e'*64,exporter_version='test',bodies=[record],
            files_sha256={p.name:sha(p) for p in self.capture.iterdir()})
        write(self.capture/'capture.json',self.receipt)

    def update(self):write(self.capture/'capture.json',self.receipt)

    def test_logical_labels_and_mm_to_m_volume_are_independent_of_cad_display_names(self):
        manifest,link=import_capture(self.capture,self.root/'imported');r,=manifest['regions']
        self.assertEqual(r['body'],'Clay');self.assertEqual(r['material_id'],1)
        self.assertAlmostEqual(r['volume_m3'],400.,places=8)
        self.assertLess(abs(r['bounds_m'][1][2]-4.),1e-8)
        self.assertEqual(link['cad_body_to_analysis_body'],{'CAD renamed earth':'Clay'})
        self.assertFalse(link['equilibrium_verified'])
        self.assertEqual(sha(self.capture/'source.solidcraft.json'),sha(self.root/'imported/source.solidcraft.json'))

    def test_resealed_ground_changes_and_bool_number_coercion_are_rejected(self):
        ground=json.loads((self.capture/'ground.json').read_text());ground['materials'][0]['id']=True
        write(self.capture/'ground.json',ground);self.receipt['files_sha256']['ground.json']=sha(self.capture/'ground.json');self.update()
        with self.assertRaisesRegex(ValueError,'differ'):validate(self.capture)
        self.assertFalse((self.root/'imported').exists())

    def test_nonfinite_json_is_rejected_before_file_import(self):
        self.receipt['unused']=float('nan');self.update()
        with self.assertRaisesRegex(ValueError,'non-finite JSON'):validate(self.capture)

    def test_overflowing_json_exponent_is_rejected_before_file_import(self):
        path=self.capture/'capture.json';data=path.read_text(encoding='utf-8')
        path.write_text(data[:-1]+',"overflow":1e309}',encoding='utf-8')
        with self.assertRaisesRegex(ValueError,'non-finite JSON'):validate(self.capture)

    def test_step_tampering_and_escaping_inventory_paths_are_rejected(self):
        old=(self.capture/'body-0001.step').read_bytes();(self.capture/'body-0001.step').write_bytes(old+b'changed')
        with self.assertRaisesRegex(ValueError,'hash mismatch'):validate(self.capture)
        (self.capture/'body-0001.step').write_bytes(old)
        self.receipt['files_sha256']['../outside.step']=self.receipt['files_sha256'].pop('body-0001.step')
        self.receipt['bodies'][0]['step']='../outside.step';self.update()
        with self.assertRaisesRegex(ValueError,'filename'):validate(self.capture)

    def test_existing_output_and_missing_completion_marker_are_preserved(self):
        output=self.root/'output';output.mkdir();(output/'keep.txt').write_text('keep')
        with self.assertRaises(FileExistsError):import_capture(self.capture,output)
        self.assertEqual((output/'keep.txt').read_text(),'keep')
        self.receipt['complete']=False;self.update()
        with self.assertRaisesRegex(ValueError,'complete'):validate(self.capture)


if __name__=='__main__':unittest.main()
