import copy,json,sys,tempfile,unittest
from unittest.mock import patch
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
import source_point_descriptor as points

class SourcePointDescriptorTests(unittest.TestCase):
 def fixture(self):
  return dict(schema='soil-fem-projected-mcc-stage-result/2',complete=True,failure_reason=None,
    integration_points=[{'element_index':0},{'element_index':1},{'element_index':0}],
    frames=[dict(name='excavate',state=dict(active_elements=[True,False],material_states=[{}, {}, {}]))])
 def test_exact_source_count_active_mask_frame_units_and_components(self):
  run=self.fixture();descriptor=points.make(run,0,'a'*64)
  self.assertEqual(descriptor,dict(schema='soil-fem-view-source-points/1',
    source_result_schema='soil-fem-projected-mcc-stage-result/2',source_result_sha256='a'*64,
    frame_index=0,frame_name='excavate',integration_point_count=3,active_integration_point_count=2,
    length_unit='m',stress_unit='kPa',stress_sign='compression_positive',stress_components=['xx','yy','zz','xy','yz','xz']))
  self.assertNotIn('fields',descriptor);self.assertNotIn('positions_m',descriptor)
  points.validate(descriptor,run,0,'a'*64)
 def test_each_binding_and_semantic_tamper_is_rejected(self):
  run=self.fixture();base=points.make(run,0,'a'*64)
  for key,value in [('source_result_sha256','b'*64),('frame_index',True),('frame_name','preload'),
    ('integration_point_count',2),('active_integration_point_count',3),('length_unit','mm'),
    ('stress_unit','Pa'),('stress_sign','tension_positive'),('stress_components',['zz','yy','xx','xy','yz','xz']),
    ('schema','soil-fem-view-points/1'),('unexpected','value')]:
   descriptor=copy.deepcopy(base);descriptor[key]=value
   with self.subTest(key=key),self.assertRaises(ValueError):points.validate(descriptor,run,0,'a'*64)
 def test_incomplete_wrong_family_invalid_hash_and_frame_are_rejected(self):
  for key,value in [('complete',False),('failure_reason','not converged'),('schema','soil-fem-stage-result/2')]:
   run=self.fixture();run[key]=value
   with self.assertRaises(ValueError):points.make(run,0,'a'*64)
  for index in [-1,1,True]:
   with self.assertRaises(ValueError):points.make(self.fixture(),index,'a'*64)
  for digest in ['a'*63,'A'*64,'../result.json']:
   with self.assertRaises(ValueError):points.make(self.fixture(),0,digest)
 def test_invalid_inventory_owner_and_empty_domain_are_rejected(self):
  run=self.fixture();run['frames'][0]['state']['material_states'].pop()
  with self.assertRaises(ValueError):points.make(run,0,'a'*64)
  for owner in [-1,2,True]:
   run=self.fixture();run['integration_points'][0]['element_index']=owner
   with self.assertRaises(ValueError):points.make(run,0,'a'*64)
  run=self.fixture();run['frames'][0]['state']['active_elements']=[False,False]
  with self.assertRaises(ValueError):points.make(run,0,'a'*64)

 def test_view_verifier_requires_explicit_version4_and_complete_source_semantics(self):
  from verify_solidcraft_view import modes
  manifest=dict(schema=points.VIEW_SCHEMA,complete=True,semantics=dict(source_result_schema='soil-fem-projected-mcc-stage-result/2',
      raw_point_storage='authoritative_source',source_point_descriptor_schema=points.SCHEMA))
  self.assertEqual(modes(manifest),(True,True,True))
  for key,value in [('raw_point_storage','averaged_cells'),('source_point_descriptor_schema','other'),('source_result_schema','soil-fem-stage-result/6')]:
   wrong=copy.deepcopy(manifest);wrong['semantics'][key]=value
   with self.subTest(key=key),self.assertRaises(ValueError):modes(wrong)
  for schema in ['soil-fem-solidcraft-view/1','soil-fem-solidcraft-view/2','soil-fem-solidcraft-view/3']:
   wrong=copy.deepcopy(manifest);wrong['schema']=schema
   with self.subTest(schema=schema),self.assertRaises(ValueError):modes(wrong)

 def test_view_verifier_rejects_resealed_descriptor_tamper_before_reading_values(self):
  from verify_solidcraft_view import source_point_fields
  run=self.fixture();descriptor=points.make(run,0,'a'*64)
  for key,value in [('source_result_sha256','b'*64),('frame_index',True),('active_integration_point_count',3),('stress_unit','Pa')]:
   wrong=copy.deepcopy(descriptor);wrong[key]=value
   with self.subTest(key=key),patch('export_solidcraft_view._scene_frame') as extract:
    with self.assertRaises(ValueError):source_point_fields({},run,0,wrong,'a'*64)
    extract.assert_not_called()
 def test_json_duplicate_key_cannot_be_hidden_by_last_value_wins(self):
  descriptor=points.make(self.fixture(),0,'a'*64)
  encoded=json.dumps(descriptor).replace('"stress_unit": "kPa"','"stress_unit": "Pa", "stress_unit": "kPa"')
  points.validate(json.loads(encoded),self.fixture(),0,'a'*64)
  with tempfile.TemporaryDirectory() as tmp:
   path=Path(tmp)/'points.json';path.write_text(encoded,'utf-8')
   with self.assertRaisesRegex(ValueError,'duplicate'):points.load(path)

 def test_json_size_nonfinite_and_nonobject_are_rejected(self):
  with tempfile.TemporaryDirectory() as tmp:
   path=Path(tmp)/'points.json'
   for text in [' '*8193,'{"frame_index":NaN}','{"frame_index":1e999}','[]']:
    path.write_text(text,'utf-8')
    with self.subTest(text=text[:30]),self.assertRaises(ValueError):points.load(path)

 def manifest_fixture(self):
  return dict(schema=points.VIEW_SCHEMA,complete=True,geometry='geometry.json',
   frames=[dict(index=0,name='seed',frame='frame-000.json',integration_points='points-000.json')],
   files_sha256={name:'a'*64 for name in ['geometry.json','frame-000.json','points-000.json']},
   job_sources_sha256={name:'b'*64 for name in ['analysis/ground/mesh.json','analysis/stage-results.json','analysis/accepted-stage-run.json',
       'analysis/solver/soil-fem-3d.exe','analysis/ground/snapshot-canonical.json','analysis/ground/native-capture.json']},
   semantics=dict(length_unit='m',stress_unit='kPa',stress_sign='compression positive',stress_components=points.COMPONENTS.copy()))

 def test_manifest_resealed_path_and_frame_inventory_tampering_is_rejected(self):
  from verify_solidcraft_view import reference_inventory
  original=self.manifest_fixture();reference_inventory(original)
  for key,value in [('index',True),('index',1),('integration_points','../points-000.json'),('frame','frame-001.json')]:
   wrong=copy.deepcopy(original);wrong['frames'][0][key]=value
   with self.subTest(key=key),self.assertRaises(ValueError):reference_inventory(wrong)
  for name in ['../outside.json','analysis/../outside.json','C:/outside.json','analysis\\outside.json']:
   wrong=copy.deepcopy(original);wrong['job_sources_sha256'][name]='a'*64
   with self.subTest(name=name),self.assertRaises(ValueError):reference_inventory(wrong)
  wrong=copy.deepcopy(original);wrong['files_sha256']['unlisted.json']='a'*64
  with self.assertRaises(ValueError):reference_inventory(wrong)

 def test_manifest_units_and_missing_source_are_rejected(self):
  from verify_solidcraft_view import reference_inventory
  original=self.manifest_fixture()
  for key,value in [('length_unit','mm'),('stress_unit','Pa'),('stress_sign','tension positive'),('stress_components',list(reversed(points.COMPONENTS)))]:
   wrong=copy.deepcopy(original);wrong['semantics'][key]=value
   with self.subTest(key=key),self.assertRaises(ValueError):reference_inventory(wrong)
  wrong=copy.deepcopy(original);del wrong['job_sources_sha256']['analysis/stage-results.json']
  with self.assertRaises(ValueError):reference_inventory(wrong)

if __name__=='__main__':unittest.main()
