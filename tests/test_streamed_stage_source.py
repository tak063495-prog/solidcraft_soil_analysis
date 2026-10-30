import io,json,sys,tempfile,unittest
from collections.abc import Mapping,Sequence
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
import streamed_stage_source as stream

def materialize(value):
 if isinstance(value,Mapping):return {k:materialize(value[k]) for k in value}
 if isinstance(value,Sequence) and not isinstance(value,(str,bytes)):return [materialize(v) for v in value]
 return value

class StreamSourceTests(unittest.TestCase):
 def setUp(self):
  self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)
  self.source=self.root/'run.json';self.cache=self.root/'cache'
 def fixture(self):
  return dict(schema='test',complete=True,case={'stages':[{'name':'a'}]},
   frames=[dict(name='a',state={'material_states':[{'void_ratio':.7}],'displacement_m':[[0.,1.,2.]]}),
           dict(name='b',state={'material_states':[{'void_ratio':.6}],'displacement_m':[[0.,2.,3.]]})],
   geostatic=dict(seed={'material_states':[{'void_ratio':.8}]},last_committed={'material_states':[{'void_ratio':.7}]},progress=1.,complete=True,attempts=[]),
   transition_audits=[{'stage':'b','prepared_state':{'material_states':[{'void_ratio':.65}]}}],
   integration_points=[dict(element_index=0,position_m=[1.,2.,3.],weight_m3=1.)],last_committed={'material_states':[{'void_ratio':.6}]})
 def write(self,value):self.source.write_text(json.dumps(value,ensure_ascii=False,allow_nan=False),'utf-8')
 def test_all_fields_frames_seed_last_and_audits_equal_with_tiny_LRU(self):
  expected=self.fixture();self.write(expected);actual=stream.index_source(self.source,self.cache,cache_payload=128)
  self.assertEqual(materialize(actual),expected);self.assertEqual(actual['frames'][-1]['name'],'b')
  self.assertEqual(materialize(actual['frames'][:1]),expected['frames'][:1]);actual.verify_source()
  self.assertLessEqual(actual._store.retained,128)
 def test_null_geostatic_and_empty_transition_array_preserved(self):
  expected=self.fixture();expected['geostatic']=None;expected['transition_audits']=[];self.write(expected)
  self.assertEqual(materialize(stream.index_source(self.source,self.cache)),expected)
 def test_duplicate_keys_rejected_in_nested_and_root_and_no_complete_index(self):
  for index,raw in enumerate(['{"a":1,"a":2}','{"frames":[{"state":{"v":1,"v":2}}]}']):
   self.source.write_text(raw,'utf-8');cache=self.root/str(index)
   with self.assertRaises(ValueError):stream.index_source(self.source,cache)
   self.assertFalse((cache/'index.json').exists())
 def test_source_byte_component_and_retained_budgets_are_independent(self):
  self.write(self.fixture())
  for name,kwargs in [('source',dict(source_limit=2)),('component',dict(component_limit=2)),('retained',dict(payload_limit=2))]:
   cache=self.root/name
   with self.assertRaises(ValueError):stream.index_source(self.source,cache,**kwargs)
   self.assertFalse((cache/'index.json').exists())
 def test_changed_component_rejected_even_if_external_index_is_resealed(self):
  self.write(self.fixture());actual=stream.index_source(self.source,self.cache);actual.clear_cache()
  index=json.loads((self.cache/'index.json').read_text('utf-8'));d=index['children']['complete'];path=self.cache/d['file']
  path.write_text('false','utf-8');d['sha256']=stream.sha(path);d['bytes']=5
  (self.cache/'index.json').write_text(json.dumps(index),'utf-8')
  with self.assertRaisesRegex(ValueError,'component changed'):actual['complete']
 def test_changed_original_source_rejected(self):
  self.write(self.fixture());actual=stream.index_source(self.source,self.cache)
  self.source.write_text('{}','utf-8')
  with self.assertRaisesRegex(ValueError,'original source changed'):actual.verify_source()
 def test_invalid_numbers_root_trailing_and_nesting_are_rejected(self):
  values=['{"a":1e400}','{"a":1e101}','[]','{"a":1} []','{"a":['*130+'0'+']}'*130]
  for i,raw in enumerate(values):
   self.source.write_text(raw,'utf-8');cache=self.root/str(i)
   with self.assertRaises(ValueError):stream.index_source(self.source,cache)
   self.assertFalse((cache/'index.json').exists())
 def test_encoded_tokens_bounded_even_across_read_chunks_and_escape_boundary(self):
  reader=stream.CheckedReader(io.BytesIO(b'{"text":"a\\\"b\\\\c"}'),100)
  while reader.read(1):pass
  self.assertFalse(reader.inside)
  reader=stream.CheckedReader(io.BytesIO(b'{"x":"'+b'a'*(1024**2+1)+b'"}'),2*1024**2)
  with self.assertRaisesRegex(ValueError,'string token'):
   while reader.read(65536):pass
  reader=stream.CheckedReader(io.BytesIO(b'{"x":'+b'1'*4097+b'}'),10000)
  with self.assertRaisesRegex(ValueError,'numeric token'):
   while reader.read(100):pass

 def test_literal_dots_item_collision_empty_and_unicode_map_keys_are_preserved(self):
  expected=self.fixture()
  expected['extra.with.dot']={'name':'地盤'}
  expected['frames.item']={'collision':True}
  expected['']={'':{'inside':3}}
  expected['geostatic']['diagnostic.with.dot']={'value':2}
  expected['geostatic']['']={'':4}
  self.write(expected);actual=stream.index_source(self.source,self.cache)
  self.assertEqual(materialize(actual),expected)
  actual.verify_source()

if __name__=='__main__':unittest.main()
