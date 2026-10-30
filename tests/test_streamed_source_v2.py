from pathlib import Path
import json,tempfile,unittest
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
import numpy as np
import streamed_source_v2 as mod

class Blocks(unittest.TestCase):
 def source(self,folder,value):
  path=Path(folder)/'source.json';path.write_text(json.dumps(value),'utf-8');return path
 def test_every_value_and_type_and_previous_history_preserved(self):
  with tempfile.TemporaryDirectory() as folder:
   state={'material_states':[{'stress':[i,2,3,0,0,0],'e':1.3} for i in range(31)],'mixed_mcc':{'theta':list(range(31))}}
   original={'schema':'test','integration_points':[{'owner':i//3} for i in range(31)],'case':{'short':[1,2]},'frames':[{'state':state},{'state':state}], 'last_committed':state,'geostatic':{'seed':state},'transition_audits':[{'prepared_state':state}]}
   source=self.source(folder,original);run=mod.index_source(source,Path(folder)/'cache',block_items=4,cache_payload=512)
   self.assertEqual(dict(run),original);self.assertIsInstance(run['frames'][0],dict)
   self.assertIsInstance(run['frames'][0]['state']['mixed_mcc'],dict)
   self.assertEqual(run['last_committed'],run['frames'][-1]['state'])
   self.assertIsInstance(run['case']['short'],list)
   self.assertTrue(np.array_equal(np.asarray(run['frames'][0]['state']['mixed_mcc']['theta']),np.arange(31)))
   self.assertLessEqual(run._store.retained,512);run.verify_source()
 def test_array_negative_slice_bounds_and_cache_eviction(self):
  with tempfile.TemporaryDirectory() as folder:
   source=self.source(folder,{'v':list(range(25))});run=mod.index_source(source,Path(folder)/'cache',block_items=4,cache_payload=512)
   values=run['v'];self.assertEqual(values[-1],24);self.assertEqual(values[2:19:3],list(range(25))[2:19:3])
   with self.assertRaises(IndexError):_=values[25]
   with self.assertRaises(IndexError):_=values[-26]
   self.assertEqual(list(values),list(range(25)));self.assertLessEqual(run._store.retained,512)
 def test_duplicate_nested_and_literal_dot_keys(self):
  for raw in ['{"a":{"v":1,"v":2}}','{"v":[{"stress":1,"stress":2}]}']:
   with self.subTest(raw=raw),tempfile.TemporaryDirectory() as folder:
    source=Path(folder)/'source.json';source.write_text(raw)
    with self.assertRaisesRegex(ValueError,'duplicate'):mod.index_source(source,Path(folder)/'cache')
    self.assertFalse((Path(folder)/'cache/index.json').exists())
  with tempfile.TemporaryDirectory() as folder:
   original={'frames.item':{'':3,'item':[1,2]},'frames':[{'state':{'a.b':[2,4,6]}}]}
   source=self.source(folder,original);self.assertEqual(dict(mod.index_source(source,Path(folder)/'cache',block_items=2)),original)
 def test_one_entry_over_payload_fails_without_complete_index(self):
  with tempfile.TemporaryDirectory() as folder:
   source=self.source(folder,{'v':[[1]*40]})
   with self.assertRaisesRegex(ValueError,'entry exceeds'):mod.index_source(source,Path(folder)/'cache',payload_limit=512)
   self.assertFalse((Path(folder)/'cache/index.json').exists())
 def test_large_array_is_split_without_raising_payload_cap(self):
  with tempfile.TemporaryDirectory() as folder:
   source=self.source(folder,{'v':list(range(100))});run=mod.index_source(source,Path(folder)/'cache',block_items=4,payload_limit=512,cache_payload=512)
   self.assertEqual(run['v'],list(range(100)))
   self.assertTrue(all(v['payload']<=512 for v in run._store.index['components']))
 def test_component_and_source_tamper_and_oversize(self):
  for target in ['component','source']:
   with self.subTest(target=target),tempfile.TemporaryDirectory() as folder:
    source=self.source(folder,{'v':list(range(20))});run=mod.index_source(source,Path(folder)/'cache',block_items=4)
    path=source if target=='source' else Path(folder)/'cache'/run._store.index['components'][0]['file']
    raw=path.read_bytes();path.write_bytes(raw.replace(b'0',b'9',1));run.clear_cache()
    with self.assertRaisesRegex(ValueError,'changed'):run.verify_source()
  with tempfile.TemporaryDirectory() as folder:
   source=self.source(folder,{'v':list(range(20))})
   with self.assertRaisesRegex(ValueError,'encoded byte'):mod.index_source(source,Path(folder)/'cache',source_limit=10)
 def test_nonfinite_trailing_and_incomplete(self):
  for raw in ['{"v":1e309}','{"v":1}{}','{"v":[1,2']:
   with self.subTest(raw=raw),tempfile.TemporaryDirectory() as folder:
    source=Path(folder)/'source.json';source.write_text(raw)
    with self.assertRaises(ValueError):mod.index_source(source,Path(folder)/'cache')
    self.assertFalse((Path(folder)/'cache/index.json').exists())
 def test_frame_and_transition_caps_remain_distinct(self):
  for name,count,allowed in [('frames',258,False),('transition_audits',258,True),('transition_audits',513,False)]:
   with self.subTest(name=name,count=count),tempfile.TemporaryDirectory() as folder:
    source=self.source(folder,{name:[{'id':i} for i in range(count)]})
    if allowed:self.assertEqual(len(mod.index_source(source,Path(folder)/'cache')[name]),count)
    else:
     with self.assertRaisesRegex(ValueError,'too many'):mod.index_source(source,Path(folder)/'cache')

if __name__=='__main__':unittest.main()
