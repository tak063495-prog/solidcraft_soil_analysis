import json,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
import export_solidcraft_view as exporter

class StreamedExportSourceTests(unittest.TestCase):
 def setUp(self):
  self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)
  self.source=self.root/'source.json'
 def fixture(self):
  return dict(schema='soil-fem-projected-mcc-stage-result/2',complete=True,failure_reason=None,
      frames=[{'name':'initial','state':{'material_states':[{'void_ratio':.8}]}}],
      integration_points=[{'element_index':0}],geostatic=None,last_committed={'marker':1})
 def test_large_input_uses_chunked_source_and_preserves_final_state(self):
  expected=self.fixture();self.source.write_text(json.dumps(expected),'utf-8')
  with patch.object(exporter,'MAX_SOURCE_RESULT',16):run=exporter.read_result(self.source,cache=self.root/'cache')
  self.assertNotIsInstance(run,dict)
  self.assertEqual(run['frames'][0],expected['frames'][0]);self.assertEqual(run['last_committed'],expected['last_committed'])
  run.verify_source()
 def test_large_input_requires_explicit_cache_and_projected_family(self):
  expected=self.fixture();self.source.write_text(json.dumps(expected),'utf-8')
  with patch.object(exporter,'MAX_SOURCE_RESULT',16),self.assertRaises(ValueError):exporter.read_result(self.source)
  expected['schema']='soil-fem-stage-result/13';self.source.write_text(json.dumps(expected),'utf-8')
  with patch.object(exporter,'MAX_SOURCE_RESULT',16),self.assertRaisesRegex(ValueError,'projected MCC'):
   exporter.read_result(self.source,cache=self.root/'cache')
  self.assertFalse((self.root/'cache').exists())
 def test_existing_small_source_keeps_original_dictionary_reader(self):
  expected=self.fixture();self.source.write_text(json.dumps(expected),'utf-8');run=exporter.read_result(self.source)
  self.assertIsInstance(run,dict);self.assertEqual(run,expected)
 def test_blocked_scene_clean_keeps_arrays_strings_null_and_missing_semantics(self):
  import numpy as np
  from streamed_source_v2 import index_source
  expected={'frame':{'values':list(range(13)),'name':'xx'},'unavailable':None}
  self.source.write_text(json.dumps(expected),'utf-8')
  run=index_source(self.source,self.root/'cache',block_items=2)
  self.assertEqual(exporter.clean(run),expected)
  self.assertEqual(exporter.clean({'x':np.asarray([1.,np.nan,np.inf]),'note':'MCC'}),{'x':[1.,None,None],'note':'MCC'})
 def test_streamed_encoded_admission_does_not_raise_component_or_display_caps(self):
  import streamed_source_v2 as blocked
  import streamed_stage_source as legacy
  expected=self.fixture();self.source.write_text(json.dumps(expected),'utf-8')
  with patch.object(exporter,'MAX_SOURCE_RESULT',16),patch.object(blocked,'MAX_STREAM_SOURCE',20),self.assertRaisesRegex(ValueError,'source limit'):
   exporter.read_result(self.source,cache=self.root/'cache')
  self.assertFalse((self.root/'cache').exists())
  self.assertEqual(exporter.MAX_FILE,256*1024**2);self.assertEqual(exporter.MAX_TOTAL,512*1024**2)
  self.assertEqual(legacy.MAX_COMPONENT,256*1024**2);self.assertEqual(legacy.MAX_PAYLOAD,512*1024**2)
  self.assertEqual(blocked.MAX_STREAM_SOURCE,32*1024**3)

if __name__=='__main__':unittest.main()
