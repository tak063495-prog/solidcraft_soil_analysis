"""Regression: delegated audits must not hold a second materialized GP result."""
from pathlib import Path
import codecs,json,sys,tempfile,unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
import verify_stage_package as package

class ReadOnce:
 def __init__(self,path):self.path=path;self.full_reads=0
 def open(self,*a,**kw):return self.path.open(*a,**kw)
 def stat(self):return self.path.stat()
 def read_text(self,*a,**kw):
  self.full_reads+=1
  if self.full_reads>1:raise AssertionError('same large GP result materialized twice across audit dispatch')
  return self.path.read_text(*a,**kw)

class ProjectedAuditDispatch(unittest.TestCase):
 def test_delegated_mixed_result_materialized_only_once_even_schema_last_or_BOM(self):
  schemas=['soil-fem-mixed-mcc-stage-result/1','soil-fem-mixed-mcc-stage-result/2','soil-fem-projected-mcc-stage-result/1','soil-fem-projected-mcc-stage-result/2']
  for schema in schemas:
   for bom in [False,True]:
    with self.subTest(schema=schema,bom=bom),tempfile.TemporaryDirectory() as folder:
     folder=Path(folder);mesh=folder/'mesh.json';ground=folder/'ground.json';result=folder/'result.json'
     mesh.write_text('{}');ground.write_text('{}')
     expected={'payload':[{'stress':[i,2.,3.,0.,0.,0.],'pc':120.} for i in range(1500)],'schema':schema}
     result.write_bytes((codecs.BOM_UTF8 if bom else b'')+json.dumps(expected).encode('utf-8'))
     counted=ReadOnce(result)
     def delegated(a,b,c,d,solver=None):
      self.assertIs(c,counted)
      value=json.loads(c.read_text(encoding='utf-8-sig'));self.assertEqual(value,expected)
      return {'checked_payload_entries':len(value['payload'])}
     with patch('verify_projected_stage_package.verify',side_effect=delegated):
      self.assertEqual(package.verify(mesh,ground,counted,folder),{'checked_payload_entries':1500})
     self.assertEqual(counted.full_reads,1)

if __name__=='__main__':unittest.main()
