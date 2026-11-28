"""The view budgets cover exactly the emitted geometry/JSON/manifest bytes."""
import json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
from test_mc_epoch_view import setup
import export_mc_epoch_view as module

class Budgets(unittest.TestCase):
    def test_actual_file_and_total_boundaries_and_no_success_on_overflow(self):
        with tempfile.TemporaryDirectory() as d:
            root,job,native,accepted=setup(d)
            with patch.object(module,'accepted',return_value=accepted):module.export(job,native,root/'control')
            actual={p.name:p.read_bytes() for p in (root/'control/view').iterdir()}
            total=sum(map(len,actual.values()));largest=max(map(len,actual.values()))
            for kind,limit,ok in (('file',largest,True),('file',largest-1,False),('total',total,True),('total',total-1,False)):
                out=root/(kind+str(limit))
                with self.subTest(kind=kind,limit=limit),patch.object(module,'accepted',return_value=accepted),patch.object(module,'MAX_FILE',limit if kind=='file' else 1024**2),patch.object(module,'MAX_TOTAL',limit if kind=='total' else 1024**2):
                    if ok:
                        module.export(job,native,out)
                        self.assertEqual({p.name:p.read_bytes() for p in (out/'view').iterdir()},actual)
                    else:
                        with self.assertRaises(ValueError):module.export(job,native,out)
                        self.assertFalse((out/'view/view-manifest.json').exists());self.assertFalse((out/'view-verification.json').exists())
    def test_every_payload_uses_same_encoded_bytes_for_size_and_write(self):
        budget_type=getattr(module,'ViewBudget')
        for name,value in (('frame',dict(a=[[0.]*6]*100)),('points',dict(a=[[1.,2.,3.]]*100)),('manifest',dict(complete=True,files=['a','b'])),('geometry',dict(nodes=[[0.,0.,0.]]*100))):
            data=(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n').encode('utf-8')
            with self.subTest(name=name),patch.object(module,'MAX_FILE',len(data)),patch.object(module,'MAX_TOTAL',len(data)):
                budget=budget_type();self.assertEqual(budget.add(data),data)
                with self.assertRaises(ValueError):budget.add(b'x')
            with patch.object(module,'MAX_FILE',len(data)-1):
                with self.assertRaises(ValueError):budget_type().add(data)

if __name__=='__main__':unittest.main()
