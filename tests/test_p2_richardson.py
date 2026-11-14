"""Exercise the actual plane-strain judge CLI with signed and guarded series."""
import json,os,subprocess,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BINARY=ROOT/'target/release'/('soil-fem-3d.exe' if os.name=='nt' else 'soil-fem-3d')
class SignedRichardson(unittest.TestCase):
    def judge(self,values,reference):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);source=folder/'series.json';out=folder/'report.json'
            source.write_text(json.dumps(dict(schema='soil-fem-plane-strain-2d-series/1',quantity='signed displacement',reference=reference,values=values,nested_h_refinement=True)))
            p=subprocess.run([str(BINARY),'plane-strain-2d','judge',str(source),str(out)],capture_output=True,text=True)
            if any(not __import__('math').isfinite(x) for x in values):self.assertNotEqual(p.returncode,0);self.assertFalse(out.exists());return None
            self.assertEqual(p.returncode,0,p.stderr);return json.loads(out.read_text())['judgment']
    def test_oscillating_positive_and_negative_are_indeterminate(self):
        for sign in (-1,1):
            report=self.judge([sign*x for x in (.001,.0013,.00115)],sign*.001)
            self.assertEqual(report['continuum_verdict'],'indeterminate');self.assertIsNone(report['reference_free_estimate'])
    def test_monotone_known_orders_and_existing_guards(self):
        for sign in (-1,1):
            for order in (1,2,3):
                values=[sign*(.001+.0002/2**(order*n)) for n in range(3)]
                report=self.judge(values,sign*.001);self.assertAlmostEqual(report['reference_free_estimate'],sign*.001);self.assertAlmostEqual(report['reference_free_order'],order);self.assertEqual(report['continuum_verdict'],'pass')
        for values in ([.002,.002,.002],[.002,.003,.005],[.001,.0012],[.001,.0012,.0012]):
            self.assertEqual(self.judge(values,.001)['continuum_verdict'],'indeterminate')
        self.assertEqual(self.judge([.001,.001,.001],.001)['continuum_verdict'],'pass')
        for value in (float('nan'),float('inf'),-float('inf')):self.judge([.001,.0011,value],.001)
if __name__=='__main__':unittest.main()
