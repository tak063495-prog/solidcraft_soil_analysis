"""CLI provenance and file-preservation contracts for material-state diagnostics."""
from pathlib import Path
import hashlib,json,math,os,subprocess,sys,tempfile,unittest
ROOT=Path(__file__).resolve().parents[1]

class Cli(unittest.TestCase):
    def inputs(self,folder):
        material=dict(critical_state_m=1.4,**{'lambda':.12},kappa=.025,poisson=.2)
        q=1.4*math.sqrt(100*300)
        state=dict(effective_stress_compression_positive_kpa=[100+q/3,100+q/3,100-2*q/3,0.,0.,0.],preconsolidation_pressure_kpa=400.,void_ratio=.9)
        mesh=dict(nodes=[[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],[0.,0.,-1.]],tetrahedra=[[0,1,2,3],[0,2,1,4]],material_ids=[1,1],region_ids=[1,1],regions=[dict(id=1,body='TestBody')])
        result=dict(complete=True,integration_points=[dict(element_index=0,position_m=[.25,.25,.25],weight_m3=2.),dict(element_index=1,position_m=[.25,.25,-.25],weight_m3=3.)],frames=[dict(name='material_inventory',state=dict(active_elements=[True,False],material_states=[state,state]))])
        values=dict(mesh=mesh,ground=dict(materials=[dict(id=1,mcc=material)]),result=result)
        paths={}
        for name,value in values.items():
            paths[name]=folder/(name+'.json');paths[name].write_text(json.dumps(value),encoding='utf-8')
        return paths

    def call(self,name,*args):
        env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONUTF8='1')
        return subprocess.run([sys.executable,str(ROOT/'tools'/name),*map(str,args)],capture_output=True,text=True,env=env)

    def test_bound_cli_reports_and_changed_input_rejection(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder=Path(temporary);paths=self.inputs(folder);before={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths.values()}
            sampled=folder/'sampled.json';exact=folder/'exact.json'
            first=self.call('mcc_acoustic_diagnostic.py',paths['mesh'],paths['ground'],paths['result'],sampled)
            self.assertEqual(first.returncode,0,first.stderr)
            self.assertEqual(json.loads(sampled.read_text())['inputs_sha256'],before)
            second=self.call('mcc_acoustic_exact_report.py',sampled,exact)
            self.assertEqual(second.returncode,0,second.stderr)
            report=json.loads(exact.read_text());self.assertEqual(report['all_direction_negative_candidates'],1)
            self.assertEqual(report['negative_volume_m3'],2.);self.assertEqual(report['cad_body_groups']['TestBody']['negative_candidates'],1)
            self.assertEqual({str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths.values()},before)
            paths['result'].write_text('{}',encoding='utf-8')
            rejected=self.call('mcc_acoustic_exact_report.py',sampled,folder/'rejected.json')
            self.assertNotEqual(rejected.returncode,0);self.assertFalse((folder/'rejected.json').exists())

    def test_existing_outputs_are_preserved_by_both_clis(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder=Path(temporary);paths=self.inputs(folder);output=folder/'existing.json';output.write_bytes(b'keep-this-file')
            result=self.call('mcc_acoustic_diagnostic.py',paths['mesh'],paths['ground'],paths['result'],output)
            self.assertNotEqual(result.returncode,0);self.assertEqual(output.read_bytes(),b'keep-this-file')
            result=self.call('mcc_acoustic_exact_report.py',folder/'not-read.json',output)
            self.assertNotEqual(result.returncode,0);self.assertEqual(output.read_bytes(),b'keep-this-file')

if __name__=='__main__':unittest.main()
