"""Committed material inventory and optimized exact-report binding regressions."""
import copy,json,os,subprocess,sys,tempfile,unittest
from pathlib import Path
import test_mcc_acoustic_cli as acoustic_cli
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from mcc_acoustic_diagnostic import diagnose

class MaterialInventory(unittest.TestCase):
    def values(self,folder):
        paths=acoustic_cli.Cli().inputs(folder);values={k:json.loads(p.read_text()) for k,p in paths.items()}
        old=values['ground']['materials'][0];new=copy.deepcopy(old);new['id']=2;old['mcc']['critical_state_m']=1.6
        values['ground']['materials'].append(new)
        before=values['result']['frames'][0];before['state']['element_material_ids']=[1,1]
        after=copy.deepcopy(before);after['name']='after_matset';after['state']['element_material_ids']=[2,1];values['result']['frames'].append(after)
        return paths,values
    def call(self,script,args,mode='normal'):
        env=dict(os.environ);env.pop('PYTHONOPTIMIZE',None)
        if mode=='env':env['PYTHONOPTIMIZE']='1'
        return subprocess.run([sys.executable,*(['-O'] if mode=='-O' else []),str(ROOT/'tools'/script),*map(str,args)],capture_output=True,text=True,env=env)
    def test_matset_before_after_and_exact_cli(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);paths,v=self.values(folder)
            report=diagnose(v['mesh'],v['ground'],v['result'])
            self.assertEqual([r['near_yield_candidates'] for r in report['frames']],[0,1])
            self.assertEqual([r['sampled_negative_acoustic_candidates'] for r in report['frames']],[0,1])
            for name,p in paths.items():p.write_text(json.dumps(v[name]))
            sampled=folder/'sampled.json';exact=folder/'exact.json'
            p=self.call('mcc_acoustic_diagnostic.py',[*paths.values(),sampled]);self.assertEqual(p.returncode,0,p.stderr)
            p=self.call('mcc_acoustic_exact_report.py',[sampled,exact]);self.assertEqual(p.returncode,0,p.stderr)
            self.assertEqual(json.loads(exact.read_text())['all_direction_negative_candidates'],1)
    def test_material_constants_same_id_and_legacy_fallback(self):
        with tempfile.TemporaryDirectory() as d:
            _,v=self.values(Path(d));last=v['result']['frames'][-1]
            for key,new in (('lambda',.15),('kappa',.03),('poisson',.25)):
                changed=copy.deepcopy(v);changed['ground']['materials'][1]['mcc'][key]=new
                expected=copy.deepcopy(changed);expected['mesh']['material_ids']=[2,1]
                for f in expected['result']['frames']:f['state'].pop('element_material_ids')
                self.assertEqual(diagnose(changed['mesh'],changed['ground'],dict(changed['result'],frames=[changed['result']['frames'][-1]]))['frames'][0],diagnose(expected['mesh'],expected['ground'],dict(expected['result'],frames=[expected['result']['frames'][-1]]))['frames'][0])
            for ids in (None,[],[1,1]):
                data=copy.deepcopy(v);data['result']['frames']=[data['result']['frames'][0]]
                if ids is None:data['result']['frames'][0]['state'].pop('element_material_ids')
                else:data['result']['frames'][0]['state']['element_material_ids']=ids
                self.assertEqual(diagnose(data['mesh'],data['ground'],data['result'])['frames'][0]['near_yield_candidates'],0)
    def test_bad_inventory_rejected_including_inactive_elements(self):
        with tempfile.TemporaryDirectory() as d:
            _,v=self.values(Path(d))
            for ids in ([2],[2,1,1],[2,999],[2,True],[2,1.0],['2',1],None,{},''):
                data=copy.deepcopy(v);data['result']['frames'][-1]['state']['element_material_ids']=ids
                with self.subTest(ids=ids),self.assertRaises(ValueError):diagnose(data['mesh'],data['ground'],data['result'])
    def test_exact_matset_legacy_and_invalid_inventory_in_all_modes(self):
        import hashlib
        for mode in ('normal','-O','env'):
            for ids,expected in (([2,1],1),([1,1],0),([],0),('missing',0),([2],None),([2,999],None),([True,1],None)):
                with self.subTest(mode=mode,ids=ids),tempfile.TemporaryDirectory() as d:
                    folder=Path(d);paths,v=self.values(folder)
                    if expected is not None:
                        if ids=='missing':v['result']['frames'][-1]['state'].pop('element_material_ids')
                        else:v['result']['frames'][-1]['state']['element_material_ids']=ids
                    for name,p in paths.items():p.write_text(json.dumps(v[name]))
                    sampled=folder/'sampled.json';out=folder/'exact.json'
                    p=self.call('mcc_acoustic_diagnostic.py',[*paths.values(),sampled]);self.assertEqual(p.returncode,0,p.stderr)
                    if expected is None:
                        v['result']['frames'][-1]['state']['element_material_ids']=ids;paths['result'].write_text(json.dumps(v['result']))
                        data=json.loads(sampled.read_text());data['inputs_sha256'][str(paths['result'].resolve())]=hashlib.sha256(paths['result'].read_bytes()).hexdigest();sampled.write_text(json.dumps(data))
                    p=self.call('mcc_acoustic_exact_report.py',[sampled,out],mode)
                    if expected is None:self.assertNotEqual(p.returncode,0);self.assertFalse(out.exists())
                    else:self.assertEqual(p.returncode,0,p.stderr);self.assertEqual(json.loads(out.read_text())['all_direction_negative_candidates'],expected)

    def test_exact_changed_input_and_implementation_rejected_in_all_modes(self):
        for mode in ('normal','-O','env'):
            for mutation in ('ground','implementation','shape','weight','candidate_count'):
                with self.subTest(mode=mode,mutation=mutation),tempfile.TemporaryDirectory() as d:
                    folder=Path(d);paths=acoustic_cli.Cli().inputs(folder);sampled=folder/'sampled.json';output=folder/'exact.json'
                    p=self.call('mcc_acoustic_diagnostic.py',[*paths.values(),sampled]);self.assertEqual(p.returncode,0,p.stderr)
                    if mutation=='ground':
                        v=json.loads(paths['ground'].read_text());v['materials'][0]['mcc']['lambda']=.13;paths['ground'].write_text(json.dumps(v))
                    else:
                        data=json.loads(sampled.read_text())
                        if mutation=='implementation':data['implementation_sha256']='0'*64
                        elif mutation=='candidate_count':data['frames'][-1]['near_yield_candidates']+=1
                        else:
                            import hashlib
                            v=json.loads(paths['result'].read_text())
                            if mutation=='shape':v['frames'][-1]['state']['active_elements']=[[True,False]]
                            else:v['integration_points'][0]['weight_m3']=-1
                            paths['result'].write_text(json.dumps(v));data['inputs_sha256'][str(paths['result'].resolve())]=hashlib.sha256(paths['result'].read_bytes()).hexdigest()
                        sampled.write_text(json.dumps(data))
                    p=self.call('mcc_acoustic_exact_report.py',[sampled,output],mode);self.assertNotEqual(p.returncode,0,p.stdout);self.assertFalse(output.exists())

if __name__=='__main__':unittest.main()
