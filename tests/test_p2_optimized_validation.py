"""Real subprocess exports/unpack preserve external-data checks under -O."""
import hashlib,json,os,subprocess,sys,tempfile,unittest,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/'interim/2026-10-27-native-mcc-selector/unpack_inputs.py'
CODE=r'''
import json,sys
from pathlib import Path
root=Path(sys.argv[1]);sys.path[:0]=[str(root/'tests'),str(root/'tools')]
from test_solidcraft_structure_fields import scene_fixture
from test_solidcraft_view import job_fixture,write
from verify_acceptance_packet import sha
from export_solidcraft_view import export
out=Path(sys.argv[2]);mutation=sys.argv[3];mesh,run=scene_fixture();frame=run['frames'][0]
if mutation=='moment':frame['beam_results']['plates'][0]['samples'][0]['bending_moment_kn'][0]+=.1
elif mutation=='axes':frame['beam_results']['members'][0]['local_axes'][0][0]+=1.
elif mutation=='force':frame['beam_results']['external_force_kn'][0]+=100.
elif mutation=='history':
 frame['state']['beams']['node_displacement_m'][3][0]+=.1
 frame['beam_results']['node_displacement_m'][3][0]+=.1
job=out.parent/'job';job_fixture(job,(mesh,run));export(job,out)
'''
class OptimizedChecks(unittest.TestCase):
    def call(self,mode,args):
        env=dict(os.environ);env.pop('PYTHONOPTIMIZE',None)
        if mode=='env':env['PYTHONOPTIMIZE']='1'
        return subprocess.run([sys.executable,*(['-O'] if mode=='-O' else []),*map(str,args)],capture_output=True,text=True,env=env)
    def test_structural_normal_and_mismatches_never_publish_manifest(self):
        for mode in ('normal','-O','env'):
            for mutation in ('normal','moment','axes','force','history'):
                with self.subTest(mode=mode,mutation=mutation),tempfile.TemporaryDirectory() as d:
                    out=Path(d)/'view';p=self.call(mode,['-c',CODE,ROOT,out,mutation])
                    if mutation=='normal':self.assertEqual(p.returncode,0,p.stderr);self.assertTrue((out/'view-manifest.json').exists())
                    else:self.assertNotEqual(p.returncode,0,p.stdout);self.assertFalse((out/'view-manifest.json').exists())
    def test_unpack_inventory_bytes_and_hash_in_all_modes(self):
        self.assertTrue(SCRIPT.is_file())
        expected=['mesh.json','ground.json','mesh-identity.json','stages-coarse.json','stages-fine.json']
        for mode in ('normal','-O','env'):
            for mutation in ('normal','hash','size','extra','duplicate'):
                with self.subTest(mode=mode,mutation=mutation),tempfile.TemporaryDirectory() as d:
                    folder=Path(d);script=folder/'unpack_inputs.py';script.write_bytes(SCRIPT.read_bytes());out=folder/'out'
                    manifest={n:dict(bytes=2,sha256=hashlib.sha256(b'{}').hexdigest()) for n in expected}
                    (folder/'input-manifest.json').write_text(json.dumps(manifest))
                    with zipfile.ZipFile(folder/'inputs.zip','w') as z:
                        for n in expected:z.writestr(n,b'[]' if n=='mesh.json' and mutation=='hash' else b'{} ' if n=='mesh.json' and mutation=='size' else b'{}')
                        if mutation=='extra':z.writestr('extra.json',b'{}')
                        if mutation=='duplicate':z.writestr('mesh.json',b'{}')
                    p=self.call(mode,[script,out])
                    if mutation=='normal':self.assertEqual(p.returncode,0,p.stderr);self.assertTrue((out/'fine/ground/mesh.json').is_file())
                    else:self.assertNotEqual(p.returncode,0,p.stdout);self.assertFalse(out.exists());self.assertNotIn('Exact verification inputs created',p.stdout)

if __name__=='__main__':unittest.main()
