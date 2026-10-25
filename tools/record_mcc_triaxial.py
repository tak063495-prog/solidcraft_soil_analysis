"""Seal raw coarse/fine reference, native paths, model source and executable."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile
import numpy as np


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))


def record(run):
    root=Path(__file__).resolve().parents[1];paths=run/'paths';report=read(paths/'verification.json')
    assert report['passed'] and report['solver_sha256']==sha(run/'solver/soil-fem-3d.exe')
    assert (run/'verification.log').read_text(encoding='utf-8').rstrip().endswith('Exit code: 0')
    for name,digest in report['sha256'].items():assert sha(paths/name)==digest,name
    assert sha(paths/'reference-code.py')==report['reference_code_sha256']
    assert {r['case'] for r in report['cases']}=={'nc-drained','oc-drained','nc-undrained','oc-undrained'}
    for row in report['cases']:
        name=row['case'];coarse=read(paths/(name+'-coarse-reference.json'));fine=read(paths/(name+'-fine-reference.json'))
        expected=np.asarray(fine['states']);refinement=np.linalg.norm(np.asarray(coarse['states'])[:,:3]-expected[:,:3])/np.linalg.norm(expected[:,:3])
        assert refinement<1e-6 and np.isclose(refinement,row['reference_refinement_error'],rtol=1e-12,atol=1e-15)
        errors=[]
        for suffix in ('coarse','fine'):
            data=read(paths/(name+'-'+suffix+'-result.json'));actual=[]
            for s in data['states'][1:]:
                v=np.asarray(s['effective_stress_compression_positive_kpa']);p=v[:3].mean();dev=v.copy();dev[:3]-=p
                q=np.sqrt(1.5*(dev[:3]@dev[:3]+2*(dev[3:]@dev[3:])))
                actual.append([p,q,s['preconsolidation_pressure_kpa'],s['void_ratio']])
            actual=np.asarray(actual);errors.append(np.linalg.norm(actual[:,:3]-expected[:,:3])/np.linalg.norm(expected[:,:3]))
            assert (paths/(name+'-'+suffix+'.log')).read_text(encoding='utf-8').rstrip().endswith('Exit code: 0')
        assert errors[1]<errors[0] and errors[1]<2e-5
        assert np.allclose(errors,row['native_relative_p_q_pc_error_coarse_fine'],rtol=1e-12,atol=1e-15)
        assert np.allclose(actual[-1],row['final_p_q_pc_e'],rtol=1e-12,atol=1e-12)
    sources=[root/p for p in ('Cargo.lock','Cargo.toml','rust-toolchain.toml','README.md')]
    sources += [p for name in ('src','tests','examples','scripts','docs') for p in (root/name).rglob('*') if p.is_file() and p.suffix in ('.rs','.py','.json','.ps1','.md')]
    sources += list((root/'tools').glob('*.py'));manifest={str(p.relative_to(root)):sha(p) for p in sorted(set(sources))}
    with zipfile.ZipFile(run/'source-snapshot.zip','x',compression=zipfile.ZIP_DEFLATED) as z:
        for name,digest in manifest.items():
            payload=(root/name).read_bytes();assert hashlib.sha256(payload).hexdigest()==digest;z.writestr(name,payload)
    (run/'source-snapshot.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    final={'schema':'soil-fem-mcc-triaxial-final-validation/1','passed':True,'solver_sha256':report['solver_sha256'],
        'cases':report['cases'],'limitations':report['limitations'],
        'source_snapshot_sha256':sha(run/'source-snapshot.zip'),
        'sha256':{str(p.relative_to(root)):sha(p) for p in sorted(run.rglob('*')) if p.is_file()}}
    with (run/'final-validation.json').open('x',encoding='utf-8') as stream:json.dump(final,stream,indent=2)
    print('MCC triaxial source and evidence sealed:',run)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args();record(a.run.resolve())
