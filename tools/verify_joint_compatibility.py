"""Check unchanged contact specimen physics after adding rigid-mode reports."""
import argparse
import hashlib
import json
from pathlib import Path


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(baseline,current,output):
    root=Path(__file__).resolve().parents[1]
    accepted_path=baseline/'final-validation.json';accepted=json.loads(accepted_path.read_text())
    assert accepted['passed'];records=[];files={str(accepted_path):sha(accepted_path),str(Path(__file__).resolve()):sha(Path(__file__))}
    for name in ('mc','mcc'):
        ap=baseline/name/'stage-results.json';bp=current/name/'stage-results.json'
        assert sha(ap)==accepted['sha256'][str(ap.relative_to(root))]
        a=json.loads(ap.read_text());b=json.loads(bp.read_text());assert a['complete'] and b['complete']
        assert a['case']==b['case'] and a['ground']==b['ground']
        for field in ('last_committed','reference_displacement_m','integration_points','joint_integration_points','excavation_audits'):
            assert a[field]==b[field],(name,field)
        frames=[]
        for f in b['frames']:
            assert f['joint_support']['unrestrained_modes']==0
            frame=f.copy()
            if 'joint_support' not in a['frames'][len(frames)]:frame.pop('joint_support')
            frames.append(frame)
        assert a['frames']==frames,(name,'frame physics changed')
        records.append({'material':name,'exact_committed_state':True,'exact_frames_except_added_support_report':True,
                        'baseline_solver_sha256':a['provenance']['solver_binary_sha256'],'current_solver_sha256':b['provenance']['solver_binary_sha256']})
        files.update({str(p):sha(p) for p in (ap,bp)})
    report={'schema':'soil-fem-joint-compatibility/1','passed':True,'comparisons':records,'sha256':files}
    with output.open('x',encoding='utf-8') as stream:json.dump(report,stream,indent=2)
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('baseline',type=Path);p.add_argument('current',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();verify(a.baseline.resolve(),a.current.resolve(),a.output.resolve())
