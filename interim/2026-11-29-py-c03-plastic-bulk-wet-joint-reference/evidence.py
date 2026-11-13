"""Read original records, select evidence with exact JSON pointers and hashes."""
import argparse
import hashlib
import json
from pathlib import Path

BASE='3b1133d3be9905dc24f62c27b5b1eb456a83f2e2'
OLD='interim/2026-11-20-wet-joint-plastic-bulk-integration/'
EXPECTED_HASHES=[
    '25ab2cb08799ad3e5a77d1286bab9083669acc49b7c0cefdfc1e1a55316ce685',
    'be75f689746c82ff1aa548f6bfe3534f724de47055917abf5cee104ae3f8686e',
    'ae886c201114b66824be2297a836a338f86448750a2d91d9981c6cd8cf5a6961',
    'b0027189fa5e66dd494e50ceae357b27e003c6ba5331c14cfd43a3ce19c8a9ab',
    '77c93e89ac56576075f7252b3e0bc9386da922cbbbe1a0e20ebbe84df8a05708']


def extract(root,out):
    paths=[OLD+'README.ja.md',OLD+'measurements.json',OLD+'conclusion.json',
           'interim/2026-11-16-issue-triage/todo.ja.md',
           'interim/2026-11-15-wet-joint-convergence/acceptance.json']
    sources={p:dict(sha256=hashlib.sha256((root/p).read_bytes()).hexdigest()) for p in paths}
    for p,expected in zip(paths,EXPECTED_HASHES):
        if sources[p]['sha256']!=expected:
            raise ValueError('original specified-commit record changed: '+p)
    m=json.loads((root/paths[1]).read_text(encoding='utf-8'))
    c=json.loads((root/paths[2]).read_text(encoding='utf-8'))
    a=json.loads((root/paths[4]).read_text(encoding='utf-8'))
    selections={}
    for name in sorted(m['native_vs_independent']['baseline']['all_fields']):
        pointer='/native_vs_independent/baseline/all_fields/'+name
        selections[pointer]=m['native_vs_independent']['baseline']['all_fields'][name]
    for key,value in m['sensitivity'].items():
        selections['/sensitivity/'+key+'/all_fields/plastic_slip']=value['all_fields']['plastic_slip']
    result=dict(base_commit=BASE,sources=sources,
                original_scope=c['scope'],original_C03=c['C03'],original_passed=c['passed'],
                limited_verified=c['limited_verified'],remaining=c['remaining'],counts=c['counts'],
                measurements_selections=selections,
                unchanged_thresholds=a['unchanged_thresholds'],
                threshold_source=paths[4]+'#/unchanged_thresholds',
                full_history_tolerance='未定義',accepted=None,
                unconfirmed=['original physical increment durations in seconds',
                             'original saved stress total/effective tag and component-order metadata in the cited records',
                             'a general MC/MCC coupled full-field tolerance',
                             'a universal material/contact partition dominance',
                             'future production output IDs/component metadata before supplied'],
                acceptance_filename_note='No acceptance.json in the C-03 supplement tree; its judgement is conclusion.json. Older endpoint gates are explicitly sourced separately.')
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2])
    p.add_argument('--out',type=Path,required=True)
    args=p.parse_args()
    extract(args.root,args.out)
