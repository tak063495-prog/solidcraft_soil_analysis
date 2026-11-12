"""P6-03 evidence/state audit; never runs a solver or changes evidence."""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
BASE='e3d830d8d1bdffd727c4ab7271755f52e556006d'
EXPECTED={'mcc-spatial':'未達','extrusion-stress':'判定不能','oscillating-monitor':'判定不能',
          'finite-depth-stress':'判定不能','pile':'未達','rainfall':'未達'}
RAIN=[(.0046205475457005325,.003,'m'),(1.348244318335049,.8,'1'),(.998450823613834,.8,'1'),
      (.0038863726982864173,.003,'m'),(.00014219441327791016,.0001,'m3')]

def sha(data):return hashlib.sha256(data).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def exact(a,b):return json.dumps(a,sort_keys=True,ensure_ascii=False)==json.dumps(b,sort_keys=True,ensure_ascii=False)
def load_old():
    path=ROOT/'interim/2026-12-06-p6-01-status-reconciliation/audit_status.py'
    spec=importlib.util.spec_from_file_location('p603_previous_audit',path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod

def audit(status=None, *, verify_files=True):
    s=read(HERE/'status.json') if status is None else status
    errors=[]
    def check(ok,code):
        if not ok:errors.append(code)
    check(s.get('base_commit')==BASE,'BASE')
    check(s.get('target_id')=='TODO:P6-03' and s.get('state')=='未達','TARGET_STATE')
    for key in ('p6_03_complete','p6_complete','nc01_complete'):
        check(s.get(key) is False,'COMPLETION:'+key)
    rows=s.get('rows',[])
    check(len(rows)==6 and {r.get('id') for r in rows}==set(EXPECTED),'COVERAGE')
    oldledger=read(ROOT/'interim/2026-12-06-p6-01-status-reconciliation/ledger.json')
    oldrows={r['id']:r for r in oldledger['rows']}
    check(oldrows['TODO:P6-03:overall']['state']=='未達','P601_OVERALL')
    check(oldrows['TODO:P6-03:indeterminate']['state']=='判定不能','P601_INDETERMINATE')
    for r in rows:
        rid=r.get('id','?')
        check(r.get('state')==EXPECTED.get(rid),'ROW_STATE:'+rid)
        check(r.get('independent_state')=='判定不能','INDEPENDENT_STATE:'+rid)
        check(r.get('completed') is False and r.get('new_overall_acceptance') is False,'ROW_COMPLETION:'+rid)
        check(r.get('pc_state')=='PC手順待ち' and bool(r.get('pc')),'PC:'+rid)
        for key in ('bounded','observed','independent','missing','sources'):
            check(bool(r.get(key)),'REQUIRED:'+rid+':'+key)
        for name in r.get('ledger',[]):
            check(name in oldrows,'P601_ROW:'+name)
            if name in oldrows:
                check(oldrows[name]['completed'] is False,'P601_COMPLETED:'+name)
        check('TODO:P6-03:overall' in r.get('ledger',[]),'P601_LINK:'+rid)
        if r.get('state')=='判定不能':
            check('TODO:P6-03:indeterminate' in r.get('ledger',[]),'P601_SCOPE:'+rid)
    mc=s.get('original_mc_column',{})
    check(exact(mc,dict(material='nonassociated_MC',tet4_fine_percent=6.828,tet4_fine_accepted=False,
                   continuum_estimate_percent=1.034,continuum_is_separate=True,estimate_accepts_tet4=False)),'MC_SEPARATION')
    rain=s.get('old_rainfall_five',[])
    check(len(rain)==5,'RAIN_COUNT')
    for i,(v,t,u) in enumerate(RAIN):
        expected=dict(id=i+1,value=v,upper_bound=t,operator='<',passed=False,unit=u)
        check(i<len(rain) and exact(rain[i],expected),'RAIN_ORIGINAL:'+str(i+1))
        check(not (v<t),'RAIN_FAILURE:'+str(i+1))
    check(exact(s.get('reference_legacy_rows'),read(HERE/'references/pr33/results/legacy-rows.json')['rows']),'REFERENCE_LEGACY_ROWS')
    pile=read(HERE/'references/pr28/results/historical-diagnosis.json')
    check(exact(s.get('pile_metrics'),{k:pile[k] for k in ('historical_tet10_counts','historical_settlement_mm','historical_relative_difference','metric_is_exact_solution_error','spatial_threshold','accepted','historical_spatial_convergence_accepted')}),'PILE_METRICS')
    check(exact(s.get('old_rainfall_source_rows'),read(HERE/'references/pr32/original-five.json')['old_rows']),'RAIN_REFERENCE_ROWS')
    links=read(HERE/'p601-row-links.json')
    for row in links['rows']:
        check(exact(row,oldrows.get(row['id'])),'P601_EXACT_ROW:'+row['id'])
    check(s.get('new_independent_comparison_count')==1,'COMPARISON_BUDGET')
    check(s.get('historical_python')==dict(p601_discovered=625,p601_passed=617,p601_unexecuted=8,p601_setup_errors=2,retained=True),'HISTORICAL_COUNTS')
    check(s.get('windows_ci',{}).get('success_claimed') is False,'CI_NOT_CONFIRMED')
    if verify_files:
        tails=read(HERE/'document-appendices.json')
        for item in read(HERE/'branch-sources.json')['files']:
            p=item['path']
            frozen=subprocess.check_output(['git','-C',str(ROOT),'show',BASE+':'+p])
            check(sha(frozen)==item['sha256'],'BASE_SHA:'+p)
            blob=hashlib.sha1(b'blob '+str(len(frozen)).encode()+b'\0'+frozen).hexdigest()
            check(blob==item['git_blob_sha1'],'BASE_BLOB:'+p)
            if p.endswith('/audit_status.py'):
                continue # Explicit strict appendix compatibility change, audited by old tests.
            expected=frozen+tails.get(p,'').encode()
            try:actual=(ROOT/p).read_bytes()
            except OSError:actual=None
            check(actual==expected,'CURRENT_BYTES:'+p)
        refs=read(HERE/'reference-sources.json')
        check({r['pr'] for r in refs['references']}=={24,33,28,32},'REFERENCE_COVERAGE')
        for item in refs['files']:
            data=(ROOT/item['local']).read_bytes()
            check(sha(data)==item['sha256'],'REFERENCE_SHA:'+item['local'])
            check(hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==item['git_blob_sha1'],'REFERENCE_BLOB:'+item['local'])
            check(item.get('modified') is False,'REFERENCE_MODIFIED:'+item['local'])
        old=load_old().audit(ROOT,oldledger)
        check(old['integrity_ok'],'P601_INTEGRITY')
        # A small independent comparison is a separate, narrower result.
        comp=read(HERE/'validation/mc-crest/comparison.json')
        check(comp['complete'] and comp['criterion_met'] and comp['comparison_count']==1,'MC_COMPARISON')
        check(comp['new_all_field_acceptance'] is False and comp['p6_03_complete'] is False,'MC_SCOPE')
        check(comp['tet4_original_failure_preserved'] is True,'MC_OLD_FAILURE')
        for r in comp['rows']:
            check(r['original_displacement_limit_m']==5e-8 and r['operator']=='<' and
                  r['absolute_error_m']<5e-8,'MC_ORIGINAL_LIMIT')
    return dict(schema='p603-audit/v1',integrity_ok=not errors,errors=errors,row_count=len(rows),
                numerical_overall_acceptance=False,p6_03_complete=False)

def table(status):
    lines=['|対象|全体/判定の状態|既存の限定受入範囲|独立参照との対応・不足|PC手順待ち|',
           '|---|---|---|---|---|']
    for r in status['rows']:
        vals=[r['target'],r['state'],r['bounded']+' '+r['observed'],r['independent']+' '+r['missing'],r['pc']]
        lines.append('|'+ '|'.join(v.replace('|','/').replace('\n',' ') for v in vals)+'|')
    return '\n'.join(lines)+'\n'

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--table',action='store_true');a=p.parse_args()
    if a.table:print(table(read(HERE/'status.json')),end='')
    else:
        result=audit();print(json.dumps(result,ensure_ascii=False,indent=2));raise SystemExit(0 if result['integrity_ok'] else 1)
