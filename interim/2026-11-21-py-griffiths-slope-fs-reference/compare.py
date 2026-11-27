"""Fail-closed interchange validator/comparator; no production solver is called."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import numpy as np

def case_fingerprint():
    value=json.loads((Path(__file__).resolve().parent/'case.json').read_text('utf-8-sig'))
    digest=hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
    if digest!='b67f2d11385204bca1bddff90a4f787cc8d3dce3eaf97277abcd81eacd24ce96':
        raise ValueError('pinned case.json changed; implementation supports this fixed case only')
    return digest

MODEL='griffiths_1999_coarse_reconstruction_D1'
UNITS={'length':'m','force':'kN','stress':'kPa','displacement':'m','strain':'1','factor':'1'}
CONVENTIONS={'x':'right','y':'up','stress':'tension_positive','displacement':'signed_global',
             'strain':'engineering_shear','kinematics':'plane_strain','pore_pressure':'dry_zero'}
REDUCTION={'cohesion':'c/F','friction':'atan(tan(phi)/F)','dilation':'fixed_zero'}
STATES={'equilibrated','certified_failure','unresolved','iteration_limit','trial_limit','timeout','cancelled','constitutive_unresolved'}
FORBIDDEN={'midpoint','fs_midpoint','point_value','fs_value','point_estimate','fs','factor_of_safety','safety_factor','factor_midpoint','estimated_fs'}


def number(v):
    if isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v):
        raise ValueError('finite numeric value required')
    return v


def reject_points(obj):
    if isinstance(obj,dict):
        if FORBIDDEN.intersection(obj): raise ValueError('Fs midpoint/point value forbidden')
        for v in obj.values(): reject_points(v)
    elif isinstance(obj,list):
        for v in obj: reject_points(v)


def validate(m):
    reject_points(m)
    expected={'schema_version','model_id','case_sha256','units','conventions','reduction','accepted','fs_tolerance',
              'trials','trial_count','factor_interval','endpoint_statuses','interval_kind','field_state'}
    if set(m)!=expected: raise ValueError('missing/unknown manifest field; point estimates are not interchangeable')
    if m.get('schema_version')!=1 or m.get('model_id')!=MODEL: raise ValueError('unknown model/schema')
    if m.get('units')!=UNITS or m.get('conventions')!=CONVENTIONS or m.get('reduction')!=REDUCTION:
        raise ValueError('units/sign/plane strain/reduction mismatch')
    if m.get('accepted','missing') is not None or m.get('fs_tolerance') is not None:
        raise ValueError('Fs acceptance tolerance is undefined: accepted must be null')
    case=m.get('case_sha256','')
    if case!=case_fingerprint(): raise ValueError('case fingerprint differs from pinned physical input')
    trials=m.get('trials')
    if not isinstance(trials,list) or not trials: raise ValueError('all explicit trials required')
    if m.get('trial_count')!=len(trials): raise ValueError('missing/implicitly omitted trials')
    seen={}
    for t in trials:
        f=number(t['factor']); state=t['status']
        if f<=0 or state not in STATES or f in seen: raise ValueError('invalid/duplicate trial')
        seen[f]=t
        if state=='equilibrated':
            if t.get('gravity_fraction')!=1.: raise ValueError('incomplete gravity trial')
            r=number(t['relative_residual']); tol=number(t['residual_tolerance'])
            if not 0<=r<=tol or not 0<tol<=1e-8: raise ValueError('equilibrated without residual evidence')
        elif state=='certified_failure':
            cert=t.get('failure_certificate')
            if not isinstance(cert,dict) or cert.get('independently_verified') is not True or not cert.get('artifact_sha256'):
                raise ValueError('nonconvergence is not a failure certificate')
            # Certificates require a verifier outside this diagnostic package.
            raise ValueError('failure-certificate verifier unimplemented; cannot certify Fs upper')
        if t.get('failure_certified',False): raise ValueError('unsupported failure certification')
    interval=m.get('factor_interval')
    if interval is not None:
        if not isinstance(interval,list) or len(interval)!=2: raise ValueError('two ordered endpoints required')
        lo,hi=map(number,interval)
        if not 0<lo<hi or lo not in seen or hi not in seen: raise ValueError('reversed/collapsed/untried interval')
        if seen[lo]['status']!='equilibrated': raise ValueError('lower endpoint is not equilibrated')
        if m.get('endpoint_statuses')!=[seen[lo]['status'],seen[hi]['status']]: raise ValueError('endpoint status mismatch')
        if m.get('interval_kind')!='computational_unresolved' or seen[hi]['status'] in {'equilibrated','certified_failure'}:
            raise ValueError('unresolved upper cannot be an excluded/physical Fs upper')
        if any(f>hi and t['status']=='equilibrated' for f,t in seen.items()):
            raise ValueError('nonmonotone trials: transition interval not justified')
    elif m.get('interval_kind') is not None or m.get('endpoint_statuses') is not None:
        raise ValueError('interval metadata without interval')
    fields=m.get('field_state')
    if fields is not None:
        f=number(fields['factor'])
        if f not in seen or seen[f]['status']!='equilibrated' or fields.get('gravity_fraction')!=1.:
            raise ValueError('fields require a complete equilibrated trial')
        if fields.get('reference')!='stress_free_gravity_turn_on' or fields.get('plastic_measure')!='engineering_tensor_components':
            raise ValueError('field reference/plastic definition mismatch')
    return seen


def table(path,columns):
    if Path(path).suffix=='.json':
        data=json.loads(Path(path).read_text('utf-8-sig'))
        if data.get('columns')!=columns or not isinstance(data.get('rows'),list):
            raise ValueError('JSON field columns/units contract mismatch')
        rows=[dict(zip(columns,row,strict=True)) for row in data['rows']]
    else:
        with Path(path).open(encoding='utf-8-sig',newline='') as f:
            reader=csv.DictReader(f)
            if reader.fieldnames!=columns: raise ValueError('CSV columns/units/sign contract mismatch')
            rows=list(reader)
    if not rows: raise ValueError('empty field')
    if any(set(row)!=set(columns) for row in rows): raise ValueError('malformed field row')
    arr=np.array([[float(row[k]) if isinstance(row[k],str) else number(row[k]) for k in columns] for row in rows])
    if not np.isfinite(arr).all() or len(set(arr[:,0]))!=len(arr): raise ValueError('nonfinite/duplicate field IDs')
    if np.any(arr[:,0]!=arr[:,0].astype(int)) or np.any(arr[:,0]<1): raise ValueError('invalid ID')
    return arr[np.argsort(arr[:,0])]

NODE=['node_id','x_m','y_m','ux_m','uy_m']
PLASTIC=['element_id','x_m','y_m','plastic_exx','plastic_eyy','plastic_ezz','plastic_gamma_xy']


def compare(candidate,reference,candidate_dir=Path('.'),reference_dir=Path('.')):
    validate(candidate); validate(reference)
    if candidate['case_sha256']!=reference['case_sha256']: raise ValueError('physical inputs differ')
    result={'accepted':None,'fs_tolerance':'未定義','fs_relation':'判定不能',
            'reason':'computational unresolved intervals are not certified Fs enclosures'}
    ci=candidate.get('factor_interval'); ri=reference.get('factor_interval')
    if ci and ri:
        result['numeric_computational_relation']=('disjoint' if ci[1]<ri[0] or ri[1]<ci[0] else
            'contains_reference' if ci[0]<=ri[0] and ci[1]>=ri[1] else 'overlap')
    a=candidate.get('field_state'); b=reference.get('field_state')
    if (a is None)!=(b is None): raise ValueError('missing paired fields')
    if a:
        if a['factor']!=b['factor'] or a['gravity_fraction']!=b['gravity_fraction']:
            raise ValueError('different factor/loading states; field comparison refused')
        for key,cols in [('nodes_csv',NODE),('plastic_csv',PLASTIC)]:
            x=table(candidate_dir/a[key],cols); y=table(reference_dir/b[key],cols)
            if x.shape!=y.shape or not np.array_equal(x[:,:3],y[:,:3]):
                raise ValueError('common IDs/physical points required; no silent interpolation/h mapping')
            numerator=float(np.linalg.norm(x[:,3:]-y[:,3:])); denominator=float(np.linalg.norm(y[:,3:]))
            value=numerator/denominator if denominator>0 else None
            result[key]={'relative_L2':value,'absolute_L2':numerator,'accepted':None}
            if key=='plastic_csv':
                result[key]['plastic_zone_rule']='any engineering plastic strain component != 0; diagnostic only'
                result[key]['plastic_zone_disagreements']=int(np.sum(np.any(x[:,3:]!=0,axis=1)!=np.any(y[:,3:]!=0,axis=1)))
            # Original 0.05 gate belongs to ORIGINAL 661-node displacement sample;
            # applying it to a new CST grid would create a new acceptance criterion.
    return result


def main():
    p=argparse.ArgumentParser(); p.add_argument('candidate',type=Path); p.add_argument('--reference',type=Path)
    p.add_argument('--lem-results',type=Path)
    p.add_argument('--output',type=Path,required=True); args=p.parse_args()
    m=json.loads(args.candidate.read_text('utf-8')); validate(m)
    r={'valid_interchange':True,'accepted':None}
    # Validate referenced fields even in single-manifest mode.
    f=m.get('field_state')
    if f:
        table(args.candidate.parent/f['nodes_csv'],NODE); table(args.candidate.parent/f['plastic_csv'],PLASTIC)
    if args.reference:
        r=compare(m,json.loads(args.reference.read_text('utf-8')),args.candidate.parent,args.reference.parent)
    if args.lem_results:
        env=json.loads((args.lem_results.parent/'environment.json').read_text('utf-8'))
        if env['case_sha256']!=m['case_sha256']: raise ValueError('LEM physical input fingerprint differs')
        interval=m.get('factor_interval')
        r['lem_numeric_diagnostics']=[]
        for item in json.loads(args.lem_results.read_text('utf-8')):
            for method,v in item['minima'].items():
                if v is not None:
                    relation=None if interval is None else ('含まれる' if interval[0]<=v['fs']<=interval[1] else '外れる')
                    r['lem_numeric_diagnostics'].append(dict(grid=item['grid'],slices=item['slices'],noncircular=item['noncircular'],
                        method=method,fs=v['fs'],numeric_computational_relation=relation,physical_judgment='判定不能',compression_only=v['compression_only'],accepted=None))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(r,ensure_ascii=False,indent=2,allow_nan=False)+'\n','utf-8')

if __name__=='__main__': main()
