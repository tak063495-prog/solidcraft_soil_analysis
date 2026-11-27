"""Strict future-output consumer; no production schema changes or interpolation.

Usage: compare.py --manifest manifest.json reference.json candidate.json --out diff.json
CSV requires an adjacent filename.csv.meta.json envelope with metadata only.
"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np


META = dict(stress_unit='kPa',length_unit='m',time_unit='s',
            stress_sign='tension_positive',strain_sign='tension_positive',
            pore_sign='positive_compression',strain_shear='tensor',
            pore_model='fixed_hydrostatic_nonnegative',
            traction_kind='effective_skeleton',
            stress_components=['xx','yy','zz','xy','yz','zx'],
            vector_components=['x','y','z'],
            traction_definition='force_on_a_per_area_with_jump_ub_minus_ua',
            plastic_slip_definition='ub_minus_ua_plastic_part')
BULK = {'stress':6,'plastic_strain':6,'pore':1}
CONTACT = {'traction':3,'plastic_slip':3}


def no_duplicates(pairs):
    result={}
    for k,v in pairs:
        if k in result:
            raise ValueError('duplicate JSON key: '+k)
        result[k]=v
    return result


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'),
                      object_pairs_hook=no_duplicates,
                      parse_constant=lambda x: (_ for _ in ()).throw(ValueError('nonfinite '+x)))


def finite(v):
    if isinstance(v,bool) or not isinstance(v,(int,float)) or not np.isfinite(v):
        raise ValueError('finite numeric value required')
    return float(v)


def validate_metadata(metadata):
    if set(metadata)!=set(META)|{'stress_kind','biot_alpha'}:
        raise ValueError('metadata fields mismatch')
    for key,value in META.items():
        if metadata[key]!=value:
            raise ValueError('unit/sign/components convention mismatch: '+key)
    if metadata['stress_kind'] not in ('effective','total'):
        raise ValueError('explicit total/effective stress kind required')
    if not 0 <= finite(metadata['biot_alpha']) <= 1:
        raise ValueError('invalid Biot alpha')


def validate_manifest(m):
    if set(m)!={'schema','metadata','schedule','bulk_ids','contact_ids'} or m['schema']!='py-c03-history/1':
        raise ValueError('manifest schema mismatch')
    validate_metadata(m['metadata'])
    keys=[]
    closed=set()
    previous_stage=None
    previous_time=-np.inf
    next_increment=1
    for row in m['schedule']:
        if set(row)!={'stage','increment','time'}:
            raise ValueError('schedule fields mismatch')
        stage,inc=row['stage'],row['increment']
        if not isinstance(stage,str) or not stage:
            raise ValueError('stage ID required')
        if stage!=previous_stage:
            if stage in closed:
                raise ValueError('noncontiguous stage')
            closed.add(stage)
            previous_stage=stage
            next_increment=1
        if type(inc) is not int or inc!=next_increment:
            raise ValueError('missing/noninteger increment in manifest')
        next_increment+=1
        t=finite(row['time'])
        if t<=previous_time:
            raise ValueError('nonincreasing physical time')
        previous_time=t
        keys.append((stage,inc,t))
    if not keys:
        raise ValueError('empty schedule')
    for kind in ('bulk','contact'):
        ids=m[kind+'_ids']
        if not isinstance(ids,list) or not ids or any(not isinstance(x,str) or not x for x in ids) or len(ids)!=len(set(ids)):
            raise ValueError('nonempty unique string GP IDs required')
    return keys


def load(path,m):
    schedule=validate_manifest(m)
    path=Path(path)
    if path.suffix.lower()=='.csv':
        envelope=read_json(str(path)+'.meta.json')
        validate_metadata(envelope['metadata'])
        if envelope!={'schema':m['schema'],'metadata':m['metadata']}:
            raise ValueError('CSV sidecar metadata mismatch')
        rows=[]
        with path.open(encoding='utf-8-sig',newline='') as f:
            reader=csv.DictReader(f)
            if reader.fieldnames is None or len(reader.fieldnames)!=len(set(reader.fieldnames)):
                raise ValueError('duplicate/missing CSV headers')
            expected={'kind','stage','increment','time','id','pore'}|{f'{field}_{i}' for field,size in (BULK|CONTACT).items() if size>1 for i in range(size)}
            if set(reader.fieldnames)!=expected:
                raise ValueError('CSV headers mismatch')
            for row in reader:
                kind=row['kind']
                fields=BULK if kind=='bulk' else CONTACT if kind=='contact' else None
                if fields is None:
                    raise ValueError('unknown CSV kind')
                out=dict(kind=kind,stage=row['stage'],increment=int(row['increment']),time=float(row['time']),id=row['id'])
                used={'kind','stage','increment','time','id'}
                for field,size in fields.items():
                    cols=[field] if size==1 else [f'{field}_{i}' for i in range(size)]
                    vals=[float(row[c]) for c in cols]
                    used.update(cols)
                    out[field]=vals[0] if size==1 else vals
                if any(row[c] not in ('',None) for c in set(row)-used):
                    raise ValueError('unexpected values for CSV record kind')
                rows.append(out)
    else:
        envelope=read_json(path)
        validate_metadata(envelope['metadata'])
        if set(envelope)!={'schema','metadata','records'} or envelope['schema']!=m['schema'] or envelope['metadata']!=m['metadata']:
            raise ValueError('JSON metadata/schema mismatch')
        rows=envelope['records']
    values={}
    for r in rows:
        kind=r.get('kind')
        fields=BULK if kind=='bulk' else CONTACT if kind=='contact' else None
        if fields is None or set(r)!=set(fields)|{'kind','stage','increment','time','id'}:
            raise ValueError('record fields/kind mismatch')
        if type(r['increment']) is not int:
            raise ValueError('integer increment required')
        key=(r['stage'],r['increment'],finite(r['time']),kind,r['id'])
        if key in values:
            raise ValueError('duplicate GP/increment')
        normalized={}
        for name,size in fields.items():
            v=r[name]
            if size==1:
                normalized[name]=[finite(v)]
                if name=='pore' and normalized[name][0]<0:
                    raise ValueError('negative pore pressure violates fixed hydrostatic convention')
            elif not isinstance(v,list) or len(v)!=size:
                raise ValueError('wrong component count: '+name)
            else:
                normalized[name]=[finite(x) for x in v]
        values[key]=normalized
    expected={(*s,kind,ident) for s in schedule for kind in ('bulk','contact') for ident in m[kind+'_ids']}
    if set(values)!=expected:
        raise ValueError(f'missing/extra ID or increment or exact time: missing={len(expected-set(values))}, extra={len(set(values)-expected)}; interpolation forbidden')
    return values


def compare(manifest,reference,candidate):
    m=read_json(manifest)
    a,b=load(reference,m),load(candidate,m)
    result=dict(accepted=None,tolerance='未定義',C03='未達',records=len(a),fields={})
    for kind,fields in (('bulk',BULK),('contact',CONTACT)):
        for name in fields:
            keys=sorted(k for k in a if k[3]==kind)
            av=np.array([a[k][name] for k in keys])
            bv=np.array([b[k][name] for k in keys])
            d=av-bv
            denom=np.linalg.norm(av)
            result['fields'][kind+'.'+name]=dict(max_abs=float(abs(d).max()),l2_absolute=float(np.linalg.norm(d)),
                                                l2_relative=None if denom==0 else float(np.linalg.norm(d)/denom),accepted=None)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest',required=True)
    p.add_argument('reference')
    p.add_argument('candidate')
    p.add_argument('--out',type=Path,required=True)
    args=p.parse_args()
    result=compare(args.manifest,args.reference,args.candidate)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
    print('比較完了; 許容値未定義; accepted=null; C-03未達')
