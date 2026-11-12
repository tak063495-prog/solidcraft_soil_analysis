"""Export the reduced column reference as explicit JSON/CSV consumer examples.

This is an example consumer contract only; it does not add a production schema.
The normalized quasistatic cycle is explicitly assigned a duration of 24 s.
"""
import argparse
import csv
import json
from pathlib import Path
from compare import META, BULK, CONTACT, load


def voigt(a):
    return [a[0][0],a[1][1],a[2][2],a[0][1],a[1][2],a[2][0]]


def export(history,out):
    out.mkdir(parents=True,exist_ok=True)
    metadata=dict(META,stress_kind='effective',biot_alpha=1.)
    schedule=[dict(stage=h['stage'],increment=h['increment'],time=24.*h['time']) for h in history]
    manifest=dict(schema='py-c03-history/1',metadata=metadata,schedule=schedule,bulk_ids=['column/B0'],contact_ids=['interface/C0'])
    records=[]
    for h,s in zip(history,schedule):
        records.append(dict(s,kind='bulk',id='column/B0',stress=voigt(h['stress']),plastic_strain=voigt(h['plastic_strain']),pore=h['pore']))
        records.append(dict(s,kind='contact',id='interface/C0',traction=[h['traction'],0.,-h['effective_contact_pressure']],plastic_slip=[h['plastic_slip'],0.,0.]))
    envelope=dict(schema=manifest['schema'],metadata=metadata,records=records)
    for name,obj in [('manifest.json',manifest),('reference.json',envelope),('reference.csv.meta.json',dict(schema=manifest['schema'],metadata=metadata))]:
        (out/name).write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
    cols=['kind','stage','increment','time','id','pore']+[f'{f}_{i}' for f,size in (BULK|CONTACT).items() if size>1 for i in range(size)]
    with (out/'reference.csv').open('w',encoding='utf-8',newline='') as file:
        w=csv.DictWriter(file,fieldnames=cols,lineterminator='\n')
        w.writeheader()
        for r in records:
            row={k:r[k] for k in ('kind','stage','increment','time','id')}
            for field,size in (BULK if r['kind']=='bulk' else CONTACT).items():
                if size==1:row[field]=r[field]
                else:
                    for i in range(size):row[f'{field}_{i}']=r[field][i]
            w.writerow(row)
    if load(out/'reference.json',manifest)!=load(out/'reference.csv',manifest):
        raise ValueError('JSON/CSV example roundtrip mismatch')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--history',type=Path,required=True)
    p.add_argument('--case',default='m8-c8')
    p.add_argument('--out',type=Path,required=True)
    args=p.parse_args()
    export(json.loads(args.history.read_text(encoding='utf-8'))[args.case],args.out)
