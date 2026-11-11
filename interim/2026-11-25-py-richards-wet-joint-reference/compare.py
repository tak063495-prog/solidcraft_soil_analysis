"""Fail-closed, proposed exchange contract, not a production schema/exporter."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
import numpy as np
from reference import PARAMETERS,write

CONVENTIONS=dict(contract='py-richards-joint-diagnostic/1-proposed',
    units=dict(head='m',flow='m3/s',time='s',pressure='kPa',gap='m',water='m3',work='kPa*m'),
    signs=dict(head='total_H=psi+z',flow_a='a_to_joint_positive',flow_b='joint_to_b_positive',
               pressure='compression_positive',gap='opening_positive',shear='positive_with_imposed_displacement'),
    time_basis='physical_seconds_from_initial_equilibrium',sampling='instantaneous_backward_euler_endpoint',
    interpolation='none',id_mapping='explicit_identity',generation_policy='fixed_generation_zero',
    pressure_policy='max(gamma_w*(H-z),0);no_joint_suction',complete=True)
CATEGORICAL=('joint_id','point_id','generation','hydraulic_open','contact_state')
KEYS=('time_s','joint_id','point_id','generation')
NUMERIC=('head_a_m','head_b_m','joint_head_m','q_a_to_joint_m3_s','q_joint_to_b_m3_s',
         'reservoir_inflow_m3_s','joint_storage_rate_m3_s','gap_m','hydraulic_aperture_m','physical_pressure_kpa','effective_pressure_kpa','shear_kpa',
         'plastic_slip_m','friction_dissipation_kpa_m','water_balance_defect_m3')
FIELDS=set(KEYS+CATEGORICAL+NUMERIC)


def packet(result,producer='independent_reference'):
    config=result['configuration']
    digest=hashlib.sha256(json.dumps(dict(configuration=config,parameters=PARAMETERS),sort_keys=True).encode()).hexdigest()
    return dict(metadata=dict(CONVENTIONS,configuration=config,parameters=PARAMETERS,input_sha256=digest,
               producer=producer,times_s=sorted(set(r['time_s'] for r in result['rows'])),
               point_ids=list(range(config['ny'])),joint_ids=['J0']),rows=result['rows'])


def read(path):
    path=Path(path)
    if path.suffix.lower()=='.json':return json.loads(path.read_text(encoding='utf-8'))
    if path.suffix.lower()!='.csv':raise ValueError('CSV or JSON required')
    metadata=json.loads(path.with_suffix('.metadata.json').read_text(encoding='utf-8'))
    with path.open(encoding='utf-8',newline='') as f:
        reader=csv.DictReader(f)
        if len(reader.fieldnames)!=len(set(reader.fieldnames)):raise ValueError('duplicate CSV column')
        rows=list(reader)
    for row in rows:
        for k in ('point_id','generation'):row[k]=int(row[k])
        if row['hydraulic_open'] not in ('true','false'):raise ValueError('explicit boolean required')
        row['hydraulic_open']=row['hydraulic_open']=='true'
        for k in ('time_s',)+NUMERIC:row[k]=float(row[k])
    return dict(metadata=metadata,rows=rows)


def validate(p):
    if set(p)!=set(('metadata','rows')):raise ValueError('packet keys differ')
    m=p['metadata']
    if set(m)!=set(CONVENTIONS)|set(('configuration','parameters','input_sha256','producer','times_s','point_ids','joint_ids')):
        raise ValueError('metadata fields differ')
    for k,v in CONVENTIONS.items():
        if m[k]!=v:raise ValueError('convention differs: '+k)
    if m['producer'] not in ('independent_reference','future_body_export','synthetic_contract_check'):
        raise ValueError('explicit producer required')
    if m['parameters']!=PARAMETERS:raise ValueError('input parameters differ')
    digest=hashlib.sha256(json.dumps(dict(configuration=m['configuration'],parameters=m['parameters']),sort_keys=True).encode()).hexdigest()
    if m['input_sha256']!=digest:raise ValueError('input digest differs')
    config=m['configuration']
    if config['mode'] not in ('ONE_WAY','TWO_WAY') or config['policy'] not in ('CONTINUOUS','PRESCRIBED'):
        raise ValueError('unknown coupling/mapping')
    ts=m['times_s'];ids=m['point_ids']
    if not ts or ts!=sorted(set(ts)) or any(type(t) not in (int,float) or not math.isfinite(t) or t<=0 for t in ts):
        raise ValueError('strictly increasing finite physical times required')
    if ts[-1]!=config['end_s'] or m['joint_ids']!=['J0'] or ids!=list(range(config['ny'])):
        raise ValueError('time/ID inventory differs')
    expected=[(t,'J0',i,0) for t in ts for i in ids]
    actual=[]
    for row in p['rows']:
        if set(row)!=FIELDS:raise ValueError('row fields missing or extra')
        if type(row['point_id']) is not int or type(row['generation']) is not int or type(row['hydraulic_open']) is not bool:
            raise ValueError('typed ID/generation/open state required')
        if row['contact_state'] not in ('open','stick','slip'):raise ValueError('explicit contact state required')
        for k in ('time_s',)+NUMERIC:
            if type(row[k]) not in (int,float) or not math.isfinite(row[k]):raise ValueError('finite numeric field required: '+k)
        if min(row[k] for k in ('gap_m','hydraulic_aperture_m','physical_pressure_kpa','effective_pressure_kpa','friction_dissipation_kpa_m'))<0:
            raise ValueError('pressure/aperture sign violated')
        ledger=row['q_a_to_joint_m3_s']-row['q_joint_to_b_m3_s']+row['reservoir_inflow_m3_s']-row['joint_storage_rate_m3_s']
        scale=sum(abs(row[k]) for k in ('q_a_to_joint_m3_s','q_joint_to_b_m3_s','reservoir_inflow_m3_s','joint_storage_rate_m3_s'))
        # Existing C-04 hydraulic flux criteria, not a new coupled acceptance gate.
        if abs(ledger)>1e-16+1e-8*scale:raise ValueError('oriented joint water ledger differs')
        if row['hydraulic_open']!=(row['hydraulic_aperture_m']>0):raise ValueError('hydraulic opening differs')
        if config['mode']=='TWO_WAY' and row['hydraulic_aperture_m']!=row['gap_m']:
            raise ValueError('TWO_WAY aperture must be the mechanical gap')
        if row['effective_pressure_kpa']>0 and row['contact_state']=='open':raise ValueError('open contact carries pressure')
        actual.append(tuple(row[k] for k in KEYS))
    if actual!=expected:raise ValueError('missing/duplicate/reordered ID/time/generation; no interpolation')
    return p


def compare(reference,candidate):
    validate(reference);validate(candidate)
    for k in reference['metadata']:
        if k!='producer' and reference['metadata'][k]!=candidate['metadata'][k]:
            raise ValueError('metadata mismatch: '+k)
    deviations=[]
    for a,b in zip(reference['rows'],candidate['rows'],strict=True):
        for k in CATEGORICAL:
            if a[k]!=b[k]:raise ValueError('categorical state mismatch: '+k)
        deviations.append(dict(key=[a[k] for k in KEYS],difference={k:b[k]-a[k] for k in NUMERIC}))
    metrics={k:dict(rms=float(np.sqrt(np.mean([d['difference'][k]**2 for d in deviations]))),
                    maximum_absolute=max(abs(d['difference'][k]) for d in deviations)) for k in NUMERIC}
    return dict(accepted=None,tolerance='未定義',A04='未達',C04='未達',
                actual_body_comparison=candidate['metadata']['producer']=='future_body_export',
                rows=len(deviations),metrics=metrics,deviations=deviations)


def main():
    p=argparse.ArgumentParser();p.add_argument('reference',type=Path);p.add_argument('candidate',type=Path);p.add_argument('report',type=Path)
    a=p.parse_args()
    try:report=compare(read(a.reference),read(a.candidate))
    except (ValueError,KeyError,TypeError) as error:
        write(a.report,dict(accepted=None,valid=False,error=str(error)));return 1
    write(a.report,report);return 2

if __name__=='__main__':sys.exit(main())
