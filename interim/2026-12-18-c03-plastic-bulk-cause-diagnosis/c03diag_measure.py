"""All-history field/GP mapping and factorial attribution, never acceptance repair."""
import argparse
from collections import Counter
from pathlib import Path
import json
import numpy as np
from c03diag_metrics import load,compare,metric,norm,cosine,budget
from c03diag_study import read,write,Model


def main():
    p=argparse.ArgumentParser();p.add_argument('--matrix',type=Path,required=True)
    p.add_argument('--scope',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();base=a.matrix;specs=read(base/'specifications.json');data={}
    for spec in specs:
        folder=base/spec['name']
        if (folder/'fields.jsonl').exists() and ((folder/'result.json').exists() or (folder/'schedule.json').exists()):
            data[spec['name']]=load(folder/'fields.jsonl')
    report=dict(schema='c03-bulk-cause/1',passed=False,comparisons={},sensitivity={},interactions={},
                partition_histograms={},gp_mapping={},mcc={},reference_assessed_as_exact=False)
    maps=[]
    for name,(rows,fields) in sorted(data.items()):
        spec=next(s for s in specs if s['name']==name)
        assert len(rows)==7*spec['steps']
        assert [(r['stage'],r['progress']) for r in rows]==[(s,i/spec['steps']) for s in ('geostatic','shear','open','recontact','reverse','reload','return') for i in range(1,spec['steps']+1)]
        if name.endswith('-native'):
            run=read(base/name/'raw/run.json');maps.append((name,run['integration_points'],run['joint_integration_points']))
            records=[json.loads(line) for line in (base/name/'raw/constitutive.jsonl').open()]
            mp=Counter(c for r in records for c in r['partitions'][:12]);cp=Counter(c for r in records for c in r['partitions'][12:])
        else:
            records=read(base/name/'result.json')['rows']
            mp=Counter(c for r in records for c in r['material_partitions']);cp=Counter(c for r in records for c in r['contact_partitions'])
        report['partition_histograms'][name]=dict(material=dict(sorted(mp.items())),contact=dict(sorted(cp.items())))
    first=maps[0]
    assert all(x[1:]==first[1:] for x in maps)
    model=Model(read(base/'mesh.json'),read(base/'ground.json'),read(base/'m0c0-native/stages.json'))
    material_error=float(np.max(np.abs(np.asarray([p['position_m'] for p in first[1]])-model.positions)))
    contact_error=float(np.max(np.abs(np.asarray([p['position_m'] for p in first[2]])-np.asarray([p['position_m'] for p in model.system.points]))))
    assert material_error<=1e-12 and contact_error<=1e-12
    assert all(p['joint_id']==q['joint_id'] and p['face_index']==q['face_index'] for p,q in zip(first[2],model.system.points))
    contact_checks={}
    for key in ('normal','shape','weight_m2'):
        contact_checks[key]=float(np.max(np.abs(np.asarray([p[key] for p in first[2]])-np.asarray([p[key] for p in model.system.points]))))
        assert contact_checks[key]<=1e-12
    for key in ('nodes_a','nodes_b','owner_a','owner_b'):
        assert [p[key] for p in first[2]]==[p[key] for p in model.system.points]
    weight_error=float(np.max(np.abs(np.asarray([p['weight_m3'] for p in first[1]])-model.weights)))
    assert weight_error<=1e-12
    # Exact ID order is fixed by element/local_index and face_index plus within-face order.
    report['gp_mapping']=dict(native_maps_exact_all_routes=True,material_ids=[(p['element_index'],p['local_index']) for p in first[1]],
          contact_ids=[(p['joint_id'],p['face_index'],i) for i,p in enumerate(first[2])],
          native_descriptors=dict(material=first[1],contact=first[2]),independent_material_position_maximum_difference=material_error,
          independent_contact_position_maximum_difference=contact_error,independent_contact_descriptor_maximum_differences=contact_checks,independent_material_weight_maximum_difference=weight_error)
    for name in sorted(data):
        if name.endswith('-native'):
            stem=name[:-7];other=stem+'-independent'
            if other in data:report['comparisons'][stem]=compare(data[name],data[other])
    for kind in ('native','independent'):
        for m in range(3):
            for c in range(3):
                x=f'm{m}c{c}-{kind}';y=f'm{m}c{c+1}-{kind}'
                if x in data and y in data:report['sensitivity'][x+'--'+y]=compare(data[x],data[y])
        for c in range(3):
            for m in range(2):
                x=f'm{m}c{c}-{kind}';y=f'm{m+1}c{c}-{kind}'
                report['sensitivity'][x+'--'+y]=compare(data[x],data[y])
        for m in (0,2):
            for field in ('u','stress','plastic_strain','plastic_slip','traction','plastic_dissipation'):
                deltas=[data[f'm{m}c{c}-{kind}'][1][field]-data[f'm{m}c{c+1}-{kind}'][1][field] for c in range(3)]
                report['interactions'][f'{kind}-m{m}-{field}-directions']=dict(adjacent_absolute_l2=[norm(d) for d in deltas],
                    adjacent_cosines=[cosine(deltas[i],deltas[i+1]) for i in range(2)])
        for field in ('u','stress','plastic_strain','plastic_slip','traction'):
            d0=data[f'm0c1-{kind}'][1][field]-data[f'm0c2-{kind}'][1][field]
            d2=data[f'm2c1-{kind}'][1][field]-data[f'm2c2-{kind}'][1][field]
            report['interactions'][f'{kind}-material-contact-{field}']=budget(d0,d2)
        fine=data[f'g2-m2c2-{kind}'];coarse=data[f'm2c2-{kind}']
        sampled=([r for i,r in enumerate(fine[0]) if i%2==1],{k:v[1::2] for k,v in fine[1].items()})
        report['sensitivity'][f'{kind}-global32--64-m2c2']=compare(coarse,sampled)
    for m in (8,32):
        x=f'fixed-m{m}c4-independent';y=f'fixed-m{m}c16-independent'
        report['sensitivity'][x+'--'+y]=compare(data[x],data[y])
    for c in (4,16):
        x=f'fixed-m8c{c}-independent';y=f'fixed-m32c{c}-independent'
        report['sensitivity'][x+'--'+y]=compare(data[x],data[y])
    for label in ('mcc-pc200','mcc-pc105','mcc-pc105-open1'):
        folder=a.scope/label
        report['mcc'][label]=compare(load(folder/'native-fields.jsonl'),load(folder/'independent/fields.jsonl'))
        fields=load(folder/'native-fields.jsonl')[1]
        report['mcc'][label]['native_maximum_plastic_strain']=float(np.max(np.abs(fields['plastic_strain'])))
        report['mcc'][label]['native_contact_status_counts']=dict(Counter(map(int,fields['contact_status'].ravel())))
        report['mcc'][label]['native_open_gp_states']=int(fields['open'].sum())
        report['mcc'][label]['native_pc_range']=[float(fields['pc'].min()),float(fields['pc'].max())]
    write(a.output,report)

if __name__=='__main__':main()
