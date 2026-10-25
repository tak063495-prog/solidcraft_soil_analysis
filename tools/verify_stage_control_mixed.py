"""Verify unchanged mixed-column histories and independently check the new strain diagnostic."""
import argparse
import hashlib
import json
from pathlib import Path
from mixed_mcc_io import validate,frame_fields
import mixed_mcc_io
import fem_reference


def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def physical(value):
    if isinstance(value,dict):return {k:physical(v) for k,v in value.items() if k not in ('elapsed_seconds','solver_binary_sha256','maximum_absolute_material_engineering_strain_from_seed')}
    if isinstance(value,list):return [physical(v) for v in value]
    return value


def verify(study,output):
    root=Path(__file__).resolve().parents[1];old=root/'runs/mixed-mcc-ordering-resolution-20261025-v1'
    cases=[('column',study/'n4/mesh.json',old/'n4/native-finest.json',study/'n4/native-finest.json'),
           ('two-material',root/'examples/hardening-column-tet10-mesh.json',old/'two-material-example.json',study/'two-material-example.json')]
    inputs=[Path(__file__),Path(mixed_mcc_io.__file__),Path(fem_reference.__file__),study/'solver/soil-fem-3d.exe']
    reports=[]
    for name,mesh_file,before,after in cases:
        mesh=read(mesh_file);a,b=read(before),read(after)
        assert a['complete'] and b['complete']
        assert b['provenance']['solver_binary_sha256']==sha(study/'solver/soil-fem-3d.exe')
        assert json.dumps(physical(a),sort_keys=True,allow_nan=False)==json.dumps(physical(b),sort_keys=True,allow_nan=False),name+' changed histories/path or non-time diagnostics'
        validate(mesh,b);maximum=[]
        for frame in b['frames']:
            fields=frame_fields(mesh,b,frame)
            value=fields['metrics']['maximum_material_engineering_strain_from_seed'];maximum.append(value)
            assert abs(value-frame['maximum_absolute_material_engineering_strain_from_seed'])<1e-11
        reports.append(dict(name=name,physical_histories_and_primary_fields_exact=True,new_material_strain_diagnostic_recomputed=maximum))
        inputs.extend([mesh_file,before,after])
    report=dict(schema='soil-fem-stage-control-mixed-compatibility/1',passed=True,goal_complete=False,p6_complete=False,
        scope='Existing horizontal confined column and heterogeneous displacement/removal case only; intentional general-geostatic domain changes are outside equivalence',
        cases=reports,excluded_fields=['elapsed_seconds','solver_binary_sha256','maximum_absolute_material_engineering_strain_from_seed'],
        added_strain_diagnostic_checked_independently=True,input_and_source_sha256={str(p):sha(p) for p in inputs})
    with output.open('x',encoding='utf-8') as stream:json.dump(report,stream,indent=2)
    print(json.dumps(dict(passed=True,cases=reports)),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('study',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();verify(a.study.resolve(),a.output.resolve())
