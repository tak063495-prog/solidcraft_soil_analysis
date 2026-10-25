"""Validate a native SolidCraft analysis capture and import its body-resolved STEP.

The archived source design is never executed or read by the older portable CLI.
Only verified STEP and explicit ground/stage inputs enter the existing mesher.
"""
import argparse
import json
import math
from pathlib import Path
import re
import shutil
from verify_acceptance_packet import sha

MAX_STEP=256*1024**2
MAX_TOTAL=512*1024**2
BASE={'ground.json','stages.json','input-snapshot.json','snapshot-canonical.json','source.solidcraft.json'}


def read(path,limit=1024*1024):
    if path.stat().st_size>limit:raise ValueError('capture JSON exceeds size limit')
    def invalid(value):raise ValueError('non-finite JSON number: '+value)
    def real(value):
        number=float(value)
        if not math.isfinite(number):invalid(value)
        return number
    return json.loads(path.read_text(encoding='utf-8'),parse_constant=invalid,parse_float=real)


def equal(a,b):
    # Python's True == 1 must not accept a different Rust JSON value.
    if type(a) is not type(b):return False
    if isinstance(a,dict):return a.keys()==b.keys() and all(equal(a[k],b[k]) for k in a)
    if isinstance(a,list):return len(a)==len(b) and all(equal(x,y) for x,y in zip(a,b))
    return a==b


def validate(folder,*,manifest_name='capture.json'):
    if manifest_name not in ('capture.json','native-capture.json'):raise ValueError('unsupported capture manifest name')
    folder=folder.resolve();capture=read(folder/manifest_name)
    if capture.get('schema')!='solidcraft-soil-capture/1' or capture.get('complete') is not True:raise ValueError('complete native analysis capture required')
    digest=capture.get('input_sha256');exporter=capture.get('exporter_binary_sha256')
    if not all(isinstance(h,str) and re.fullmatch('[0-9a-f]{64}',h) for h in (digest,exporter)):raise ValueError('missing input/exporter fingerprint')
    records=capture.get('bodies',[])
    if not isinstance(records,list) or not 1<=len(records)<=256:raise ValueError('expected 1..256 captured bodies')
    files=capture.get('files_sha256',{})
    if not isinstance(files,dict) or not BASE.issubset(files) or len(files)!=len(records)+len(BASE):raise ValueError('capture file inventory is incomplete')
    total=0
    for name,h in files.items():
        if name not in BASE and not re.fullmatch(r'body-[0-9]{4}\.step',name):raise ValueError('invalid capture filename')
        path=folder/name
        if not path.resolve().is_relative_to(folder) or not path.is_file():raise ValueError('capture path escapes its directory')
        size=path.stat().st_size;total+=size
        limit=MAX_STEP if name.endswith('.step') else (64*1024**2 if name=='source.solidcraft.json' else 1024**2)
        if size>limit or total>MAX_TOTAL:raise ValueError('capture exceeds bounded size')
        if sha(path)!=h:raise ValueError('capture file hash mismatch: '+name)
    if sha(folder/'snapshot-canonical.json')!=digest:raise ValueError('canonical snapshot fingerprint mismatch')
    wrapped=read(folder/'input-snapshot.json');payload=read(folder/'snapshot-canonical.json')
    if wrapped.get('input_sha256')!=digest or not equal(wrapped.get('snapshot'),payload):raise ValueError('snapshot wrapper differs from canonical payload')
    if payload.get('schema')!='solidcraft-soil-snapshot/1' or payload.get('cad_length_unit')!='mm' or payload.get('analysis_units')!='m-kN-kPa':raise ValueError('capture units/schema mismatch')
    ground=read(folder/'ground.json');stages=read(folder/'stages.json')
    if not equal(ground,payload.get('ground')) or not equal(stages,payload.get('stages')):raise ValueError('ground/stages differ from the captured snapshot')
    if ground.get('schema')!='soil-fem-ground/1' or ground.get('units')!='m-kN-kPa' or stages.get('schema') not in ('soil-fem-stages/1','soil-fem-projected-mcc-stages/1'):raise ValueError('unsupported FEM input schema')
    field='bodies' if ground.get('initialization',{}).get('method')=='equilibrate' else 'layers'
    entries=ground.get(field,[]);mapping={x['body']:x['material_id'] for x in entries};tools=ground.get('excavation_bodies',[])
    if not mapping or len(mapping)!=len(entries) or len(set(tools))!=len(tools) or set(mapping)&set(tools):raise ValueError('duplicate/missing material or excavation labels')
    expected=set(mapping)|set(tools);bindings=payload.get('bindings',[])
    aliases=[r.get('analysis_body') for r in records];cad=[r.get('cad_body') for r in records];steps=[r.get('step') for r in records]
    if len(set(aliases))!=len(records) or set(aliases)!=expected or len(set(cad))!=len(records) or len(set(steps))!=len(records):raise ValueError('body mapping is missing, duplicate or ambiguous')
    if len(bindings)!=len(records):raise ValueError('snapshot body bindings are incomplete')
    bound={b['analysis_body']:b['cad_body'] for b in bindings}
    if len(bound)!=len(bindings) or bound!={r['analysis_body']:r['cad_body'] for r in records}:raise ValueError('body mapping differs from snapshot')
    if set(files)!=BASE|set(steps):raise ValueError('STEP inventory differs from body mapping')
    for r in records:
        if r.get('cad_length_unit')!='mm' or files.get(r['step'])!=r.get('source_sha256'):raise ValueError('body STEP units/hash mismatch')
    return capture,ground,stages,mapping


def import_capture(capture_folder,output):
    capture_folder=capture_folder.resolve();output=output.resolve()
    manifest_hash=sha(capture_folder/'capture.json')
    capture,ground,stages,mapping=validate(capture_folder)
    output.mkdir(parents=True,exist_ok=False)
    for name,h in capture['files_sha256'].items():
        shutil.copyfile(capture_folder/name,output/name)
        if sha(output/name)!=h:raise ValueError('input changed while copying capture')
    shutil.copyfile(capture_folder/'capture.json',output/'native-capture.json')
    if sha(output/'native-capture.json')!=manifest_hash or not equal(read(output/'native-capture.json'),capture):raise ValueError('capture manifest changed during import')
    # Gmsh is an out-of-process mesher dependency of the FEM project. SolidCraft
    # neither links this kernel nor depends on Python as an application runtime.
    import gmsh
    gmsh.initialize()
    try:
        gmsh.option.setNumber('General.Terminal',0)
        gmsh.option.setString('Geometry.OCCTargetUnit','MM')
        regions=[];tools=[]
        for i,r in enumerate(sorted(capture['bodies'],key=lambda v:v['analysis_body']),1):
            gmsh.clear();entities=gmsh.model.occ.importShapes(str(output/r['step']))
            volumes=[x for x in entities if x[0]==3]
            if len(volumes)!=1:raise ValueError('each captured body must contain exactly one closed STEP solid')
            volume=gmsh.model.occ.getMass(*volumes[0])*1e-9
            if not 0<volume<float('inf'):raise ValueError('invalid STEP volume')
            bounds=gmsh.model.occ.getBoundingBox(*volumes[0]);bounds=[[v*1e-3 for v in bounds[:3]],[v*1e-3 for v in bounds[3:]]]
            row=dict(id=i,body=r['analysis_body'],step=r['step'],source_sha256=r['source_sha256'],bounds_m=bounds,volume_m3=volume)
            if r['analysis_body'] in mapping:row['material_id']=mapping[r['analysis_body']];regions.append(row)
            else:tools.append(row)
    finally:gmsh.finalize()
    manifest=dict(schema='soil-fem-cad-bodies/1',units='m-kN-kPa',ground_sha256=sha(output/'ground.json'),regions=regions,
        volume_measurement=dict(analysis='STEP B-rep',source='native SolidCraft capture; independent Gmsh volume'))
    if tools:manifest['excavation_tools']=tools
    link=dict(schema='soil-fem-native-capture-import/1',passed=True,input_sha256=capture['input_sha256'],
        exporter_binary_sha256=capture['exporter_binary_sha256'],source_capture_sha256=sha(output/'native-capture.json'),
        source_files_sha256=capture['files_sha256'],cad_body_to_analysis_body={r['cad_body']:r['analysis_body'] for r in capture['bodies']},
        source_design_execution='never executed; archived bytes only',equilibrium_verified=False)
    for name,value in (('bodies.json',manifest),('native-capture-import.json',link)):
        with (output/name).open('x',encoding='utf-8') as f:json.dump(value,f,indent=2,allow_nan=False)
    print('Native SolidCraft capture imported:',len(regions),'material bodies,',len(tools),'excavation tools')
    return manifest,link


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('capture',type=Path);p.add_argument('output',type=Path);a=p.parse_args();import_capture(a.capture,a.output)
