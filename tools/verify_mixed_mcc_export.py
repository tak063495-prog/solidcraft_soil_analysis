"""Roundtrip actual mixed construction output and reject corrupted primary data."""
import argparse
import copy
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import zipfile
import meshio
import numpy as np
from export_stages import export


def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,d):
    with p.open('x',encoding='utf-8') as f:json.dump(d,f,indent=2)


def verify(output,base=None,example=None):
    root=Path(__file__).resolve().parents[1]
    base=base or root/'runs/mixed-mcc-independent-20261025-v2'
    example=example or root/'runs/mixed-mcc-export-20261025-v1/two-material-example.json'
    cases=[('column',base/'n4/mesh.json',base/'n4/native-finest.json'),
        ('two-material',root/'examples/hardening-column-tet10-mesh.json',example)]
    legacy=root/'runs/mixed-mcc-legacy-20261025-v1'
    archive=root/'runs/mixed-mcc-development-20261025-v1/verification-sources/source-snapshot.zip'
    legacy_cases=[('legacy-mcc',root/'runs/mcc-tet10-cycle-acceptance-20261025-v1/inputs/mesh.json',legacy/'legacy-mcc.json'),
        ('legacy-mc',root/'runs/2d-3d-column-acceptance-20261025-v1/data/mesh-3d.json',legacy/'legacy-mc.json')]
    inputs=[Path(__file__),root/'tools/export_stages.py',root/'tools/high_order_io.py',root/'tools/mixed_mcc_io.py',root/'tools/fem_reference.py',archive]
    inputs.extend(p for _,m,r in cases+legacy_cases for p in (m,r))
    fingerprints={str(p):sha(p) for p in inputs}
    output.mkdir(exist_ok=False);reports=[]
    for tag,mesh_path,result_path in cases:
        mesh=read(mesh_path);result=read(result_path);target=output/tag
        export(mesh_path,result_path,target)
        counts=[];volume_difference=0.;jumps=[]
        for i,frame in enumerate(result['frames']):
            primary=meshio.read(target/f'mean-pressure-dofs-{i:03}.vtu')
            p=frame['state']['mixed_mcc']['pressure_dofs'];counts.append(len(p))
            assert np.array_equal(primary.point_data['mean_effective_stress_compression_positive_kpa'],[d['compression_positive_kpa'] for d in p])
            by_node={}
            for row in p:by_node.setdefault(row['node_index'],[]).append(row)
            for rows in by_node.values():
                if len(rows)>1:jumps.append(max(d['compression_positive_kpa'] for d in rows)-min(d['compression_positive_kpa'] for d in rows))
            raw=meshio.read(target/f'integration-points-{i:03}.vtu');ids=raw.point_data['integration_point_id']
            mat=np.asarray(frame['state']['mixed_mcc']['cumulative_volume_strain_compression_positive'])[ids]
            assert np.array_equal(raw.point_data['cumulative_material_volume_strain_compression_positive'],mat)
            tensor=raw.point_data['cumulative_material_engineering_strain_compression_positive']
            assert np.allclose(tensor[:,:3].sum(axis=1),mat,rtol=0.,atol=1e-15)
            volume_difference=max(volume_difference,float(abs(raw.point_data['material_minus_kinematic_volume_strain']).max()))
            with (target/f'material-history-{i:03}.csv').open(encoding='utf-8-sig',newline='') as f:history=list(csv.DictReader(f))
            assert len(history)==len(result['integration_points'])
            inactive=[h for h in history if h['active']=='0']
            assert all(h['primary_mean_effective_stress_kpa']=='' for h in inactive)
        assert counts[-1]<counts[0] and volume_difference>1e-12
        if tag=='two-material':assert max(jumps)>1.,'benchmark must contain a real material pressure jump'
        reports.append(dict(name=tag,frames=len(result['frames']),pressure_dof_counts=counts,
            tangent_method=read(target/'export-verification.json')['tangent_method'],
            maximum_preserved_interface_pressure_jump_kpa=max(jumps),
            maximum_material_kinematic_volume_difference=volume_difference,
            unsmoothed_primary_pressure_and_all_gp_history_roundtrip=True))
    # Run the exact previous exporter from its frozen source, under the same
    # runtime, to compare every legacy output byte rather than only a screenshot.
    old=output/'previous-exporter';old.mkdir()
    with zipfile.ZipFile(archive) as z:
        names={n.replace('\\','/'):n for n in z.namelist()}
        for name in ('export_stages.py','high_order_io.py'):
            (old/name).write_bytes(z.read(names['tools/'+name]))
    for tag,mesh_path,result_path in legacy_cases:
        target=output/tag;previous=output/(tag+'-previous')
        export(mesh_path,result_path,target)
        subprocess.run([sys.executable,'-X','utf8',str(old/'export_stages.py'),str(mesh_path),str(result_path),str(previous)],check=True,cwd=root)
        a={str(p.relative_to(target)):sha(p) for p in target.rglob('*') if p.is_file()}
        b={str(p.relative_to(previous)):sha(p) for p in previous.rglob('*') if p.is_file()}
        assert a==b,tag+' changed existing exported fields/files'
        reports.append(dict(name=tag,all_exported_files_byte_exact_to_previous_exporter=True,files=len(a)))
    mesh_path=cases[0][1];original=read(cases[0][2]);refused=[]
    def refuse(name,change,source=None,mesh=mesh_path):
        candidate=copy.deepcopy(source or original);change(candidate)
        path=output/(name+'.json');target=output/(name+'-must-not-exist');write(path,candidate)
        try:export(mesh,path,target)
        except ValueError as error:
            assert not target.exists(),name+' created output before validation'
            refused.append(dict(name=name,reason=str(error),output_not_created=True))
        else:raise AssertionError(name+' was accepted')
    refuse('missing-pressure',lambda d:d['frames'][1]['state']['mixed_mcc']['pressure_dofs'].pop())
    refuse('duplicate-pressure',lambda d:d['frames'][1]['state']['mixed_mcc']['pressure_dofs'].append(copy.deepcopy(d['frames'][1]['state']['mixed_mcc']['pressure_dofs'][0])))
    def pressure(d):d['frames'][1]['state']['mixed_mcc']['pressure_dofs'][0]['compression_positive_kpa']+=1.
    refuse('pressure-state-mismatch',pressure)
    def volume(d):d['frames'][1]['state']['mixed_mcc']['cumulative_volume_strain_compression_positive'][0]+=.001
    refuse('volume-void-mismatch',volume)
    def weak(d):
        s=d['frames'][1]['state'];s['mixed_mcc']['cumulative_volume_strain_compression_positive'][0]+=.001
        s['material_states'][0]['void_ratio']=(1+s['mixed_mcc']['reference_void_ratio'][0])*np.exp(-s['mixed_mcc']['cumulative_volume_strain_compression_positive'][0])-1
    refuse('weak-volume-mismatch',weak)
    def residual(d):d['frames'][1]['state']['mixed_mcc']['relative_pressure_consistency_residual']=.1
    refuse('false-residual',residual)
    def geometry(d):d['integration_points'][0]['position_m'][0]+=.01
    refuse('changed-geometry',geometry)
    def removed(d):
        s=d['frames'][-1]['state'];i=next(i for i,p in enumerate(d['integration_points']) if not s['active_elements'][p['element_index']])
        s['material_states'][i]['plastic_strain_compression_positive'][0]+=.001
    refuse('changed-removed-history',removed)
    refuse('wrong-schema',lambda d:d.update(schema='soil-fem-stage-result/2'))
    refuse('incomplete-result',lambda d:d.update(complete=False))
    def method(d):d['frames'][1]['state']['mixed_mcc']['tangent_method']='unrecognized-method'
    refuse('unknown-tangent-method',method)
    def switch(d):
        m=d['frames'][0]['state']['mixed_mcc'].get('tangent_method','finite_difference_fixed_partition')
        d['frames'][1]['state']['mixed_mcc']['tangent_method']='implicit_branch_fixed_partition' if m=='finite_difference_fixed_partition' else 'finite_difference_fixed_partition'
    refuse('changed-tangent-method',switch)
    def final(d):d['last_committed']['mixed_mcc']['reference_void_ratio'][0]+=.001
    refuse('different-final-state',final)
    def strain_cache(d):d['frames'][1]['maximum_absolute_material_engineering_strain_from_seed']=.5
    refuse('false-maximum-material-strain',strain_cache)
    def swap(d):
        p=d['frames'][0]['state']['mixed_mcc']['pressure_dofs']
        a,b=next((a,b) for a in range(len(p)) for b in range(a+1,len(p)) if p[a]['node_index']==p[b]['node_index'])
        p[a]['compression_positive_kpa'],p[b]['compression_positive_kpa']=p[b]['compression_positive_kpa'],p[a]['compression_positive_kpa']
    refuse('interface-pressure-swap',swap,read(example),cases[1][1])
    assert fingerprints=={str(p):sha(p) for p in inputs}
    report=dict(schema='soil-fem-mixed-mcc-export-verification/1',passed=True,goal_complete=False,p6_complete=False,
        scope='Output invariants, roundtrips and legacy export compatibility only; no independent physics acceptance',
        cases=reports,rejected_inputs=refused,input_and_source_sha256=fingerprints,
        evidence_sha256={str(p.relative_to(output)):sha(p) for p in output.rglob('*') if p.is_file()})
    write(output/'verification.json',report)
    print(json.dumps(dict(passed=True,cases=reports,rejected_inputs=len(refused))),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path)
    p.add_argument('--mixed-base',type=Path);p.add_argument('--example',type=Path)
    a=p.parse_args();verify(a.output.resolve(),a.mixed_base.resolve() if a.mixed_base else None,a.example.resolve() if a.example else None)
