"""Input-only 1D FV comparison and separate space/time refinement of 3D RT0.

The column is 0.01 x 0.01 x 1 m. This verifies a symmetry-reduced infiltration
benchmark using actual tetrahedra, not arbitrary finite-width 3D seepage.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import numpy as np
import fem_reference as fem
from richards_column_reference import solve as reference,retention
from richards_io import geometry
from verify_checkpoint import invoke
from verify_richards import audit
from water_reference import column_mesh


def write(p,value):
    with p.open('x',encoding='utf-8') as f:json.dump(value,f,indent=2)
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def mesh_column(n,width=.01):
    mesh=column_mesh(n,False);topology=fem.topology(mesh)
    mesh['face_sets']={name:[list(f) for f,o in topology.items() if len(o)==1 and all(abs(mesh['nodes'][i][2]-z)<1e-12 for i in f)] for name,z in [('bottom',0.),('top',1.)]}
    for point in mesh['nodes']:point[0]*=width;point[1]*=width
    for region in mesh['regions']:
        for p in region['bounds_m']:p[0]*=width;p[1]*=width
        region['volume_m3']*=width**2
    return mesh


def endpoints(run):return {name:next(f for f in reversed(run['frames']) if f['name']==name) for name in ('rain','pause','rain_again')}
def profile(frame,n):return np.asarray(frame['state']['relative_cell_head_m']).reshape(n,6).mean(axis=1)
def rms(a):return float(np.linalg.norm(a)/np.sqrt(len(a)))
def front(z,change):
    wet=np.flatnonzero(change>=.005)
    if not len(wet):return 0.
    i=wet[0]
    height=z[0] if i==0 else z[i-1]+(.005-change[i-1])/(change[i]-change[i-1])*(z[i]-z[i-1])
    return float(1-height)


def verify(solver,output):
    root=Path(__file__).resolve().parents[1];output.mkdir(parents=True,exist_ok=False)
    physical=read(root/'examples/richards-rain-pause.json');write(output/'reference-case.json',physical)
    refs={}
    for n,dt in [(200,12.5),(400,12.5),(800,12.5),(400,25.),(400,50.)]:
        r=reference(physical,n,dt);r['inputs_sha256']={str(output/'reference-case.json'):sha(output/'reference-case.json'),str(root/'tools/richards_column_reference.py'):sha(root/'tools/richards_column_reference.py')}
        refs[n,dt]=r;write(output/f'fv-{n}-dt{dt:g}.json',r)
    runs={};meshes={};audits={}
    for n,dt in [(32,50.),(64,50.),(128,50.),(64,25.),(64,12.5)]:
        folder=output/f'rt0-{n}-dt{dt:g}';folder.mkdir();mesh=mesh_column(n);case=copy.deepcopy(physical)
        case['settings']=dict(initial_time_step_s=dt,maximum_time_step_s=dt,absolute_flow_tolerance_m3_s=1e-22,store_integration_points=False)
        write(folder/'mesh.json',mesh);write(folder/'case.json',case)
        invoke([str(solver),'solve-richards',str(folder/'mesh.json'),str(folder/'case.json'),str(folder/'result.json')],folder/'solve.log')
        run=read(folder/'result.json');assert run['complete'] and run['solver_binary_sha256']==sha(solver)
        checked=audit(mesh,run);write(folder/'verification.json',checked);audits[n,dt]=checked;runs[n,dt]=run;meshes[n]=mesh
    rows=[];ref=refs[800,12.5];space=[];time=[];reference_refinements=[]
    for stage_index,rf in enumerate(ref['frames']):
        stage=rf['name'];zref=np.asarray(rf['z_m']);href=np.asarray(rf['head_m'])
        # Compare refinement differences on a common fine sample grid.
        sampled=[]
        for n in (32,64,128):
            f=endpoints(runs[n,50.])[stage];h=profile(f,n);z=(np.arange(n)+.5)/n
            sampled.append(np.interp(zref,z,h));rms_error=rms(h-np.interp(z,zref,href))
            row=dict(stage=stage,layers=n,time_step_s=50.,head_rms_error_m=rms_error,
                storage_offset_per_area_m=sum(f['state']['cell_storage_offset_m3'])/1e-4,
                storage_error_per_area_m=abs(sum(f['state']['cell_storage_offset_m3'])/1e-4-rf['storage_offset_m']))
            if n==128:
                mesh=meshes[n];law=physical['materials'][0]['retention'];theta=[];initial=[]
                for e in range(len(mesh['tetrahedra'])):
                    g=geometry(mesh,e);psi=f['state']['relative_cell_head_m'][e]-g['position'][:,2]
                    theta.append(float(np.average(retention(psi,law)[0],weights=g['weight'])));initial.append(float(np.average(retention(-g['position'][:,2],law)[0],weights=g['weight'])))
                change=(np.asarray(theta)-initial).reshape(n,6).mean(axis=1);ref_change=np.asarray(rf['water_content'])-retention(-zref,law)[0]
                row.update(front_depth_m=front(z,change),reference_front_depth_m=front(zref,ref_change))
                assert abs(row['front_depth_m']-row['reference_front_depth_m'])<2/n,row
                assert row['head_rms_error_m']<.002 and row['storage_error_per_area_m']<5e-5,row
            rows.append(row)
        coarse=rms(sampled[0]-sampled[1]);fine=rms(sampled[1]-sampled[2]);assert fine<.8*coarse,(stage,coarse,fine)
        space.append(dict(stage=stage,rt0_32_to_64_head_rms_m=coarse,rt0_64_to_128_head_rms_m=fine))
        hs=[profile(endpoints(runs[64,dt])[stage],64) for dt in (50.,25.,12.5)]
        coarse=rms(hs[0]-hs[1]);fine=rms(hs[1]-hs[2]);assert fine<.8*coarse,(stage,coarse,fine)
        time.append(dict(stage=stage,rt0_dt50_to_25_head_rms_m=coarse,rt0_dt25_to_12_5_head_rms_m=fine))
        spatial_ref=[np.interp(zref,refs[n,12.5]['frames'][stage_index]['z_m'],refs[n,12.5]['frames'][stage_index]['head_m']) for n in (200,400,800)]
        temporal_ref=[np.asarray(refs[400,dt]['frames'][stage_index]['head_m']) for dt in (50.,25.,12.5)]
        sr=[rms(spatial_ref[0]-spatial_ref[1]),rms(spatial_ref[1]-spatial_ref[2])];tr=[rms(temporal_ref[0]-temporal_ref[1]),rms(temporal_ref[1]-temporal_ref[2])]
        assert sr[1]<.8*sr[0] and tr[1]<.8*tr[0],(stage,sr,tr)
        reference_refinements.append(dict(stage=stage,fv_spatial_rms_m=sr,fv_temporal_rms_m=tr))
    report=dict(schema='soil-fem-richards-column-convergence/1',passed=True,solver_sha256=sha(solver),comparisons=rows,spatial_refinement=space,
        temporal_refinement=time,independent_reference_refinement=reference_refinements,
        maximum_native_integrated_mass_error_m3=max(v['integrated_mass_error_m3'] for v in audits.values()),
        maximum_reference_integrated_mass_error_m=max(v['integrated_mass_error_m'] for v in refs.values()),
        scope='3D RT0 thin-column rain/pause/restart vs independent 1D FV. Separate space and time refinement. This does not certify arbitrary finite-width 3D front/seepage lengths or mechanical suction coupling.')
    write(output/'verification.json',report);print('Richards space/time and independent FV comparisons passed')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args();verify(a.solver.resolve(),a.output_dir.resolve())
