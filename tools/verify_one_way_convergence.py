"""Independent FV water + continuous elastic column, and separate refinements."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import argparse
import copy
import json
from pathlib import Path
import numpy as np
from richards_column_reference import solve as flow_reference,retention
from verify_richards_convergence import mesh_column
from one_way_reference import audit
from verify_one_way_cases import read,write
from verify_checkpoint import invoke
from verify_acceptance_packet import sha


def continuum(z,head,law):
    """Input-only u(z) from stress equilibrium with chi=Se and no exterior water.

    Integrate each FV cell separately and split at the target height. The
    independent 1D water heads, not native RT0 states, generate this response.
    """
    n=len(head);edges=np.unique(np.r_[np.linspace(0,1,n+1),z]);x,w=np.polynomial.legendre.leggauss(8)
    pos=(edges[:-1,None]+edges[1:,None])/2+(edges[1:]-edges[:-1])[:,None]*x/2
    weight=np.broadcast_to((edges[1:]-edges[:-1])[:,None]*w/2,pos.shape);idx=np.minimum((pos*n).astype(int),n-1)
    psi=np.asarray(head)[idx]-pos;theta,_,_=retention(psi,law);se=(theta-law['residual_water_content'])/(law['saturated_water_content']-law['residual_water_content'])
    water=theta+law['saturated_specific_storage_per_m']*np.maximum(psi,0);effective=se*10*psi
    return float(-(10*np.sum(weight*np.minimum(pos,z)*water)-np.sum(weight*(pos<z)*effective))/12000.)


def rms(v):return float(np.sqrt(np.mean(np.asarray(v)**2)))
def endpoint(run,name):return next(f for f in reversed(run['frames']) if f['phase']=='transient' and f['name']==name)
def profile(mesh,run,frame):
    z=np.asarray(mesh['nodes'])[:,2];u=np.asarray(frame['state']['mechanical']['displacement_m'])[:,2];ref=np.asarray(run['initial_reference_displacement_m'])[:,2];height=np.unique(z)
    return height,np.array([u[z==s].mean() for s in height]),np.array([(u-ref)[z==s].mean() for s in height])


def verify(solver,output):
    root=Path(__file__).resolve().parents[1];output.mkdir(parents=True,exist_ok=False)
    ground=read(root/'examples/ground-one-way-mc.json');case=read(root/'examples/one-way-mc.json');law=case['flow']['materials'][0]['retention'];refs={}
    for n,dt in [(200,12.5),(400,12.5),(800,12.5),(400,25.),(400,50.)]:
        value=flow_reference(case['flow'],n,dt);value['inputs_sha256']={str(root/'examples/one-way-mc.json'):sha(root/'examples/one-way-mc.json'),str(root/'tools/richards_column_reference.py'):sha(root/'tools/richards_column_reference.py')}
        refs[n,dt]=value;write(output/f'fv-{n}-dt{dt:g}.json',value)
    runs={};meshes={};checks={}
    for n,dt in [(16,50.),(32,50.),(64,50.),(32,25.),(32,12.5)]:
        folder=output/f'rt0-{n}-dt{dt:g}';folder.mkdir();mesh=mesh_column(n);value=copy.deepcopy(case)
        value['store_transfer_integration_points']=False
        value['flow']['settings'].update(initial_time_step_s=dt,maximum_time_step_s=dt,absolute_flow_tolerance_m3_s=1e-22)
        value['mechanics']['settings'].update(relative_force_tolerance=1e-11,relative_linear_tolerance=1e-11)
        for label,v in [('mesh',mesh),('ground',ground),('case',value)]:write(folder/f'{label}.json',v)
        invoke([str(solver),'solve-one-way',str(folder/'mesh.json'),str(folder/'ground.json'),str(folder/'case.json'),str(folder/'result.json')],folder/'solve.log')
        run=read(folder/'result.json');assert run['complete'] and run['solver_binary_sha256']==sha(solver)
        checked=audit(mesh,run);write(folder/'verification.json',checked);runs[n,dt]=run;meshes[n]=mesh;checks[n,dt]=checked
        print('ONE_WAY column finished',n,dt,flush=True)
    comparisons=[];space=[];time=[];refinements=[]
    for j,period in enumerate(case['flow']['analysis']['periods']):
        stage=period['name'];zfine=np.linspace(0,1,65);reference=refs[800,12.5]['frames'][j]
        exact=np.array([continuum(z,reference['head_m'],law) for z in zfine]);initial=np.array([continuum(z,np.zeros(800),law) for z in zfine]);increment=exact-initial
        sampled=[]
        for n in (16,32,64):
            run=runs[n,50.];f=endpoint(run,stage);z,u,du=profile(meshes[n],run,f)
            expected=np.interp(z,zfine,exact);expected_increment=np.interp(z,zfine,increment)
            row=dict(stage=stage,layers=n,time_step_s=50.,absolute_displacement_rms_error_m=rms(u-expected),incremental_displacement_rms_error_m=rms(du-expected_increment),maximum_increment_difference_m=float(np.max(np.abs(du-expected_increment))),maximum_reference_displacement_increment_m=float(np.max(np.abs(expected_increment))))
            if n==64:assert row['incremental_displacement_rms_error_m']<1e-6 and row['maximum_increment_difference_m']<2e-6,row
            comparisons.append(row);sampled.append(np.interp(zfine,z,du));print(row,flush=True)
        gaps=[rms(sampled[0]-sampled[1]),rms(sampled[1]-sampled[2])];assert gaps[1]<.85*gaps[0],(stage,gaps)
        space.append(dict(stage=stage,increment_rms_gaps_m=gaps))
        values=[]
        for dt in (50.,25.,12.5):
            z,_,du=profile(meshes[32],runs[32,dt],endpoint(runs[32,dt],stage));values.append(du)
        gaps=[rms(values[0]-values[1]),rms(values[1]-values[2])];assert gaps[1]<.85*gaps[0],(stage,gaps)
        time.append(dict(stage=stage,increment_rms_gaps_m=gaps))
        # Reference water space and time are refined independently too.
        spatial=[np.array([continuum(z,refs[n,12.5]['frames'][j]['head_m'],law) for z in zfine]) for n in (200,400,800)]
        temporal=[np.array([continuum(z,refs[400,dt]['frames'][j]['head_m'],law) for z in zfine]) for dt in (50.,25.,12.5)]
        sg=[rms(spatial[0]-spatial[1]),rms(spatial[1]-spatial[2])];tg=[rms(temporal[0]-temporal[1]),rms(temporal[1]-temporal[2])]
        assert sg[1]<.85*sg[0] and tg[1]<.85*tg[0],(stage,sg,tg)
        refinements.append(dict(stage=stage,independent_spatial_displacement_rms_gaps_m=sg,independent_temporal_displacement_rms_gaps_m=tg))
    report=dict(schema='soil-fem-one-way-column-convergence/1',passed=True,solver_sha256=sha(solver),comparisons=comparisons,native_spatial_refinement=space,native_temporal_refinement=time,independent_reference_refinement=refinements,
        maximum_hydraulic_integrated_mass_error_m3=max(v['hydraulic_audit']['integrated_mass_error_m3'] for v in checks.values()),
        scope='Actual 3D thin tetrahedral column with tied horizontal planes vs independent 1D FV water and continuous elastic stress-equilibrium integrals. Independent water/time and native space/time refinements. Does not certify arbitrary 3D plastic deformation, finite-width rainfall or slope failure.')
    write(output/'verification.json',report)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args();verify(a.solver.resolve(),a.output_dir.resolve())
