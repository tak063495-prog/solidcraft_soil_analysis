"""Independent 3D steady FV / RT0 comparison, including finite-width inflow.

The benchmarks are unit homogeneous cubes. Full-width reservoir is a 2D
symmetry check; the central half-width reservoir produces genuine q_y.
No slope mechanics, excavation-domain change or arbitrary geology is implied.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from richards_dam_mesh import make
from richards_dam_reference import solve as reference
from richards_dam_audit import audit
from verify_checkpoint import invoke


def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,v):
    with p.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2)
def rms(v):return float(np.sqrt(np.mean(np.asarray(v)**2)))


def resample(a,n):
    """Trilinear cell-centre sampling, clamped only outside centre coverage."""
    a=np.asarray(a);m=a.shape[0];assert a.shape==(m,m,m)
    s=np.clip(m*(np.arange(n)+.5)/n-.5,0,m-1)
    low=np.minimum(np.floor(s).astype(int),m-2);fraction=s-low
    out=np.zeros((n,n,n))
    for dz in (0,1):
        for dy in (0,1):
            for dx in (0,1):
                wz=fraction if dz else 1-fraction;wy=fraction if dy else 1-fraction;wx=fraction if dx else 1-fraction
                out+=a[(low+dz)[:,None,None],(low+dy)[None,:,None],(low+dx)[None,None,:]]*wz[:,None,None]*wy[None,:,None]*wx[None,None,:]
    return out


def waterline(head):
    n=head.shape[0];z=(np.arange(n)+.5)/n;out=np.zeros((n,n))
    for j in range(n):
        for i in range(n):
            psi=head[:,j,i]-z;dry=np.flatnonzero(psi<=0)
            if not len(dry):out[j,i]=1.;continue
            k=int(dry[0]);out[j,i]=z[0] if k==0 else z[k-1]-psi[k-1]/(n*(psi[k]-psi[k-1]))
    return out


def verify(solver,output):
    root=Path(__file__).resolve().parents[1];output.mkdir(parents=True,exist_ok=False)
    case=read(root/'examples/richards-dam.json');write(output/'case.json',case)
    # Independent solver has an intentionally bounded benchmark contract.
    assert case['initial_head']=={'head_m':.5,'gradient':[-.25,0.,0.]}
    assert case['analysis']=={'kind':'steady','boundaries':[{'face_set':'reservoir','condition':{'kind':'head','value':{'head_m':.5}}},{'face_set':'air','condition':{'kind':'atmospheric_no_ponding','rainfall_vertical_m_s':0.}}]}
    assert case['materials']==[dict(id=i,saturated_conductivity_m_s=(np.eye(3)*1e-6).tolist(),retention=dict(residual_water_content=.08,saturated_water_content=.45,alpha_per_m=1.2,n=1.6,saturated_specific_storage_per_m=1e-5)) for i in (1,2)]
    all_rows=[];refinements=[];native_refinements=[]
    for finite_width,label in [(False,'full'),(True,'strip')]:
        refs={}
        for n in (16,32,64):
            ref=reference(n,finite_width);ref['input_sha256']={str(output/'case.json'):sha(output/'case.json'),str(root/'tools/richards_dam_reference.py'):sha(root/'tools/richards_dam_reference.py')}
            write(output/f'fv-{label}-{n}.json',ref);refs[n]=ref
        samples=[resample(refs[n]['head_m'],32) for n in (16,32,64)]
        head_diffs=[rms(samples[0]-samples[1]),rms(samples[1]-samples[2])]
        flow_diffs=[abs(refs[16]['right_outflow_m3_s']-refs[32]['right_outflow_m3_s']),abs(refs[32]['right_outflow_m3_s']-refs[64]['right_outflow_m3_s'])]
        assert head_diffs[1]<.85*head_diffs[0] and flow_diffs[1]<.85*flow_diffs[0],(label,head_diffs,flow_diffs)
        if finite_width:assert refs[64]['transverse_flux_rms_m_s']>1e-8
        else:assert refs[64]['transverse_flux_rms_m_s']<1e-14
        refinements.append(dict(case=label,fv16_to_32_head_rms_m=head_diffs[0],fv32_to_64_head_rms_m=head_diffs[1],fv16_to_32_right_flow_difference_m3_s=flow_diffs[0],fv32_to_64_right_flow_difference_m3_s=flow_diffs[1],finest_reference=refs[64]['right_outflow_m3_s'],finest_reference_seepage_area_m2=refs[64]['right_seepage_area_m2']))
        coarse_samples=[]
        for n in (4,8,16,24):
            folder=output/f'rt0-{label}-{n}';folder.mkdir();mesh=make(n,finite_width);write(folder/'mesh.json',mesh)
            invoke([str(solver),'solve-richards',str(folder/'mesh.json'),str(output/'case.json'),str(folder/'result.json')],folder/'solve.log')
            result=read(folder/'result.json');assert result['complete'] and result['solver_binary_sha256']==sha(solver)
            checked=audit(mesh,result);write(folder/'verification.json',checked)
            head=np.asarray(result['last_committed']['relative_cell_head_m']).reshape(n,n,n,6).mean(axis=3)
            reference_head=resample(refs[64]['head_m'],n);h_error=rms(head-reference_head);phreatic_error=rms(waterline(head)-waterline(reference_head))
            # Recover the *volume* mean q_y from the four Piola basis values at
            # each straight-tet centroid. No nodal flux smoothing is involved.
            points=np.asarray(mesh['nodes'])[np.asarray(mesh['tetrahedra'])];j=np.swapaxes(points[:,1:]-points[:,0,None],1,2);det=np.linalg.det(j)
            reference_vertices=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
            basis=np.einsum('eij,fj->efi',j,2*(.25-reference_vertices))/det[:,None,None]
            q=np.einsum('efi,ef->ei',basis,np.asarray(result['last_committed']['cell_outflow_m3_s']));q_y_rms=rms(q[:,1])
            flow_error=abs(checked['right_outflow_m3_s']-refs[64]['right_outflow_m3_s'])/refs[64]['right_outflow_m3_s']
            area_error=abs(checked['seepage_area_m2']-refs[64]['right_seepage_area_m2'])
            row=dict(case=label,subdivisions_per_axis=n,tetrahedra=len(mesh['tetrahedra']),head_rms_error_m=h_error,phreatic_height_rms_difference_m=phreatic_error,right_flow_relative_difference=flow_error,right_outflow_m3_s=checked['right_outflow_m3_s'],right_seepage_area_m2=checked['seepage_area_m2'],seepage_area_difference_m2=area_error,transverse_cell_mean_flux_rms_m_s=q_y_rms,maximum_native_local_scaled_residual=result['frames'][-1]['balance']['maximum_scaled_residual'],independent_global_mass_error_m3_s=checked['independent_boundary_mass_error_m3_s'])
            all_rows.append(row);coarse_samples.append(resample(head,24));print(label,n,row,flush=True)
            if n==24:
                assert h_error<.0015 and flow_error<.015,(label,row)
                # Endpoint support is quantized by face size; report its full
                # mesh uncertainty, not fictitious decimal-place convergence.
                assert area_error<2/n+2/64 and phreatic_error<2/n,(label,row)
                if finite_width:assert q_y_rms>1e-8
                else:assert q_y_rms<2e-8  # finite tetrahedral discretization error, not exact symmetry
        selected=[r for r in all_rows if r['case']==label]
        assert all(selected[i+1]['head_rms_error_m']<.85*selected[i]['head_rms_error_m'] for i in range(3)),selected
        assert all(selected[i+1]['right_flow_relative_difference']<.85*selected[i]['right_flow_relative_difference'] for i in range(3)),selected
        native_refinements.append(dict(case=label,head_rms_differences_m=[rms(coarse_samples[i]-coarse_samples[i+1]) for i in range(3)],seepage_area_error_bound_m2=2/24+2/64,scope='All three axes refined together; face support quantization is retained, not asserted strictly monotone.'))
    report=dict(schema='soil-fem-richards-dam-verification/1',passed=True,solver_sha256=sha(solver),comparisons=all_rows,independent_reference_refinements=refinements,native_refinements=native_refinements,
        scope='Steady homogeneous isotropic unit-cube full-width and central-half-width reservoir benchmarks. Finite-width case has genuine transverse flow. Independent 3D FV and all-cell local RT0 audit; no mechanical coupling, excavation domain transfer or general geology certification.')
    write(output/'verification.json',report)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args();verify(a.solver.resolve(),a.output_dir.resolve())
