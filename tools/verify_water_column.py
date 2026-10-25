"""Resolve FE-space error separately from cut-volume quadrature error."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from water_reference import column_mesh,continuum,q2_column
from verify_checkpoint import invoke


def write(p,d):
    with p.open('x') as f:json.dump(d,f,indent=2)


def verify(solver,output):
    output.mkdir(parents=True,exist_ok=False);root=Path(__file__).resolve().parents[1]
    ground=root/'examples/ground-water-column-mc.json';case=json.loads((root/'examples/stages-water-column.json').read_text())
    case['stages']=case['stages'][:1];case['verification_reject_once_stages']=[];case['stages'][0]['initial_increment']=1.
    sp=output/'stages.json';write(sp,case);rows=[]
    for layers in (2,4,8,16):
        for order in (3,4,5):
            label=f'n{layers}-q{order}';mesh=column_mesh(layers,True,order);mp=output/(label+'-mesh.json');rp=output/(label+'-result.json');write(mp,mesh)
            invoke([str(solver),'solve-stages',str(mp),str(ground),str(sp),str(rp)],output/(label+'.log'));d=json.loads(rp.read_text());assert d['complete']
            z=np.asarray(mesh['nodes'])[:,2];native=np.asarray(d['last_committed']['displacement_m'])[:,2];exact=continuum(z,.37)
            zz,fe=q2_column(layers,.37);expected=fe[np.rint(z*2*layers).astype(int)]
            error=float(np.max(abs(native-exact)));quadrature=float(np.max(abs(native-expected)))
            rows.append({'layers':layers,'quadrature_order':order,'maximum_nodal_continuum_error_m':error,'maximum_nodal_exact_q2_fe_error_m':quadrature,
                'exact_q2_space_error_m':float(np.max(abs(fe-continuum(zz,.37)))),'top_displacement_m':float(native[np.argmax(z)])})
            print(label,error,quadrature,flush=True)
    fine=[r for r in rows if r['quadrature_order']==5];assert fine[-1]['maximum_nodal_continuum_error_m']<fine[0]['maximum_nodal_continuum_error_m']/16
    assert fine[-1]['maximum_nodal_exact_q2_fe_error_m']<1e-8
    assert fine[-1]['maximum_nodal_continuum_error_m']<1e-7
    report={'schema':'soil-fem-water-column-refinement/1','passed':True,'solver_sha256':hashlib.sha256(solver.read_bytes()).hexdigest(),'cases':rows,
        'scope':'Prescribed hydrostatic drawdown, homogeneous MC elastic branch, tied-plane 3D Tet10. Continuum integral and exact Q2 FE load split at the water horizon distinguish spatial interpolation error from material-volume quadrature error. Fixed total unit weight; no Richards front or consolidation time.'}
    write(output/'verification.json',report);print('Water-column continuum and exact-FE refinement passed')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();verify(a.solver.resolve(),a.output_dir.resolve())
