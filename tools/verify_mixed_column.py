"""Two-material body-loaded column, pressure jump and rotated independent solve."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
from verify_incompressibility import mesh as cube
from verify_mixed_mechanics import verify
import fem_reference as ref


def run(solver,out):
    root=Path(__file__).resolve().parents[1];solver=solver.resolve();out.mkdir(exist_ok=False)
    digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    fingerprints={str(p):digest(p) for p in (solver,Path(__file__),root/'tools/verify_mixed_mechanics.py',root/'tools/verify_incompressibility.py',root/'tools/fem_reference.py')}
    reports=[]
    for rotated in (False,True):
        name='rotated' if rotated else 'column';m=cube(2,True);xyz=np.asarray(m['nodes']);m['region_ids']=[];m['material_ids']=[]
        for t in m['tetrahedra']:
            mid=1 if xyz[t,2].mean()<.5 else 2;m['material_ids'].append(mid);m['region_ids'].append(mid)
        m['regions']=[dict(id=i,body=f'Layer{i}',material_id=i,bounds_m=[[0.,0.,(i-1)*.5],[1.,1.,i*.5]],volume_m3=.5,
            source_sha256=hashlib.sha256(f'analytic layer {i}'.encode()).hexdigest()) for i in (1,2)]
        m['node_sets']['all']=list(range(len(xyz)));m['node_sets']['bottom']=np.flatnonzero(xyz[:,2]==0.).tolist()
        m['face_sets']={'top':[list(f) for f,owners in ref.topology(m).items() if len(owners)==1 and np.all(xyz[list(f),2]==1.)]}
        angle=.37;tilt=.29
        q=np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1.]])@np.array([[1.,0,0],[0,np.cos(tilt),-np.sin(tilt)],[0,np.sin(tilt),np.cos(tilt)]]) if rotated else np.eye(3)
        m['nodes']=(xyz@q.T).tolist()
        for region in m['regions']:
            nodes=np.unique(np.asarray(m['tetrahedra10'])[np.asarray(m['material_ids'])==region['id']])
            region['bounds_m']=[np.min(np.asarray(m['nodes'])[nodes],axis=0).tolist(),np.max(np.asarray(m['nodes'])[nodes],axis=0).tolist()]
        case=dict(schema='soil-fem-mechanics/1',units='m-kN-kPa',materials=[dict(id=i,young_kpa=1000.,poisson=nu,body_force_kn_m3=(q@np.array([0.,0.,-float(i)])).tolist()) for i,nu in ((1,.2),(2,.4999))],
            boundary=dict(constraints=[dict(node_set=nodes,direction=q[:,a].tolist(),value_m=0.) for nodes,a in (('all',0),('all',1),('bottom',2))],
                surface_loads=[dict(kind='pressure',face_set='top',pressure_kpa=3.)]),
            relative_linear_tolerance=1e-11,maximum_free_dofs=2000,linear_solver=dict(method='sparse_gmres_ilu0',options=dict(restart=100,maximum_iterations=5000,pivot_tolerance=1e-14)))
        mp=out/(name+'-mesh.json');cp=out/(name+'-case.json');rp=out/(name+'-result.json');vp=out/(name+'-independent.json')
        for p,value in ((mp,m),(cp,case)):
            with p.open('x') as f:json.dump(value,f)
        subprocess.run([sys.executable,str(root/'tools/run_logged.py'),str(out/(name+'.log')),str(solver),'solve-mixed-mechanics',str(mp),str(cp),str(rp)],cwd=root,check=True)
        native=json.loads(rp.read_text());actual=np.asarray(native['displacement_m'])@q
        compliance=lambda nu:(1+nu)*(1-2*nu)/(1000*(1-nu))
        z=xyz[:,2];lo=np.minimum(z,.5);hi=np.maximum(z,.5)
        exact=-(4.5*lo-.5*lo**2)*compliance(.2)-(5*hi-hi**2-2.25)*compliance(.4999)
        np.testing.assert_allclose(actual[:,2],exact,rtol=1e-7,atol=1e-11)
        np.testing.assert_allclose(actual[:,:2],0.,atol=1e-11)
        pressure_errors=[]
        for p in native['pressure_dofs']:
            nu=.2 if p['material_id']==1 else .4999;z=xyz[p['node_index'],2]
            vertical=4.5-z if p['material_id']==1 else 5-2*z
            pressure_errors.append(abs(p['compression_positive_kpa']-(1+nu)/(3*(1-nu))*vertical))
        assert max(pressure_errors)<1e-7,max(pressure_errors)
        verify(mp,cp,rp,solver,vp,5)
        independent=json.loads(vp.read_text())
        reports.append(dict(name=name,maximum_displacement_error_m=float(np.max(np.abs(actual[:,2]-exact))),
            maximum_pressure_error_kpa=float(max(pressure_errors)),independent=independent))
    assert fingerprints=={p:digest(Path(p)) for p in fingerprints},'source/solver changed during column verification'
    record=dict(passed=True,goal_complete=False,p6_complete=False,scope='Finite-bulk linear two-material column and rotated copy; independent fields and analytic pressure jump',
                cases=reports,input_and_source_sha256=fingerprints,evidence_sha256={p.name:digest(p) for p in out.iterdir() if p.is_file()})
    with (out/'verification.json').open('x') as f:json.dump(record,f,indent=2)
    print(json.dumps(dict(passed=True,cases=[{k:c[k] for k in ('name','maximum_displacement_error_m','maximum_pressure_error_kpa')} for c in reports])),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('solver',type=Path);p.add_argument('output',type=Path);a=p.parse_args();run(a.solver,a.output.resolve())
