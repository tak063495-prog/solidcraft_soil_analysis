"""Independent ray-root, analytic and bounded native wet-face comparison."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from hydrostatic_reference import face_force
from verify_checkpoint import invoke


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,d):
    with p.open('x') as f:json.dump(d,f,indent=2)


def verify(solver,output):
    output.mkdir(parents=True,exist_ok=False);natural=[(0,0),(1,0),(0,1),(.5,0),(.5,.5),(0,.5)]
    cases={'planar':([[0.,r,s] for r,s in natural],.37),
        'curved':([[0.,r,.2*s+.8*s*s] for r,s in natural],.37),
        'island':([[r,s,(r-1/3)**2+(s-1/3)**2-.01] for r,s in natural],0.),
        'spatial':([[r+.12*r*s,s+.08*r*(1-r-s),.8*s+.15*r*r-.2*r*s] for r,s in natural],.31)}
    rows={}
    for name,(points,h) in cases.items():
        reference,diagnostics=face_force(points,h,10.);write(output/(name+'-reference.json'),{'nodal_force_kn':reference.tolist(),**diagnostics})
        levels=[]
        for tol in (1e-7,1e-10):
            label=f'{name}-{tol:.0e}';inp=output/(label+'-input.json');out=output/(label+'-native.json')
            spec={'schema':'soil-fem-hydrostatic-face/1','points_m':points,'water_level_m':h,'water_unit_weight_kn_m3':10.,
                'integration':{'relative_tolerance':tol,'absolute_tolerance_kn':tol/10,'maximum_cells':50000}}
            write(inp,spec);invoke([str(solver),'hydrostatic-face',str(inp),str(out)],output/(label+'.log'));native=json.loads(out.read_text())
            error=float(np.linalg.norm(np.asarray(native['nodal_force_kn'])-reference));bound=native['domain_error_bound_kn']
            assert error<=bound+5e-12,(name,tol,error,bound)
            assert bound<=native['permitted_domain_error_kn']
            levels.append({'relative_tolerance':tol,'absolute_nodal_force_l2_error_kn':error,'domain_error_bound_kn':bound,'visited_cells':native['visited_cells']})
        if name!='planar':assert levels[-1]['absolute_nodal_force_l2_error_kn']<levels[0]['absolute_nodal_force_l2_error_kn']/10
        rows[name]={'levels':levels,'reference':diagnostics}
    bad=spec.copy();bad['integration']={**spec['integration'],'maximum_cells':1};inp=output/'unresolved-input.json';out=output/'unresolved-result.json';write(inp,bad)
    invoke([str(solver),'hydrostatic-face',str(inp),str(out)],output/'unresolved.log',1);assert not out.exists()
    bad={'schema':'soil-fem-hydrostatic-face/1','points_m':[[0.,0.,0.],[1e100,0.,0.],[0.,1e100,0.]],'water_level_m':-1.,'water_unit_weight_kn_m3':10.}
    inp=output/'overflow-input.json';out=output/'overflow-result.json';write(inp,bad);before=sha(inp)
    invoke([str(solver),'hydrostatic-face',str(inp),str(out)],output/'overflow.log',1);assert not out.exists() and sha(inp)==before
    report={'schema':'soil-fem-independent-hydrostatic-faces/1','passed':True,'solver_sha256':sha(solver),'cases':rows,'negative_cases':['unresolved_budget','nonfinite_area_scale'],
        'scope':'Native Bernstein-bound reference-patch clipping versus independent inner quadratic roots and adaptive outer Gauss integration. Actual Tri6 mapping and all nodal forces; planar, curved, interior wet island and spatial surface; unresolved budget rejects.'}
    write(output/'verification.json',report);print(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();verify(a.solver.resolve(),a.output_dir.resolve())
