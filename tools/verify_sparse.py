"""Independent Tet4 assembly/dense solve plus measured sparse grid scaling."""
import argparse
import csv
import ctypes
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import numpy as np


def monitored(command,log):
    """Observe the working set of this child process only; no external jobs."""
    peak=None;started=time.perf_counter()
    if os.name=='nt':
        from ctypes import wintypes
        class Counters(ctypes.Structure):
            _fields_=[('cb',wintypes.DWORD),('PageFaultCount',wintypes.DWORD)]+[
                (name,ctypes.c_size_t) for name in ('PeakWorkingSetSize','WorkingSetSize',
                    'QuotaPeakPagedPoolUsage','QuotaPagedPoolUsage','QuotaPeakNonPagedPoolUsage',
                    'QuotaNonPagedPoolUsage','PagefileUsage','PeakPagefileUsage')]
        query=ctypes.WinDLL('psapi',use_last_error=True).GetProcessMemoryInfo
        query.argtypes=[wintypes.HANDLE,ctypes.POINTER(Counters),wintypes.DWORD];query.restype=wintypes.BOOL
    with log.open('x',encoding='utf-8') as stream:
        stream.write('Command: '+repr(command)+'\n');stream.flush()
        process=subprocess.Popen(command,stdout=stream,stderr=subprocess.STDOUT)
        while True:
            if os.name=='nt':
                counters=Counters();counters.cb=ctypes.sizeof(counters)
                if query(int(process._handle),ctypes.byref(counters),counters.cb):
                    peak=max(peak or 0,counters.PeakWorkingSetSize)
            code=process.poll()
            if code is not None:break
            time.sleep(.02)
        stream.write(f'\nExit code: {code}\n')
    if code:raise RuntimeError(f'child failed ({code}): {log}')
    return {'wall_seconds_including_io':time.perf_counter()-started,
            'observed_peak_working_set_bytes':peak,'memory_method':'Windows GetProcessMemoryInfo child-process peak working set' if os.name=='nt' else 'unavailable'}


def independent_small(report):
    case=report['case'];n=case['subdivisions'];w=n+1
    xyz=np.array([[x/n,y/n,z/n] for z in range(w) for y in range(w) for x in range(w)])
    lookup=np.full(xyz.shape,-1,dtype=int);expected=[]
    for i,(x,y,z) in enumerate(xyz):
        fixed=[x==0 or x==1,y==0 or y==1,z==0]
        u=[x*(1-x)*(1+y)*z,y*(1-y)*(1+z)*x,z*(1+x*y)]
        for a in range(3):
            if not fixed[a]:lookup[i,a]=len(expected);expected.append(u[a])
    free=len(expected);a=np.zeros((free,free));mu=case['young_kpa']/(2*(1+case['poisson']))
    lam=case['young_kpa']*case['poisson']/((1+case['poisson'])*(1-2*case['poisson']))
    d=np.zeros((6,6));d[:3,:3]=lam
    for i in range(3):d[i,i]+=2*mu;d[i+3,i+3]=mu
    for z in range(n):
        for y in range(n):
            for x in range(n):
                ids=[x+dx+w*(y+dy+w*(z+dz)) for dx,dy,dz in
                     ((0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1))]
                for corners in ((0,1,2,6),(0,2,3,6),(0,3,7,6),(0,7,4,6),(0,4,5,6),(0,5,1,6)):
                    tet=np.array([ids[i] for i in corners]);p=xyz[tet];j=(p[1:]-p[0]).T
                    v=np.linalg.det(j)/6;g=np.linalg.inv(j);g=np.vstack((-g.sum(axis=0),g))
                    b=np.zeros((6,12))
                    for i,(gx,gy,gz) in enumerate(g):
                        k=3*i;b[0,k]=gx;b[1,k+1]=gy;b[2,k+2]=gz
                        b[3,k]=gy;b[3,k+1]=gx;b[4,k+1]=gz;b[4,k+2]=gy;b[5,k]=gz;b[5,k+2]=gx
                    ke=v*b.T@d@b;dofs=lookup[tet].reshape(-1);keep=np.flatnonzero(dofs>=0)
                    a[np.ix_(dofs[keep],dofs[keep])]+=ke[np.ix_(keep,keep)]
    system=report['system'];rust=np.zeros_like(a)
    for i in range(free):
        lo,hi=system['row_ptr'][i:i+2];rust[i,np.array(system['columns'][lo:hi])]=system['values'][lo:hi]
    matrix_error=np.linalg.norm(a-rust)/np.linalg.norm(a)
    expected=np.array(expected);rhs=a@expected
    assert np.allclose(rhs,system['rhs'],rtol=1e-12,atol=1e-10)
    assert np.allclose(expected,system['expected_solution'],rtol=1e-14,atol=1e-14)
    exact=np.linalg.solve(a,rhs);x=np.array(system['computed_solution'])
    difference=np.linalg.norm(exact-x)/np.linalg.norm(exact)
    assert matrix_error<1e-13 and difference<1e-8
    return {'matrix_relative_difference':float(matrix_error),'solution_relative_difference_from_numpy_dense':float(difference)}


def verify(solver,output):
    output.mkdir(parents=True,exist_ok=False);rows=[]
    options={'restart':60,'maximum_iterations':3000,'pivot_tolerance':1e-14}
    for n in (3,6,12,26):
        case={'subdivisions':n,'young_kpa':1000.,'poisson':.3,'relative_linear_tolerance':1e-10,'options':options,'export_system':n==3}
        case_path=output/f'grid-{n}.json';report_path=output/f'grid-{n}-result.json';case_path.write_text(json.dumps(case,indent=2))
        timing=monitored([str(solver),'benchmark-tet4',str(case_path),str(report_path)],output/f'grid-{n}.log')
        report=json.loads(report_path.read_text());linear=report['linear']
        assert report['tetrahedra']==6*n**3 and linear['true_relative_residual']<=1e-10 and report['relative_solution_error']<1e-7
        row={'subdivisions':n,'tetrahedra':report['tetrahedra'],'free_dofs':report['free_dofs'],
             'csr_entries':linear['stored_nonzeros'],'symbolic_seconds':report['symbolic_seconds'],
             'assembly_seconds':report['assembly_seconds'],'linear_solve_seconds':linear['elapsed_seconds'],
             'iterations':linear['iterations'],'true_linear_residual':linear['true_relative_residual'],
             'relative_solution_error':report['relative_solution_error'],'matrix_storage_bytes':linear['matrix_storage_bytes'],
             'workspace_storage_bytes':linear['workspace_storage_bytes'],**timing}
        if n==3:small=independent_small(report)
        rows.append(row);print(json.dumps(row),flush=True)
    assert rows[-1]['tetrahedra']>=100_000
    summary={'schema':'soil-fem-sparse-verification/1','scope':'linear Tet4 FEM backend; not full nonlinear 100000-element excavation performance',
             'independent_small':small,'measurements':rows,'solver_sha256':hashlib.sha256(solver.read_bytes()).hexdigest(),
             'inputs_and_reports_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in output.glob('grid-*.json')}}
    (output/'sparse-verification.json').write_text(json.dumps(summary,indent=2))
    with (output/'scaling.csv').open('x',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print('Independent sparse FEM verification passed',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,default=Path('target/release/soil-fem-3d.exe'));p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();verify(a.solver.resolve(),a.output_dir.resolve())
