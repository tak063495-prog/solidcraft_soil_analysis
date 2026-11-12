"""Independent free structural patch diagnostic; not the soil-coupled stiffness or acceptance."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import argparse
from pathlib import Path
import numpy as np
from b01_cases import read,write
from plate6_offset_reference import System


def run(formal,out):
    rows=[]
    for rec in read(formal/'measurements.json')['cases']:
        if rec['n'] != 2:continue
        f=formal/rec['case'];m,c=read(f/'mesh.json'),read(f/'stages.json');model=System(m,c['beam_system'])
        state=model.initial_state();state['active_plates']=[True]*len(model.plates);k=model.assemble(state,0.)[0]
        modes=[]
        for a in range(6):
            d=np.zeros((len(model.nodes),6))
            if a<3:d[:,a]=1.
            else:d[:,3+a-3]=1.;d[:,:3]=np.cross(np.eye(3)[a-3],model.positions)
            modes.append(dict(basis=a,maximum_force_residual=float(np.max(np.abs(k@d.ravel()))),
                normalized_force_residual=float(np.linalg.norm(k@d.ravel())/max(np.linalg.norm(k,2)*np.linalg.norm(d),1.)),energy_kn_m=float(d.ravel()@k@d.ravel()/2)))
        adjacency={i:{i} for i in range(len(model.nodes))}
        for nodes,_,_ in model.plate_data:
            for i in nodes:adjacency[i].update(nodes)
        unvisited=set(adjacency);components=[]
        while unvisited:
            todo=[next(iter(unvisited))];found=set()
            while todo:
                i=todo.pop()
                if i in found:continue
                found.add(i);todo.extend(adjacency[i]-found)
            unvisited-=found;components.append(sorted(found))
        eigen=np.linalg.eigvalsh(k);threshold=1e-10*eigen.max()
        rows.append(dict(case=rec['case'],components=components,structural_components=len(components),
                         expected_free_modes_from_explicit_components=6*len(components),
                         rank_threshold_from_existing_offset_element_test=threshold,positive_rank=int((eigen>threshold).sum()),
                         modes_below_existing_rank_threshold=int((eigen<=threshold).sum()),eigenvalues=eigen.tolist(),global_six_rigid_basis=modes))
    write(out,dict(passed=False,b01_complete=False,accepted=None,cases=rows,
          scope='Input-only independent free structural patches, before soil constraints. Separate IDs do not enforce a rigid corner; soil coupling stabilizes translations separately. Existing rank diagnostic, no new acceptance gate.'))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('formal',type=Path);p.add_argument('output',type=Path);a=p.parse_args();run(a.formal.resolve(),a.output.resolve())
