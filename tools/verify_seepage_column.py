"""Input-only constant-modulus column theory for flow and excavation paths."""
import argparse
import json
from pathlib import Path
import numpy as np
import fem_reference as ref


def verify(run):
    mesh=json.loads((run/'ground/mesh.json').read_text());data=json.loads((run/'stage-results.json').read_text());assert data['complete']
    assert all('mohr_coulomb' in m for m in data['ground']['materials'])
    xyz=np.asarray(mesh['nodes']);errors=[]
    for index,frame in enumerate(data['frames']):
        used=np.unique(ref.cells(mesh)[np.asarray(frame['state']['active_elements'])]);z=xyz[used,2]
        expected=[np.zeros_like(z),4*(z-z*z/2)/12000,-4*(z-z*z/2)/12000,np.zeros_like(z),5*z/12000,(z+4*z*z)/12000][index]
        error=float(np.max(abs(np.asarray(frame['state']['displacement_m'])[used,2]-expected)));assert error<2e-11,error
        errors.append({'stage':frame['name'],'maximum_nodal_displacement_error_m':error})
    with (run/'column-theory.json').open('x') as stream:json.dump({'schema':'soil-fem-seepage-column-theory/1','passed':True,'frames':errors},stream,indent=2)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run',type=Path);verify(p.parse_args().run.resolve())
