"""Project existing worst-case rotation differences into original element frames."""
import argparse
from pathlib import Path
import numpy as np
from b01follow_common import mm, read, rotate, write


def observe(formal,measurements,output):
    m=read(measurements);kernel={r['id']:r for r in read(formal/'kernel-responses.json')};rows=[]
    for selected in m['precision_accounting']['cases']:
        name,field=selected['id'],selected['field'];base=name.split('--')[0]+'--identity'
        b=read(formal/'assemblies'/(base+'.json'));r=read(formal/'assemblies'/(name+'.json'))
        db=np.array(next(s['displacement_rotation'] for s in b['solves'] if s['name']==field)).reshape(-1,6)
        dr=rotate(np.array(next(s['displacement_rotation'] for s in r['solves'] if s['name']==field)).reshape(-1,6),np.array(r['rotation_matrix']).T)
        delta=dr-db;local=[]
        for e,nodes in enumerate(b['elements']):
            axes=np.array(kernel[base+'--e'+str(e)]['axes']);local.extend(mm(delta[nodes,3:],axes.T).tolist())
        local=np.array(local)
        rows.append(dict(id=name,field=field,local_tangent_rotation_max_rad=float(np.max(np.abs(local[:,:2]))),
             local_normal_rotation_max_rad=float(np.max(np.abs(local[:,2]))),
             interpretation='Projection uses each original element frame only; no averaged nodal director introduced. Normal rotation is the drilling direction, not a new director state.'))
    write(output,dict(diagnostic_only=True,passed=False,b01_complete=False,cases=rows))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--formal',type=Path,required=True)
    p.add_argument('--measurements',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();observe(a.formal.resolve(),a.measurements.resolve(),a.output.resolve())
