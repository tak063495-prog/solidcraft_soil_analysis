"""Post-solve GP correspondence audit; this is not another fresh solve."""
import argparse
from pathlib import Path
import numpy as np
from study import read, write
import fem_reference as fem
from joint_reference import System


def audit(folder):
    mesh=read(folder/'mesh.json');case=read(folder/'case.json');native=read(folder/'native/run.json')
    contact=System(mesh,case['joint_system']);tol=case['joint_system']['geometric_tolerance_m']
    errors={};ids=0
    def add(name,a,b):
        a,b=np.asarray(a),np.asarray(b)
        assert a.shape==b.shape
        errors[name]=max(errors.get(name,0.),float(np.max(np.abs(a-b),initial=0.)))
    assert len(native['joint_integration_points'])==len(contact.points)
    for a,b in zip(native['joint_integration_points'],contact.points):
        for key in ('joint_id','face_index','nodes_a','nodes_b','owner_a','owner_b'):
            assert a[key]==b[key],key
        ids+=1
        for key in ('position_m','normal','weight_m2','shape'):add('contact_'+key,a[key],b[key])
    offset=0
    for e in range(len(mesh['tetrahedra'])):
        q=fem.volume(mesh,e)
        for k in range(len(q['weight'])):
            a=native['integration_points'][offset+k]
            assert (a['element_index'],a['local_index'])==(e,k)
            for left,right in (('natural_coordinates','natural'),('position_m','position'),('weight_m3','weight'),('shape_functions','shape')):
                add('material_'+left,a[left],q[right][k])
        offset+=len(q['weight'])
    assert offset==len(native['integration_points'])
    assert errors['contact_position_m']<=tol and errors['material_position_m']<=tol
    return dict(postsolve_audit_not_fresh_solve=True,contact_ids_exact=ids,material_ids_exact=offset,
        position_gate_m=tol,positions_within_existing_geometry_tolerance=True,maximum_absolute_differences=errors,
        shape_weight_normal_differences_are_observations_only=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True);a=p.parse_args()
    write(a.run/'geometry-audit.json',{name:audit(a.run/name) for name in ('transfer_force','stress_free','stress_free-cycle')})


if __name__=='__main__':main()
