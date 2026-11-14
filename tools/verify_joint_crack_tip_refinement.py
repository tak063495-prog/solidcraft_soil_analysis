"""Small local h series: displacement, energy and stress away from the tip.

No singular tip-point stress comparison or general spatial-convergence claim.
"""
import argparse
import copy
import json
from pathlib import Path
import numpy as np
import fem_reference as ref
from joint_tip_fixture import fixture
from verify_joint_crack_tip import ROOT, write, solve, reassembly, compare_reference


def main(solver, output):
    output.mkdir(parents=True,exist_ok=False)
    grids=[dict(x=(0.,1.,2.),y=(0.,1.),z=(0.,.5,1.)),
           dict(x=(0.,.5,1.,1.5,2.),y=(0.,1.),z=(0.,.25,.5,.75,1.)),
           dict(x=(0.,.5,.75,1.,1.25,1.5,2.),y=(0.,.5,1.),z=(0.,.25,.5,.75,1.))]
    ground=json.loads((ROOT/'examples/ground-joint-crack-tip.json').read_text());template=json.loads((ROOT/'examples/stages-joint-crack-tip.json').read_text())
    rows=[]
    for level,grid in enumerate(grids):
        _,_,mesh,cfg=fixture(**grid);case=copy.deepcopy(template);case['joint_system']=cfg;cfg['maximum_increment_m']=.05
        case['settings']['maximum_free_dofs']=2000
        for c,v in zip(case['boundary']['constraints'][-3:],[1e-5,0.,-1e-4]):c['value_m']=v
        case['stages']=[dict(kind='boundary',name='compress',target_factor=1.,initial_increment=.5)]
        run,_=solve(solver,output/f'level-{level}',mesh,ground,case)
        assembly=reassembly(mesh,ground,run['case'],run);independent=compare_reference(mesh,ground,run['case'],run)
        state=run['last_committed'];u=np.asarray(state['displacement_m']);xyz=np.asarray(mesh['nodes']);upper=set(ref.cells(mesh)[np.asarray(mesh['region_ids'])==2].flat)
        monitor=next(n for n in upper if np.linalg.norm(xyz[n]-[0.,.5,.5])<1e-12)
        energy=0.;stress_sum=np.zeros(6);volume=0.;start=0
        mat={m['id']:m['mohr_coulomb'] for m in ground['materials']}
        for e in range(len(mesh['tetrahedra'])):
            q=ref.volume(mesh,e);m=mat[mesh['material_ids'][e]];strain=np.einsum('qai,i->qa',q['b'],u[q['nodes']].ravel());d=ref.elastic_tensor(m['young_kpa'],m['poisson'])
            energy+=.5*np.einsum('qa,ab,qb,q->',strain,d,strain,q['weight'])
            stress=np.array([p['effective_stress_compression_positive_kpa'] for p in state['material_states'][start:start+len(q['weight'])]])
            # The observation domain x >= 1.5 m stays at least 0.5 m from
            # the x=1, z=.5 tip line and has mesh-aligned boundaries at all h.
            if np.min(xyz[mesh['tetrahedra'][e],0])>=1.5-1e-12:
                stress_sum+=np.einsum('qa,q->a',stress,q['weight']);volume+=q['weight'].sum()
            start+=len(q['weight'])
        # Coarsest cells cross x=1.5: integrate the same domain independently
        # with a Tet10 polynomial evaluated on an observation submesh.
        if volume==0:
            from joint_tip_fixture import conformal
            observation=conformal(x=(1.5,2.),y=(0.,1.),z=(0.,.5,1.))
            for e in range(len(observation['tetrahedra'])):
                oq=ref.volume(observation,e,independent=True)
                for point,weight in zip(oq['position'],oq['weight']):
                    for ce,t in enumerate(mesh['tetrahedra']):
                        corners=xyz[t];r=np.linalg.solve((corners[1:]-corners[0]).T,point-corners[0]);l=np.r_[1-r.sum(),r]
                        if l.min()>=-1e-10:
                            _,gradient=ref.basis(r[None,:],True);j=xyz[mesh['tetrahedra10'][ce]].T@gradient[0]
                            g=gradient[0]@np.linalg.inv(j);disp=u[mesh['tetrahedra10'][ce]];a=disp.T@g
                            strain=np.array([a[0,0],a[1,1],a[2,2],a[0,1]+a[1,0],a[1,2]+a[2,1],a[2,0]+a[0,2]])
                            m=mat[mesh['material_ids'][ce]];s=np.array([100.,100.,100.,0.,0.,0.])-ref.elastic_tensor(m['young_kpa'],m['poisson'])@strain
                            stress_sum+=weight*s;volume+=weight;break
                    else:raise AssertionError('observation point not contained')
        assert abs(volume-.5)<1e-12
        energy+=run['frames'][-1]['joint_results']['elastic_energy_kn_m']
        rows.append(dict(level=level,elements=len(mesh['tetrahedra']),nodes=len(mesh['nodes']),upper_crack_mouth_displacement_m=u[monitor].tolist(),elastic_increment_energy_kn_m=float(energy),away_domain_mean_stress_kpa=(stress_sum/volume).tolist(),assembly=assembly,independent=independent))
    differences=[]
    for a,b in zip(rows,rows[1:]):
        differences.append(dict(levels=[a['level'],b['level']],displacement_difference_m=float(np.linalg.norm(np.array(a['upper_crack_mouth_displacement_m'])-b['upper_crack_mouth_displacement_m'])),energy_difference_kn_m=abs(a['elastic_increment_energy_kn_m']-b['elastic_increment_energy_kn_m']),away_stress_difference_kpa=float(np.linalg.norm(np.array(a['away_domain_mean_stress_kpa'])-b['away_domain_mean_stress_kpa']))))
    write(output/'measurements.json',dict(schema='soil-fem-joint-crack-tip-local-h/1',passed=True,scope='three small local h meshes; original accuracy gates pass at each h; differences are observations, no general convergence threshold invented',meshes=rows,differences=differences))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();main(a.solver.resolve(),a.output.resolve())
