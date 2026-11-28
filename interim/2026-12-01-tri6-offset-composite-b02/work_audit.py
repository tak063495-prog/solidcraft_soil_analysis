"""Work from the two fresh solved histories; postprocessing, not another solve.

Global F.d and trapezoidal accepted-endpoint path work are kept separate.
Birth endpoints use each history's own installation reference, never the other
solver's state. No work tolerance or physical continuous-path claim is added.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from cases import read,write
from plate6_offset_reference import System


def calculate(folder,kind):
    m,c=read(folder/'mesh.json'),read(folder/'stages.json');model=System(m,c['beam_system'])
    nd=3*len(m['nodes']);nb=6*len(model.nodes)
    old_u=np.zeros(nd);old_d=np.zeros(nb);old_soil=np.zeros(nd);old_beam=np.zeros(nb);old_active=np.zeros(len(model.nodes),bool)
    cumulative=np.zeros(2);rows=[]
    source=folder/('native/accepted.jsonl' if kind=='native' else 'independent.jsonl')
    with source.open() as f:
        for line in f:
            record=json.loads(line);s=record['state'] if kind=='native' else record
            bs=s['beams'];factor=s.get('boundary_factor')
            if factor is None:
                # The independent accepted path is constructed solely from its
                # input stage targets and own progress, not native observations.
                si=record['stage_index'];factor=0.
                for stage in c['stages'][:max(si-1,0)]:factor=stage.get('target_factor',factor)
                if si and c['stages'][si-1]['kind']=='boundary':factor+=record['progress']*(c['stages'][si-1]['target_factor']-factor)
            u=np.asarray(s['displacement_m']).ravel();d=np.c_[bs['node_displacement_m'],bs['node_rotation_rad']].ravel()
            _,_,fb,eff,*_=model.assemble(bs,factor)
            soil=np.asarray(s['applied_nodal_load_kn']).copy()
            for i,n in enumerate(model.attach):
                if n is not None:soil[3*n:3*n+3]-=eff[i]
            active=model.active_nodes(bs);birth=active & ~old_active;previous=old_d.copy()
            if birth.any():
                refs={}
                for i,(nodes,*_) in enumerate(model.plate_data):
                    if bs['active_plates'][i]:
                        for j,n in enumerate(nodes):refs[n]=bs['plate_reference_global_kinematics'][i][6*j:6*j+6]
                for i,(nodes,*_) in enumerate(model.data):
                    if bs['active_members'][i]:
                        for j,n in enumerate(nodes):refs[n]=bs['reference_global_kinematics'][i][6*j:6*j+6]
                for n in np.flatnonzero(birth):previous[6*n:6*n+6]=refs[n]
            increment=np.array([.5*(soil+old_soil)@(u-old_u),.5*(fb+old_beam)@(d-previous)])
            cumulative+=increment
            rows.append(dict(stage_index=record['stage_index'],progress=record['progress'],
                             soil_increment_work_kn_m=float(increment[0]),structure_increment_work_kn_m=float(increment[1]),
                             cumulative_soil_work_kn_m=float(cumulative[0]),cumulative_structure_work_kn_m=float(cumulative[1]),
                             cumulative_total_work_kn_m=float(cumulative.sum()),global_F_dot_d_kn_m=float(soil@u+fb@d)))
            old_u=u.copy();old_d=d.copy();old_soil=soil.copy();old_beam=fb.copy();old_active=active.copy()
    return rows


def run(formal,output):
    result={};maximum={}
    for case in read(formal/'measurements.json')['cases']:
        folder=formal/case['case'];native=calculate(folder,'native');independent=calculate(folder,'independent');delta=[]
        for n,p in zip(native,independent,strict=True):
            assert (n['stage_index'],n['progress'])==(p['stage_index'],p['progress'])
            row={k:n[k]-p[k] for k in n if k not in ('stage_index','progress')}
            for k,v in row.items():maximum[k]=max(maximum.get(k,0.),abs(v))
            delta.append(row)
        result[case['case']]=dict(native=native,independent=independent,differences=delta)
    write(output,dict(passed=False,numerical_acceptance=None,b02_complete=False,cases=result,maximum_absolute_difference=maximum,
        definition='Accepted-endpoint trapezoid of full soil RHS with interface loads removed, plus explicit structure loads; new nodes begin at each history own birth reference. Includes construction holding. Global F.d is separate.',
        limitations='Derived observations from the fresh solved fields, not additional solves. Birth/removal quadrature is a discrete diagnostic, not exact continuous external work or a new energy acceptance criterion.'))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('formal',type=Path);p.add_argument('output',type=Path);a=p.parse_args();run(a.formal.resolve(),a.output.resolve())
