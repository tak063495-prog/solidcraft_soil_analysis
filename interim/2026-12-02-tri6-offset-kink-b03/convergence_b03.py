"""Common physical probes, true recovery fields and work-conjugate support restriction."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import argparse,json
from pathlib import Path
import numpy as np
from b03_cases import read,write
from convergence import shape_at,differences
from plate6_offset_reference import System
from plate6_reference import SAMPLES
from work_audit import calculate
import fem_reference as ref


def tensor(v,axes):
    return (axes.T@np.array([[v[0],v[2],0.],[v[2],v[1],0.],[0.,0.,0.]])@axes).ravel()


def sample(case,coarse_case,kind):
    m,g,c=[read(case/(k+'.json')) for k in ('mesh','ground','stages')]
    coarse=read(coarse_case/'mesh.json');cc=read(coarse_case/'stages.json')
    model=System(m,c['beam_system']);cm=System(coarse,cc['beam_system']);xyz=np.asarray(m['nodes'])
    ancestry=np.arange(len(coarse['tetrahedra']))
    for a in read(case/'refinement.json'):ancestry=ancestry[np.asarray(a['parent_of_element'])]
    soil_ops=[shape_at(m,q,np.flatnonzero(ancestry==ei)) for ei in range(len(coarse['tetrahedra'])) for q in ref.volume(coarse,ei)['position']]
    materials={v['id']:v['mohr_coulomb'] for v in g['materials']}
    dmat={k:ref.elastic_tensor(v['young_kpa'],v['poisson']) for k,v in materials.items()}
    supports=m['node_sets']['bottom'];restriction=np.zeros((len(coarse['nodes']),len(supports)))
    for j,n in enumerate(supports):
        _,ids,shape,_=shape_at(coarse,xyz[n]);restriction[ids,j]=shape
    plate_ops=[];orientation_metadata=[]
    for ci,(cn,_,ck) in enumerate(cm.plate_data):
        for ll in SAMPLES:
            query=np.asarray(ll)@cm.positions[cn][:3]
            for pi,(nodes,section,kernel) in enumerate(model.plate_data):
                points=model.positions[nodes][:3];edge=(points[1:]-points[0]).T
                r=np.linalg.lstsq(edge,query-points[0],rcond=None)[0];l=np.r_[1-r.sum(),r]
                if l.min()>=-1e-8 and np.linalg.norm(edge@r-(query-points[0]))<=1e-8:
                    orientation=1. if kernel.normal@ck.normal>0 else -1.
                    plate_ops.append((pi,nodes,kernel,l,orientation))
                    orientation_metadata.append(dict(coarse_plate_id=cm.plates[ci]['id'],fine_plate_id=model.plates[pi]['id'],
                        attachment_query_m=query.tolist(),fine_barycentric=l.tolist(),fine_normal=kernel.normal.tolist(),
                        coarse_normal=ck.normal.tolist(),fine_e_m=kernel.e,coarse_e_m=ck.e,oriented_section_factor=orientation))
                    break
            else:raise ValueError(('no owning planar face',query))
    work=calculate(case,kind);rows={}
    source=case/('native/accepted.jsonl' if kind=='native' else 'independent.jsonl')
    with source.open() as f:
        for index,line in enumerate(f):
            r=json.loads(line);key=(r['stage_index'],r['progress'])
            if r['progress'] not in (.25,.5,.75,1.):continue
            s=r['state'] if kind=='native' else r;u=np.asarray(s['displacement_m']);bs=s['beams'];report=r['beam_results']
            stress=[]
            for ei,ids,shape,grad in soil_ops:
                du=u[ids].T@grad;eps=np.array([du[0,0],du[1,1],du[2,2],du[0,1]+du[1,0],du[1,2]+du[2,1],du[0,2]+du[2,0]])
                stress.append(-dmat[m['material_ids'][ei]]@eps)
            fields=dict(soil_common_nodes_m=u[:len(coarse['nodes'])].ravel(),soil_probe_stress_kpa=np.ravel(stress),
                        global_F_dot_d_kn_m=np.array([work[index]['global_F_dot_d_kn_m']]),
                        cumulative_total_work_kn_m=np.array([work[index]['cumulative_total_work_kn_m']]))
            d=np.c_[bs['node_displacement_m'],bs['node_rotation_rad']]
            if kind=='native':
                internal=np.asarray(r['soil_internal_kn']);internal[np.asarray(r['interface_dofs'],int)]+=r['interface_internal_kn']
                reaction=(internal-np.asarray(s['applied_nodal_load_kn'])).reshape(-1,3)[supports]
            else:reaction=np.asarray(r['support_reaction_by_node_kn'])[supports]
            fields['support_generalized_reactions_kn']=(restriction@reaction)[coarse['node_sets']['bottom']].ravel()
            if any(bs['active_plates']):
                vals={k:[] for k in ('plate_attachment_u_m','plate_rotation_rad','plate_membrane_global_kn_m','plate_bending_global_kn','plate_attachment_moment_global_kn','plate_shear_global_kn_m')}
                for pi,nodes,kernel,l,orientation in plate_ops:
                    kin=d[nodes].ravel()-np.asarray(bs['plate_reference_global_kinematics'][pi]);local=kernel.neutral.transform@kernel.map@kin
                    bm,bb,bq,_=kernel.neutral.operators(l)
                    nf=kernel.neutral.membrane@(bm@local);moment=kernel.neutral.bending@(bb@local);shear=kernel.neutral.shear*(bq@local)
                    if kind=='native':
                        plate=next(p for p in report['plates'] if p['id']==model.plates[pi]['id'])
                        xy=np.asarray([np.asarray(v['barycentric'])@kernel.neutral.xy[:3] for v in plate['samples']]);fit=np.c_[np.ones(len(xy)),xy];at=np.r_[1.,l@kernel.neutral.xy[:3]]
                        nf=at@np.linalg.lstsq(fit,[v['membrane_force_kn_m'] for v in plate['samples']],rcond=None)[0]
                        moment=at@np.linalg.lstsq(fit,[v['bending_moment_kn'] for v in plate['samples']],rcond=None)[0]
                        basis=np.vstack([kernel.neutral.basis(*v) for v in xy]);coef=np.linalg.lstsq(basis,np.ravel([v['transverse_shear_force_kn_m'] for v in plate['samples']]),rcond=None)[0]
                        shear=kernel.neutral.basis(*(l@kernel.neutral.xy[:3]))@coef
                    shape=kernel.neutral.shapes(*(l@kernel.neutral.xy[:3]))[0]
                    vals['plate_attachment_u_m'].extend(shape@d[nodes,:3]);vals['plate_rotation_rad'].extend(shape@d[nodes,3:])
                    vals['plate_membrane_global_kn_m'].extend(tensor(nf,kernel.axes));vals['plate_bending_global_kn'].extend(orientation*tensor(moment,kernel.axes))
                    vals['plate_attachment_moment_global_kn'].extend(orientation*tensor(moment-kernel.e*nf,kernel.axes));vals['plate_shear_global_kn_m'].extend(orientation*(kernel.axes.T@np.r_[shear,0.]))
                fields.update({k:np.asarray(v) for k,v in vals.items()})
            rows[key]=fields
    return rows,dict(plate_probe_orientation_factors=[v[4] for v in plate_ops],
        reference='Each original coarse face oriented normal. Fine input n/e preserved; membrane tensor is invariant, bending/attachment moment and transverse shear are oriented section resultants and multiply by sign(n_f.n_c). No input normalization.',
        plate_probe_count=len(plate_ops),soil_probe_count=len(soil_ops),orientation_metadata=orientation_metadata)


def run(formal,out):
    result={};works={};maxwork={}
    cases=read(formal/'measurements.json')['cases']
    for case in cases:
        folder=formal/case['case'];n=calculate(folder,'native');p=calculate(folder,'independent');d=[]
        for aa,bb in zip(n,p,strict=True):
            assert (aa['stage_index'],aa['progress'])==(bb['stage_index'],bb['progress'])
            delta={k:aa[k]-bb[k] for k in aa if k not in ('stage_index','progress')};d.append(delta)
            for k,v in delta.items():maxwork[k]=max(maxwork.get(k,0.),abs(v))
        works[case['case']]=dict(native=n,independent=p,differences=d)
    for family in sorted({c['kind'] for c in cases}-{'planar'}):
        result[family]={}
        for series,names in [('space',[f'{family}-h{level}-s4' for level in range(max(c['level'] for c in cases if c['kind']==family)+1)]),('load',[f'{family}-h0-s{s}' for s in (4,8,16)])]:
            result[family][series]={}
            for kind in ('native','independent'):
                sampled=[sample(formal/n,formal/f'{family}-h0-s4',kind) for n in names]
                values=[r[0] for r in sampled]
                pairs=[differences(a,b) for a,b in zip(values,values[1:])]
                trend={}
                if len(pairs)==2:
                    for k,first in pairs[0]['absolute_l2_all_points'].items():
                        second=pairs[1]['absolute_l2_all_points'][k];trend[k]=dict(first=first,second=second,decreasing=second<first,observed_order=float(np.log2(first/second)) if first>0 and second>0 else None)
                result[family][series][kind]=dict(cases=names,common_accepted_points=len(values[0]),observation_metadata=[r[1] for r in sampled],pairs=pairs,trends=trend,
                    limitations='2 spatial levels: no spatial convergence order; purely elastic load fields can be roundoff dominated; no new gate')
    write(out,dict(passed=False,b03_complete=False,series=result))
    write(out.with_name('work-audit.json'),dict(passed=False,b03_complete=False,numerical_acceptance=None,cases=works,maximum_absolute_difference=maxwork,
        definition='Full soil RHS minus attached structure loads, plus explicit structure loads; accepted-endpoint trapezoid; own installation reference for newly active nodes. Construction holding included.',
        note='Derived from fresh solved fields, not additional solves or exact continuous-path external work'))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('formal',type=Path);p.add_argument('output',type=Path);a=p.parse_args();run(a.formal.resolve(),a.output.resolve())
