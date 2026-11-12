"""Common physical probes and virtual-work restriction; no new acceptance gate."""
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import argparse
import json
from pathlib import Path
import numpy as np
from cases import read, write
import fem_reference as ref
from plate6_offset_reference import System
from plate6_reference import SAMPLES


def bary(points, query):
    r=np.linalg.solve((points[1:4]-points[0]).T, query-points[0])
    return r


def shape_at(mesh, query, choices=None):
    xyz=np.asarray(mesh['nodes']); cells=ref.cells(mesh)
    for ei in range(len(cells)) if choices is None else choices:
        cell=cells[ei];points=xyz[cell];r=bary(points,query)
        if min(1-r.sum(),*r) < -1e-8:continue
        shape,gradient=ref.basis(r[None],True)
        grad=gradient[0]@np.linalg.inv((points[1:4]-points[0]).T)
        return ei,cell,shape[0],grad
    raise ValueError(('no containing parent/child',query))


def tensor(v,axes):
    return (axes.T@np.array([[v[0],v[2],0.],[v[2],v[1],0.],[0.,0.,0.]])@axes)[:2,:2].ravel()


def sample(case, coarse, kind):
    m,g,c=[read(case/(x+'.json')) for x in ('mesh','ground','stages')]
    model=System(m,c['beam_system']);xyz=np.asarray(m['nodes']);base_xyz=np.asarray(coarse['nodes'])
    ancestry=np.arange(len(coarse['tetrahedra']))
    for audit in read(case/'refinement.json'):ancestry=ancestry[np.asarray(audit['parent_of_element'])]
    soil_queries=[(ei,q) for ei in range(len(coarse['tetrahedra'])) for q in ref.volume(coarse,ei)['position']]
    soil_ops=[shape_at(m,q,np.flatnonzero(ancestry==ei)) for ei,q in soil_queries]
    mat={x['id']:x['mohr_coulomb'] for x in g['materials']}
    dmat={k:ref.elastic_tensor(v['young_kpa'],v['poisson']) for k,v in mat.items()}
    # Restrict reactions by coarse virtual displacements on the exact retained
    # support trace, rather than comparing incompatible nodal force inventories.
    supports=m['node_sets']['bottom'];restriction=np.zeros((len(coarse['nodes']),len(supports)))
    for j,n in enumerate(supports):
        _,ids,shape,_=shape_at(coarse,xyz[n]);restriction[ids,j]=shape
    # Fixed physical plate points are the original 2x6 recovery points.
    coarse_case=read(case.parent/'beam-positive/stages.json');coarse_model=System(coarse,coarse_case['beam_system'])
    plate_queries=[np.asarray(l)@coarse_model.positions[nodes][:3] for nodes,_,_ in coarse_model.plate_data for l in SAMPLES]
    plate_ops=[]
    for query in plate_queries:
        for pi,(nodes,section,kernel) in enumerate(model.plate_data):
            points=model.positions[nodes][:3];r=np.linalg.solve((points[1:,:2]-points[0,:2]).T,query[:2]-points[0,:2]);l=np.r_[1-r.sum(),r]
            if l.min()>=-1e-8:
                plate_ops.append((pi,nodes,kernel,l));break
        else:raise ValueError('no plate owner')
    source=case/('native/accepted.jsonl' if kind=='native' else 'independent.jsonl')
    rows={}
    with source.open() as f:
        for line in f:
            row=json.loads(line);si=row['stage_index'];progress=row['progress']
            if progress not in (.25,.5,.75,1.):continue
            s=row['state'] if kind=='native' else row
            u=np.asarray(s['displacement_m']);bs=s['beams'];report=row['beam_results']
            strain=[];stress=[];uq=[]
            for ei,ids,shape,grad in soil_ops:
                du=u[ids].T@grad;eps=np.array([du[0,0],du[1,1],du[2,2],du[0,1]+du[1,0],du[1,2]+du[2,1],du[0,2]+du[2,0]])
                strain.append(eps);stress.append(-dmat[m['material_ids'][ei]]@eps);uq.append(shape@u[ids])
            d=np.c_[bs['node_displacement_m'],bs['node_rotation_rad']]
            fields=dict(soil_common_nodes_m=u[:len(coarse['nodes'])].ravel(),soil_probe_u_m=np.ravel(uq),
                        soil_probe_stress_kpa=np.ravel(stress),soil_probe_strain=np.ravel(strain))
            if kind=='native':
                internal=np.asarray(row['soil_internal_kn']);internal[np.asarray(row['interface_dofs'],int)]+=row['interface_internal_kn']
                residual=(internal-np.asarray(s['applied_nodal_load_kn'])).reshape(-1,3)
                reaction=residual[supports]
                fb=model.assemble(bs,s['boundary_factor'])[2];effective=model.assemble(bs,s['boundary_factor'])[3]
                attached=sum(float(effective[i]@u[n]) for i,n in enumerate(model.attach) if n is not None)
                work=float(np.asarray(s['applied_nodal_load_kn'])@u.ravel()+fb@d.ravel()-attached)
            else:
                reaction=np.asarray(row['support_reaction_by_node_kn'])[supports];work=row['load_work_kn_m']
            fields['support_generalized_reactions_kn']=(restriction@reaction)[coarse['node_sets']['bottom']].ravel()
            fields['support_resultant_kn']=reaction.sum(axis=0);fields['load_work_kn_m']=np.array([work])
            if any(bs['active_plates']):
                values={k:[] for k in ('plate_attachment_u_m','plate_rotation_rad','plate_membrane_global_kn_m','plate_bending_global_kn','plate_attachment_moment_global_kn','plate_shear_global_kn_m')}
                for pi,nodes,kernel,l in plate_ops:
                    kin=d[nodes].ravel()-np.asarray(bs['plate_reference_global_kinematics'][pi])
                    local=kernel.neutral.transform@kernel.map@kin
                    bm,bb,bq,_=kernel.neutral.operators(l)
                    nforce=kernel.neutral.membrane@(bm@local);moment=kernel.neutral.bending@(bb@local);shear=kernel.neutral.shear*(bq@local)
                    if kind == 'native':
                        # Use actual production recovery values. N/M are linear
                        # polynomials; Q is in the existing 8-dimensional
                        # Nedelec space. This is common-point observation only.
                        plate = next(v for v in report.get('plates',[]) if v['id']==model.plates[pi]['id'])
                        xy = np.asarray([np.asarray(v['barycentric'])@kernel.neutral.xy[:3] for v in plate['samples']])
                        fit = np.c_[np.ones(len(xy)),xy]
                        at = np.r_[1.,np.asarray(l)@kernel.neutral.xy[:3]]
                        nforce = at@np.linalg.lstsq(fit,np.asarray([v['membrane_force_kn_m'] for v in plate['samples']]),rcond=None)[0]
                        moment = at@np.linalg.lstsq(fit,np.asarray([v['bending_moment_kn'] for v in plate['samples']]),rcond=None)[0]
                        basis = np.vstack([kernel.neutral.basis(*v) for v in xy])
                        coeff = np.linalg.lstsq(basis,np.asarray([v['transverse_shear_force_kn_m'] for v in plate['samples']]).ravel(),rcond=None)[0]
                        shear = kernel.neutral.basis(*(np.asarray(l)@kernel.neutral.xy[:3]))@coeff
                    shape=kernel.neutral.shapes(*(l@kernel.neutral.xy[:3]))[0]
                    values['plate_attachment_u_m'].extend(shape@d[nodes,:3]);values['plate_rotation_rad'].extend(shape@d[nodes,3:])
                    values['plate_membrane_global_kn_m'].extend(tensor(nforce,kernel.axes));values['plate_bending_global_kn'].extend(tensor(moment,kernel.axes))
                    values['plate_attachment_moment_global_kn'].extend(tensor(moment-kernel.e*nforce,kernel.axes))
                    values['plate_shear_global_kn_m'].extend(kernel.axes.T@np.r_[shear,0.])
                fields.update({k:np.asarray(v) for k,v in values.items()})
            for key in ('membrane_energy_kn_m','bending_energy_kn_m','transverse_shear_energy_kn_m','drilling_penalty_energy_kn_m'):
                plates=report.get('plates',[]);plates=plates.values() if isinstance(plates,dict) else plates
                fields['sum_'+key]=np.array([sum(p[key] for p in plates)])
            tip=model.node_ids['tip'];fields['beam_tip_kinematics']=d[tip]
            members=report.get('members',[]);members=members.values() if isinstance(members,dict) else members
            fields['beam_section_forces']=np.ravel([p['section_forces_at_start_mid_end'] for p in members])
            rows[(si,progress)]=fields
    return rows,dict(soil_probe_count=len(soil_ops),plate_probe_count=len(plate_ops),common_soil_nodes=len(coarse['nodes']),
                     reaction_restriction='coarse quadratic virtual displacement transpose, exact support trace',
                     stress_observation='elastic strain from actual native/independent quadratic displacement, not stored stress interpolation',
                     plate_observation='native actual six-point N/M/Q recovered as existing polynomial fields; independent fields from its independently solved kinematics')


def differences(a,b):
    rows={};aggregate={}
    assert a.keys()==b.keys()
    for point in a:
        assert a[point].keys()==b[point].keys()
        fields={}
        for key in a[point]:
            delta=np.asarray(a[point][key])-b[point][key]
            fields[key]=dict(absolute_l2=float(np.linalg.norm(delta)),maximum_absolute=float(np.max(np.abs(delta),initial=0.)))
            aggregate[key]=aggregate.get(key,0.)+float(np.sum(delta**2))
        rows[str(point)]=fields
    return dict(all_common_points=rows,absolute_l2_all_points={k:float(np.sqrt(v)) for k,v in aggregate.items()})


def run(formal,output):
    coarse=read(formal/'beam-positive/mesh.json');result={}
    for series,names in [('space',['beam-positive','mesh-1','mesh-2']),('load',['beam-positive','load-8','load-16'])]:
        result[series]={}
        for kind in ('native','independent'):
            values=[];meta=None
            for name in names:
                r,meta=sample(formal/name,coarse,kind);values.append(r)
            pairs=[differences(values[i],values[i+1]) for i in range(2)]
            orders={}
            for field,first in pairs[0]['absolute_l2_all_points'].items():
                second=pairs[1]['absolute_l2_all_points'][field]
                orders[field]=dict(coarse_to_medium=first,medium_to_fine=second,
                                   observed_order=float(np.log2(first/second)) if first>0 and second>0 else None,
                                   decreasing=second<first,
                                   caution='roundoff and zero differences are not an asymptotic order certificate' if series=='load' else 'finite h-series, no new convergence gate')
            result[series][kind]=dict(names=names,common_points=len(values[0]),metadata=meta,pairs=pairs,trends=orders)
    write(output,dict(passed=False,b02_complete=False,series=result,
                      note='Space uses common physical probes and work-conjugate support restriction. Load series is fixed mesh at common accepted progress; purely elastic paths are roundoff dominated. No extrapolation or new tolerance.'))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('formal',type=Path);p.add_argument('output',type=Path);a=p.parse_args();run(a.formal.resolve(),a.output.resolve())
