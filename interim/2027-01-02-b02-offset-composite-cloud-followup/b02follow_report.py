"""Physical-ID fixed-order observations; no tolerance or acceptance replacement."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from b02follow_common import HERE, legacy, quadrature, write, sha, dot, norm, metric, common
import b02cause_report as prior


def fixed_soil(mesh, mapping, row):
    """Reassemble actual saved stresses at unchanged production GP, by physical ID."""
    terms = [[] for _ in range(3*len(mesh['nodes']))]
    state = row.get('state',row)
    stress = ([s['effective_stress_compression_positive_kpa'] for s in state['material_states']]
              if 'state' in row else row['stress_kpa'])
    offset=0
    blocks=[]
    for ei,original in enumerate(mapping['tetrahedra_declared_to_original']):
        qp=legacy['cases'].ref.volume(mesh,ei)
        blocks.append((original,qp,np.asarray(stress[offset:offset+len(qp['weight'])])))
        offset += len(qp['weight'])
    if offset != len(stress):
        raise ValueError('stress GP inventory mismatch')
    for original,qp,sig in sorted(blocks,key=lambda v:v[0]):
        for q in range(len(qp['weight'])):
            for local,node in enumerate(qp['nodes']):
                for a in range(3):
                    terms[3*node+a].extend(-float(qp['weight'][q])*float(qp['b'][q,j,3*local+a])*float(sig[q,j]) for j in range(6))
    return np.array([math.fsum(v) for v in terms])


def cumulative_work(folder, kind):
    m,c=[legacy['cases'].read(folder/(s+'.json')) for s in ('mesh','stages')]
    model=prior.System(m,c['beam_system']);nd=3*len(m['nodes']);nb=model.size
    old_u=np.zeros(nd);old_d=np.zeros(nb);old_soil=np.zeros(nd);old_load=np.zeros(nb)
    old_active=np.zeros(len(model.nodes),bool);increments=[[],[]];out=[]
    for row in prior.rows(folder,kind):
        s=row.get('state',row);bs=s['beams'];u=np.asarray(s['displacement_m']).ravel()
        d=np.c_[bs['node_displacement_m'],bs['node_rotation_rad']].ravel()
        _,_,load,eff,*_=model.assemble(bs,prior.factor(c,row))
        soil=np.asarray(s['applied_nodal_load_kn']).copy()
        for i,node in enumerate(model.attach):
            if node is not None:soil[3*node:3*node+3]-=eff[i]
        active=model.active_nodes(bs);previous=old_d.copy();refs={}
        # Same old endpoint definition. Sort stable IDs; incompatible refs are errors.
        for i in sorted(range(len(model.plates)),key=lambda i:model.plates[i]['id']):
            if bs['active_plates'][i]:
                for j,node in enumerate(model.plate_data[i][0]):refs[node]=bs['plate_reference_global_kinematics'][i][6*j:6*j+6]
        for i in sorted(range(len(model.members)),key=lambda i:model.members[i]['id']):
            if bs['active_members'][i]:
                for j,node in enumerate(model.data[i][0]):refs[node]=bs['reference_global_kinematics'][i][6*j:6*j+6]
        for node in np.flatnonzero(active & ~old_active):previous[6*node:6*node+6]=refs[node]
        inc=[dot(.5*(soil+old_soil),u-old_u),dot(.5*(load+old_load),d-previous)]
        for target,value in zip(increments,inc,strict=True):target.append(value)
        out.append(dict(stage_index=row['stage_index'],progress=row['progress'],
                        soil_increment_kn_m=inc[0],structure_increment_kn_m=inc[1],
                        cumulative_soil_kn_m=math.fsum(increments[0]),cumulative_structure_kn_m=math.fsum(increments[1]),
                        cumulative_total_kn_m=math.fsum(math.fsum(v) for v in increments),
                        instantaneous_kn_m=math.fsum([dot(soil,u),dot(load,d)])))
        old_u=u.copy();old_d=d.copy();old_soil=soil.copy();old_load=load.copy();old_active=active.copy()
    return out


def matrices(folder, mapping, order):
    """Same local matrices: declared += versus physical-ID fsum; no new solve."""
    mesh,ground,case=[legacy['cases'].read(folder/(s+'.json')) for s in ('mesh','ground','stages')]
    nd=3*len(mesh['nodes']);terms=[];seq=np.zeros((nd,nd))
    mats={s['id']:s['mohr_coulomb'] for s in ground['materials']}
    for ei,original in enumerate(mapping['tetrahedra_declared_to_original']):
        q=legacy['cases'].ref.volume(mesh,ei,independent=True)
        mat=mats[mesh['material_ids'][ei]];tensor=legacy['cases'].ref.elastic_tensor(mat['young_kpa'],mat['poisson'])
        idx=(q['nodes'][:,None]*3+np.arange(3)).ravel()
        local=np.einsum('qai,ab,qbj,q->ij',q['b'],tensor,q['b'],q['weight'],optimize=True)
        seq[np.ix_(idx,idx)]+=local;terms.append((original,idx,local))
    def accumulate(size,values):
        entries={}
        for _,idx,local in sorted(values,key=lambda v:v[0]):
            for a,i in enumerate(idx):
                for b,j in enumerate(idx):entries.setdefault((int(i),int(j)),[]).append(float(local[a,b]))
        fixed=np.zeros((size,size))
        for (i,j),v in sorted(entries.items()):fixed[i,j]=math.fsum(v)
        return fixed
    fixed=accumulate(nd,terms)
    model=prior.System(mesh,case['beam_system']);s=model.initial_state()
    s['active_plates']=[True]*len(model.plates);s['active_members']=[True]*len(model.members)
    full=model.assemble(s,0.)[0]
    # Existing beam stiffness; isolate all plate terms without changing a kernel.
    pterms=[];plate_sum=np.zeros_like(full)
    for i,(nodes,_,kernel) in enumerate(model.plate_data):
        idx=np.array([6*n+a for n in nodes for a in range(6)])
        plate_sum[np.ix_(idx,idx)]+=kernel.stiffness
        pterms.append((model.plates[i]['id'],idx,kernel.stiffness))
    s['active_plates']=[False]*len(model.plates)
    beam=model.assemble(s,0.)[0]
    fixed_plate=accumulate(model.size,pterms)
    fixed_structure=np.array([[math.fsum([float(beam[i,j]),float(fixed_plate[i,j])]) for j in range(model.size)] for i in range(model.size)])
    path=folder/f'assembly-gl{order}.npz'
    np.savez_compressed(path,soil_sequential=seq,soil_fixed=fixed,structure_sequential=full,
                        structure_fixed=fixed_structure,plate_sequential=plate_sum,plate_fixed=fixed_plate)
    digest=lambda a:hashlib.sha256(np.asarray(a,dtype='<f8').tobytes(order='C')).hexdigest()
    return dict(soil_fixed_sha256=digest(fixed),structure_fixed_sha256=digest(fixed_structure),
                soil_sequential_to_fixed=metric(seq,fixed),structure_sequential_to_fixed=metric(full,fixed_structure),
                raw_matrix_file=path.name,raw_matrix_sha256=sha(path))


def observe(folder, coarse, kind, mapping):
    mesh,case=[legacy['cases'].read(folder/(s+'.json')) for s in ('mesh','stages')]
    rawrows=prior.rows(folder,kind)
    details=prior.force_work_energy(folder,kind)
    probes=prior.plate_probes(folder,coarse,kind)
    work=cumulative_work(folder,kind)
    fixed={3*n+a for n in mesh['node_sets']['bottom'] for a in range(3)}
    free=[i for i in range(3*len(mesh['nodes'])) if i not in fixed]
    records=[]
    model=prior.System(mesh,case['beam_system'])
    for row,detail,w in zip(rawrows,details,work,strict=True):
        s=row.get('state',row);soil=fixed_soil(mesh,mapping,row)
        mapped=np.asarray(detail['force_and_moment_vectors']['mapped_structure_internal'])
        applied=np.asarray(s['applied_nodal_load_kn'])
        residual=np.array([math.fsum([float(a),float(b),-float(c)]) for a,b,c in zip(soil,mapped,applied,strict=True)])
        reaction=np.zeros(len(soil));reaction[sorted(fixed)]=residual[sorted(fixed)]
        detail.update(soil_fixed_order_internal_kn=soil.tolist(),fixed_free_residual_kn=residual[free].tolist(),
                      fixed_free_residual_l2_kn=norm(residual[free]),
                      fixed_free_residual_relative_algebraic=norm(residual[free])/max(norm(applied),1.),
                      support_reaction_fixed_order_kn=reaction.tolist())
        if (row['stage_index'],row['progress']) not in probes:
            continue
        vectors=detail['force_and_moment_vectors']
        struct=np.asarray(vectors['plate_structure'])+np.asarray(vectors['beam_structure'])-np.asarray(vectors['structure_load'])
        basis=model.space(s['beams'],np.asarray(s['displacement_m']))[2]
        struct_free=(basis@(basis.T@struct)).reshape(-1,6)
        records.append(dict(stage_index=row['stage_index'],progress=row['progress'],
                            translation_m=s['displacement_m'],rotation_rad=s['beams']['node_rotation_rad'],
                            reaction_kn=reaction.tolist(),force_residual_kn=residual[free].tolist(),
                            structure_residual_force_kn=struct_free[:,:3].tolist(),structure_residual_moment_kn_m=struct_free[:,3:].tolist(),q_kn_m=probes[(row['stage_index'],row['progress'])]['q'],
                            instantaneous_work_kn_m=w['instantaneous_kn_m'],cumulative_work_kn_m=w['cumulative_total_kn_m']))
    write(folder/(kind+'-observations.json'),dict(records=records,details=details,work=work))
    return records,details


def compare(a,b):
    if [(r['stage_index'],r['progress']) for r in a] != [(r['stage_index'],r['progress']) for r in b]:
        raise ValueError('path mismatch')
    result={}
    for key in ('translation_m','rotation_rad','reaction_kn','force_residual_kn','structure_residual_force_kn','structure_residual_moment_kn_m','q_kn_m','instantaneous_work_kn_m','cumulative_work_kn_m'):
        av=[];bv=[]
        for ar,br in zip(a,b,strict=True):
            aa=np.asarray(ar[key]).ravel();bb=np.asarray(br[key]).ravel()
            if aa.shape!=bb.shape:raise ValueError('activity/component mismatch')
            av.extend(aa);bv.extend(bb)
        result[key]=metric(av,bv)
    return result


def run(formal, output):
    plan=legacy['cases'].read(formal/'plan.json');data={};cases={};assembly={};maxima={}
    for rec in plan['cases']:
        name=rec['name'];folder=formal/name;coarse=formal/f"tet12-h{rec['scale']:g}-default"
        data[name]={};assembly[name]={}
        for kind,order in [('native',5),('gl4',4),('gl5',5),('gl6',6)]:
            with quadrature(order):
                values,details=observe(folder,coarse,kind,rec['mapping'])
                data[name][kind]=values
                if kind!='native':assembly[name][kind]=matrices(folder,rec['mapping'],order)
            for d in details:
                for key in ('fixed_free_residual_l2_kn','fixed_free_residual_relative_algebraic','explicit_free_structure_force_l2_kn','explicit_free_structure_moment_l2_kn_m'):
                    maxima[key]=max(maxima.get(key,0.),d[key])
        cases[name]=dict(native_vs_gl5=compare(data[name]['native'],data[name]['gl5']),
                         quadrature={k:compare(data[name][k],data[name]['gl5']) for k in ('gl4','gl6')},
                         legacy_existing_gates=legacy['cases'].read(folder/'existing-gate-result.json')['existing_field_gates_passed'])
        print(name,'observed',flush=True)
    order_effects={};thickness={};spatial={};factorial=[]
    for rec in plan['cases']:
        name=rec['name'];base=f"tet{12*8**rec['level']}-h{rec['scale']:g}-default"
        order_effects[name]={k:compare(data[name][k],data[base][k]) for k in data[name]}
        thickness[name]={k:compare(data[name][k],data[f"tet{12*8**rec['level']}-h1-{rec['mode']}"][k]) for k in data[name]}
        if rec['level']==1:
            coarse=f"tet12-h{rec['scale']:g}-{rec['mode']}"
            spatial[name]={}
            for kind in data[name]:
                # Full meshes have different DOF inventories; only shared physical probes/work.
                keys=('q_kn_m','instantaneous_work_kn_m','cumulative_work_kn_m')
                spatial[name][kind]={}
                for key in keys:
                    av=np.concatenate([np.asarray(v[key]).ravel() for v in data[coarse][kind]])
                    bv=np.concatenate([np.asarray(v[key]).ravel() for v in data[name][kind]])
                    spatial[name][kind][key]=metric(av,bv)
        for kind in ('gl4','gl5','gl6'):
            same={key:assembly[name][kind][key]==assembly[base][kind][key] for key in ('soil_fixed_sha256','structure_fixed_sha256')}
            if not all(same.values()):raise ValueError('physical-ID fixed assembly differs across declaration orders')
            factorial.append(dict(case=name,rule=kind,canonical_input_same=rec['mapping']['canonical_input_sha256']==next(v['mapping']['canonical_input_sha256'] for v in plan['cases'] if v['name']==base),fixed_matrix_exact=same))
    # Full signed inclusion-exclusion contrasts, at identical physical path coordinates.
    effects=[]
    keys=('translation_m','rotation_rad','reaction_kn','force_residual_kn','structure_residual_force_kn','structure_residual_moment_kn_m','q_kn_m','instantaneous_work_kn_m','cumulative_work_kn_m')
    for rec in plan['cases']:
        level=12*8**rec['level'];h=rec['scale'];mode=rec['mode']
        def response(hh,oo,qq,key):
            return np.concatenate([np.asarray(v[key]).ravel() for v in data[f'tet{level}-h{hh:g}-{oo}'][f'gl{qq}']])
        for q in (4,5,6):
            raw={}
            for key in keys:
                y=lambda hh,oo,qq:response(hh,oo,qq,key)
                base=y(1.,'default',5)
                def signed(items):
                    return [math.fsum(float(a[i])*sign for a,sign in items) for i in range(len(base))]
                raw[key]=dict(
                    thickness=signed([(y(h,'default',5),1),(base,-1)]),
                    declaration_order=signed([(y(1.,mode,5),1),(base,-1)]),
                    quadrature=signed([(y(1.,'default',q),1),(base,-1)]),
                    thickness_order=signed([(y(h,mode,5),1),(y(h,'default',5),-1),(y(1.,mode,5),-1),(base,1)]),
                    thickness_quadrature=signed([(y(h,'default',q),1),(y(h,'default',5),-1),(y(1.,'default',q),-1),(base,1)]),
                    order_quadrature=signed([(y(1.,mode,q),1),(y(1.,mode,5),-1),(y(1.,'default',q),-1),(base,1)]),
                    triple=signed([(y(h,mode,q),1),(y(h,mode,5),-1),(y(h,'default',q),-1),(y(1.,mode,q),-1),
                                   (y(h,'default',5),1),(y(1.,mode,5),1),(y(1.,'default',q),1),(base,-1)]))
            effects.append(dict(case=rec['name'],rule=q,signed_effects=raw))
    write(formal/'signed-factorial-effects.json',effects)
    result=dict(schema='b02follow-measurements/1',passed=False,b02_complete=False,diagnostic_only=True,
                numerical_acceptance_all_fields=None,single_cause_identified=False,
                native_histories=18,independent_histories=54,accepted_states=sum(len(prior.rows(formal/r['name'],k)) for r in plan['cases'] for k in ('native','gl4','gl5','gl6')),
                cases=cases,order_effects=order_effects,thickness_effects=thickness,spatial_12_to_96=spatial,
                fixed_order_matrix_controls=factorial,assembly=assembly,maximum=maxima,
                unsupported_combined_equilibrium=None,residual_768='quotation only: 5e3d191; no new calculation',
                order='physical element ID, production GP, component; physical plate ID; node/DOF; math.fsum',
                mixed_dof_observations='direct_linear_residual in raw independent rows and all structure matrix norms are algebraic observations; force/moment units are separated above',
                source_sha256=sha(Path(__file__)))
    write(output,result)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.run.resolve(),a.output.resolve())
