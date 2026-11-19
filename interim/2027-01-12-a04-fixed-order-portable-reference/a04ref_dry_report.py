"""Fixed-state replay, physical point correspondence and exact dyadic reduction."""
import argparse
from collections import Counter
from fractions import Fraction
import json
import math
from pathlib import Path
import sys
from types import SimpleNamespace

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'2026-12-31-a04-richards-wet-joint-cloud-followup'))
import a04follow_joint_order as joint
from a04follow_common import np,read,write
from verify_joint_coupling import contact


def peak(a,b):
    return float(np.max(np.abs(np.asarray(a)-np.asarray(b)),initial=0.))


def l2(values):
    return math.sqrt(math.fsum(float(x)**2 for x in np.ravel(values,order='C')))


def fixed(rows,nd):
    bins=[[] for _ in range(nd)]
    for key,ids,force in sorted(rows,key=lambda r:r[0]):
        for i,v in zip(ids,force,strict=True):
            bins[i].append(float(v))
    # Explicit loop, including on Python versions whose builtin sum is compensated.
    def ordered(v):
        s=0.
        for x in v:
            s+=x
        return s
    def exact(v):
        return float(sum(map(Fraction.from_float,v),Fraction()))
    return {name:np.array([func(v) for v in bins]) for name,func in
            [('fixed',ordered),('fsum',math.fsum),('exact_dyadic',exact)]}


def face(p):
    return (tuple(sorted(p['nodes_a'])),tuple(sorted(p['nodes_b'])),p['owner_a'],p['owner_b'])


def mapping(native,reference):
    # Match only within the exact side/owner pair using all six nodal shape values.
    # Require mutual and unique assignments, record the actual geometric differences.
    def distances(a,b):
        if face(a)!=face(b):
            return math.inf
        aa=dict(zip(a['nodes_a'],a['shape']))
        return max(abs(aa[n]-w) for n,w in zip(b['nodes_a'],b['shape']))
    d=np.array([[distances(a,b) for a in native] for b in reference])
    order=d.argmin(axis=1)
    assert len(native)==len(reference)==len(set(order))
    assert np.array_equal(d.argmin(axis=0)[order],np.arange(len(order)))
    nearest=np.sort(d,axis=1)
    assert np.all(np.isfinite(nearest[:,0]) & (nearest[:,0]<nearest[:,1]))
    return order,dict(paired=len(order),max_shape_difference=float(nearest[:,0].max()),
                      min_runner_up_distance=float(nearest[:,1].min()),
                      max_position_difference=max(peak(native[i]['position_m'],p['position_m'])
                          for i,p in zip(order,reference)))


def replay(model,u,old,counts):
    s=model.system
    jump=(s.c@u).reshape(-1,3)
    rest=np.array([v['birth']['elastic_offset_m'] for v in old])
    before=np.array([v['jump_m'] for v in old])-rest
    after=jump-rest
    prior=np.array([v['traction_kpa'] for v in old])
    traction=np.zeros_like(jump)
    status=np.empty(len(jump),dtype=object)
    for ct in np.unique(counts):
        mask=counts==ct
        sub=SimpleNamespace(**{k:getattr(s,k)[mask] for k in ('normal','kn','kt','mu')})
        a,b=before[mask],after[mask]
        q=prior[mask]
        traction[mask]=contact(sub,a,q,b,int(ct))[0][:,0]
        prev=a+(b-a)*(int(ct)-1)/int(ct)
        prevq=contact(sub,a,q,prev,int(ct)-1)[0][:,0] if ct>1 else q
        flow=contact(sub,prev,prevq,b,1)[1][:,0]
        status[mask]=np.where((b*sub.normal).sum(axis=1)>0,'open',
                              np.where(np.any(flow!=0,axis=1),'slip','stick'))
    return traction,status


def observations(model,raw,cells):
    bulk=[];contacts=[];info=[]
    for term in raw['rows']:
        if term[0]=='bulk':
            _,e,_,xyz,ids,v=term
            bulk.append(((0,cells[e],tuple(xyz)),ids,v))
        else:
            _,i,p,coeff,v=term
            assert len(info)==i
            info.append(p)
            key=(1,tuple(sorted(p['nodes_a'])),tuple(p['position_m']))
            ids=[3*n+a for n,_ in coeff for a in range(3)]
            contacts.append((key,ids,np.ravel(v)))
    order,match=mapping(info,model.system.points)
    u=np.asarray(raw['u']).ravel()
    old=[raw['old_joint'][i] for i in order]
    counts=np.array(raw['partitions'])[order]
    tr,status=replay(model,u,old,counts)
    native_status=[raw['joint'][i]['status'] for i in order]
    ref_terms=[]
    for e,(op,ids,d,s0,*_) in enumerate(model.ops):
        strains=op['b']@u[ids]
        stress=s0-strains@d.T
        force=-(op['b'].transpose(0,2,1)@stress[:,:,None])[:,:,0]*op['weight'][:,None]
        for xyz,v in zip(op['position'],force):
            ref_terms.append(((0,cells[e],tuple(xyz)),ids,v))
    contact_terms=[]
    for i,p in enumerate(model.system.points):
        c=model.system.c[i]
        ids=np.flatnonzero(np.any(c!=0,axis=0))
        contact_terms.append(((1,tuple(sorted(p['nodes_a'])),tuple(p['position_m'])),ids,
                          (tr[i]@c[:,ids])*model.system.weight[i]))
    native=fixed(bulk+contacts,model.nd)
    independent=fixed(ref_terms+contact_terms,model.nd)
    fifth_rule=None
    if raw['iteration']==5:
        # The self-solved reference stiffness uses the independent order-five rule.
        import fem_reference as fem
        terms=list(contact_terms)
        for e,(_,ids,d,s0,*_) in enumerate(model.ops):
            op=fem.volume(model.mesh,e,independent=True,independent_order=5)
            stress=s0-(op['b']@u[ids])@d.T
            force=-(op['b'].transpose(0,2,1)@stress[:,:,None])[:,:,0]*op['weight'][:,None]
            terms.extend(((0,cells[e],tuple(xyz)),ids,v) for xyz,v in zip(op['position'],force))
        fifth_rule=fixed(terms,model.nd)
    true_native=np.asarray(raw['internal'])
    load=np.asarray(raw['load'])
    t=model.spaces[False][0]
    def reduced(force):
        net=force-load
        projected=[math.fsum(float(col[i])*float(net[i]) for i in np.flatnonzero(col)) for col in t.T]
        return l2(projected)/raw['scale']
    return dict(iteration=raw['iteration'],relative=raw['relative'],linear=raw['linear'],
        native=native,independent_same_state=independent,independent_GL5_same_state=fifth_rule,point_mapping=match,
        states_native=native_status,states_reference=status.tolist(),state_counts=dict(Counter(native_status)),
        state_mismatches=int(np.count_nonzero(status!=native_status)),
        max_traction_difference_kpa=peak(tr,[raw['joint'][i]['traction_kpa'] for i in order]),
        original_to_fixed_kn=peak(true_native,native['fixed']),
        original_to_exact_kn=peak(true_native,native['exact_dyadic']),
        same_state_native_to_reference_kn=peak(native['exact_dyadic'],independent['exact_dyadic']),
        reduced_on_native={k:reduced(f) for k,f in native.items()},
        reduced_on_reference={k:reduced(f) for k,f in independent.items()},
        u=u,load=load,reaction=native['exact_dyadic']-load)


def report(root,out):
    _,ref,_=joint.binding('dry')
    records={}
    for folder in sorted(root.glob('e?f?q?')):
        mesh,ground,case=[read(folder/(k+'.json')) for k in ('mesh','ground','case')]
        cells=[tuple(sorted(el)) for el in mesh['tetrahedra10']]
        with joint.controls(folder.name.endswith('q1')):
            model=ref.Model(mesh,ground,case)
            assert not(model.frame or model.tendon)
            with (folder/'iterates.jsonl').open() as stream:
                rows=[observations(model,json.loads(line),cells) for line in stream]
        own=next(v for v in map(json.loads,(folder/'independent-fields.jsonl').read_text().splitlines())
                 if v['row']['stage']=='open' and v['row']['progress']==.25)
        records[folder.name]=dict(iterates=rows,independent_endpoint_reaction=own['fields']['reaction'],
            own_endpoint_difference_kn=peak(rows[-1]['reaction'],np.ravel(own['fields']['reaction'])))
    assert len(records)==8
    baseline=records['e0f0q0']['iterates']
    for path in records.values():
        for row in path['iterates']:
            b=baseline[row['iteration']]
            row['vs_base_iterate']=dict(
               u_m=peak(row['u'],b['u']), load_kn=peak(row['load'],b['load']),
               reaction_kn=peak(row['reaction'],b['reaction']))
        path['vs_base_endpoint_kn']=peak(path['iterates'][-1]['reaction'],baseline[-1]['reaction'])
    write(out,records)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--raw',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    report(a.raw,a.output)
