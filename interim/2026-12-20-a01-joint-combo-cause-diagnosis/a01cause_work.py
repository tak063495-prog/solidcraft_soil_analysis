"""Read-only work accounting on saved paths; never a global equilibrium solve."""
from types import SimpleNamespace
import math
import numpy as np
from a01cause_binding import original_reference as old
from a01cause_metrics import norm,metric,budget


def identity(rows,a,b):
    if a['jump'].shape!=b['jump'].shape:raise ValueError('GP correspondence differs')
    n,q,_=a['jump'].shape
    traction=np.zeros((n,q));jump=np.zeros_like(traction)
    rules={side:{rule:np.zeros((n,q)) for rule in ('left','trapezoid','right')} for side in ('native','reference')}
    for i in range(n):
        live=bool(np.all(a['joint_active'][i]));assert live==bool(np.all(b['joint_active'][i]))
        reset=i==0 or not np.array_equal(a['generation'][i],a['generation'][i-1])
        if i and not reset:
            traction[i]=traction[i-1];jump[i]=jump[i-1]
            for values in rules.values():
                for v in values.values():v[i]=v[i-1]
        if not live:continue
        changes=[]
        for side,f in (('native',a),('reference',b)):
            j0=f['reference_jump'][i] if reset else f['jump'][i-1]
            t0=f['initial_traction'][i] if reset else f['traction'][i-1]
            dj=f['jump'][i]-j0;changes.append(dj)
            for rule,t in (('left',t0),('trapezoid',.5*(t0+f['traction'][i])),('right',f['traction'][i])):
                rules[side][rule][i]+=np.sum(t*dj,axis=1)
        da,db=changes;ta,tb=a['traction'][i],b['traction'][i]
        traction[i]+=np.sum(.5*(ta-tb)*(da+db),axis=1)
        jump[i]+=np.sum(.5*(ta+tb)*(da-db),axis=1)
    actual=a['endpoint_work']-b['endpoint_work']
    endpoint=rules['native']['right']-rules['reference']['right']
    remainder=actual-traction-jump
    assert np.max(np.abs(endpoint-traction-jump),initial=0)<1e-12
    return dict(traction_term=traction,jump_term=jump,within_increment_remainder=remainder,
                identity_residual_max=float(np.max(np.abs(actual-traction-jump-remainder),initial=0)),
                budgets={key:budget(actual,v) for key,v in [('traction',traction),('jump',jump),('within_increment_remainder',remainder)]},
                rules={side:{rule:metric(v,a['endpoint_work'] if side=='native' else b['endpoint_work']) for rule,v in values.items()} for side,values in rules.items()},
                actual_difference=metric(a['endpoint_work'],b['endpoint_work']),
                stage_endpoints=metric(a['endpoint_work'][[i for i,r in enumerate(rows) if r['progress']==1]],b['endpoint_work'][[i for i,r in enumerate(rows) if r['progress']==1]]))


def replay(f,s,count):
    """Common contact subdivision on a saved displacement path. No equilibrium."""
    n,q,_=f['jump'].shape;work=np.zeros((n,q));traction=np.zeros((q,3));jump=np.zeros_like(traction)
    total=np.zeros(q);rest=np.zeros_like(jump)
    for i in range(n):
        reset=i==0 or not np.array_equal(f['generation'][i],f['generation'][i-1])
        if reset:
            traction=f['initial_traction'][i].copy();jump=f['reference_jump'][i].copy()
            rest=f['elastic_offset'][i].copy();total=np.zeros(q)
        if np.all(f['joint_active'][i]):
            target=f['jump'][i]
            total+=old.contact_work(s,jump-rest,traction,target-rest,count)[2]
            traction=old.contact(s,jump-rest,traction,target-rest,count)[0][:,0]
            jump=target.copy()
        work[i]=total
    return work


def geometry(s):
    ids=[];faces={};area={};shared={}
    for p in s.points:
        face=(p['joint_id'],p['face_index']);ids.append((face,tuple(p['position_m'])))
        faces[face]=(tuple(sorted(p['nodes_a'])),tuple(sorted(p['nodes_b'])),p['owner_a'],p['owner_b'])
        area[p['joint_id']]=area.get(p['joint_id'],0.)+p['weight_m2']
    if len(ids)!=len(set(ids)):raise ValueError('duplicate GP')
    physical=[(x[0],x[1]) for x in faces.values()]
    if len(physical)!=len(set(physical)):raise ValueError('double-counted face')
    names=sorted(area)
    for i,x in enumerate(names):
        a=set(v for p in s.points if p['joint_id']==x for v in p['nodes_a']+p['nodes_b'])
        for y in names[i+1:]:
            b=set(v for p in s.points if p['joint_id']==y for v in p['nodes_a']+p['nodes_b'])
            shared[x+' / '+y]=len(a&b)
    return dict(gauss_points=len(ids),unique_gauss_points=len(set(ids)),unique_paired_faces=len(faces),
                duplicate_faces=0,area_m2_by_joint={k:math.fsum(p['weight_m2'] for p in s.points if p['joint_id']==k) for k in names},
                shared_nodes=shared,matching='joint_id/face_index/fixed rule ordinal; independent System checks node/owner/shape/normal/position/weight against native layout; no nearest-point remapping')


def analyze(folder):
    from a01cause_binding import fixtures,fields
    from a01cause_metrics import load
    m=old.Model(*(fixtures.read(folder/(k+'.json')) for k in ('mesh','ground','case')))
    rows,a=load(folder/'native-fields.jsonl');rr,b=load(folder/'independent-fields.jsonl')
    assert [(r['stage'],r['progress']) for r in rows]==[(r['stage'],r['progress']) for r in rr]
    layout=fixtures.read(folder/'native/run.json')['joint_integration_points']
    if layout is not None:
        from joint_reference import System
        System(m.mesh,m.case['joint_system'],layout)
    result=identity(rows,a,b)
    # Retain signed per-GP arrays, not just a ratio or a sum of nonadditive norms.
    result['common_subdivision_saved_path_replays']={str(c):metric(replay(a,m.system,c),replay(b,m.system,c)) for c in (1,2,4,8,16)}
    actual=a['endpoint_work']-b['endpoint_work']
    result['weighted_total_work_difference_kn_m']=[math.fsum(float(w)*float(v) for w,v in zip(m.system.weight,row)) for row in actual]
    result['geometry']=geometry(m.system)
    result['local_replay_is_global_equilibrium']=False
    return result
