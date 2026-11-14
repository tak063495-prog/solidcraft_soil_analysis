"""Explicit independent diagnostic interventions; no production/default change."""
import argparse
from pathlib import Path
import time
import numpy as np
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import mapping as m
c=m.c

def contact_rule(order):
    # Isolated class copy: only the declared surface quadrature order changes.
    import joint_reference,joint_water_reference,inspect
    source=inspect.getsource(joint_reference.System)
    assert source.count('ref.jacobi(4,1)')==1 and source.count('ref.jacobi(4,0)')==1
    source=source.replace('ref.jacobi(4,1)',f'ref.jacobi({order},1)').replace('ref.jacobi(4,0)',f'ref.jacobi({order},0)')
    namespace=dict(vars(joint_reference));exec(compile(source,'<explicit-independent-contact-quadrature>','exec'),namespace)
    joint_water_reference.System=namespace['System']
    return source

def solve(out,order,g):
    source=contact_rule(order);folder=out/f'contact-q{order}-g{g}';folder.mkdir()
    (folder/'independent-system-used.py').write_text(source)
    started=time.monotonic();writer=c.Writer(folder/'independent')
    mesh=c.read(out/'mesh-h0.json');ground=c.read(out/'ground.json');case=c.read(out/'stages.json')
    value=c.old.probe._reference(mesh,ground,case,g,1,7,False,True,writer,maximum_m=5e-5,maximum_dofs=3000)
    writer.close();c.write(folder/'result.json',value)
    c.write(folder/'audit.json',c.old.audit(folder/'independent',mesh,ground,case,'cut',g))
    system=c.old.Water(mesh,ground,case).contact
    c.write(folder/'layout.json',system.points)
    c.write(folder/'timing.json',dict(elapsed_seconds=time.monotonic()-started,contact_order=order,input_only=True,production_or_default_changed=False))

def lagrange(x,nodes):
    result=np.ones((len(x),len(nodes)))
    for j in range(len(nodes)):
        for k in range(len(nodes)):
            if j!=k:result[:,j]*=(x-nodes[k])/(nodes[j]-nodes[k])
    return result

def interpolate(order,points,target):
    u=target[:,0];v=target[:,1]/(1-u)
    return (lagrange(u,c.fem.jacobi(order,1)[0])[:,:,None]*lagrange(v,c.fem.jacobi(order,0)[0])[:,None,:]).reshape(len(target),order**2)

def water(out,target,maps):
    ground=c.read(out/'ground.json');case=c.read(out/'stages.json');level=ground['initialization']['water_level_m'];gamma=ground['initialization']['water_unit_weight_kn_m3']
    values=[]
    for h,mesh in enumerate(maps.meshes):
        cut=np.zeros((len(mesh['nodes']),3));gj=cut.copy()
        for e in range(len(mesh['tetrahedra'])):
            q=c.fem.volume(mesh,e)
            cut[q['nodes']]+=c.old.cut_force(mesh,e,level,gamma)['force']
            pore=gamma*np.maximum(level-q['position'][:,2],0)
            gj[q['nodes']]+=-np.einsum('qai,q,q->i',q['b'][:,:3],pore,q['weight']).reshape(-1,3)
        interface=c.old.Water(mesh,ground,case).hydrostatic(np.ones(len(mesh['tetrahedra']),bool),level)['load']
        values.append(dict(cut=cut,gj=gj,interface=interface))
    results={}
    fine=values[-1]
    for h,v in enumerate(values):
        p=maps.operators[h]['u'];results[f'h{h}']={}
        for key in ('cut','gj','interface'):
            restricted=p.T@fine[key];d=v[key]-restricted
            results[f'h{h}'][key]=dict(coarse_vs_restricted_finest_l2_kn=float(np.linalg.norm(d)),maximum_absolute_kn=float(abs(d).max()),resultant_kn=v[key].sum(axis=0),moment_kn_m=np.cross(np.asarray(maps.meshes[h]['nodes']),v[key]).sum(axis=0))
        results[f'h{h}']['same_mesh_gj_minus_cut_l2_kn']=float(np.linalg.norm(v['gj']-v['cut']))
    c.write(target/'water-affine-cut-diagnostic.json',dict(physical_pressure_time_independent=True,force_comparison_uses_dual_FE_restriction=True,fields=results,exact_continuum_solution=False))

def measure(out,target,g):
    from measure_c02 import Stats
    mesh=c.read(out/'mesh-h0.json');case=c.read(out/'stages.json');ground=c.read(out/'ground.json')
    systems=[]
    for order in (4,6,8):
        if order!=4:contact_rule(order)
        systems.append(c.old.Water(mesh,ground,case).contact)
    # All original GP locations, separated by root face; no nearest matching.
    xyz=np.asarray(mesh['nodes']);coords=np.concatenate([np.array([p['position_m'] for p in s.points]) for s in systems]);roots=np.concatenate([np.array([p['face_index'] for p in s.points]) for s in systems])
    matrices=[]
    for order,s in zip((4,6,8),systems):
        mat=np.zeros((len(coords),len(s.points)))
        for f in (0,1):
            indices=np.array([i for i,p in enumerate(s.points) if p['face_index']==f]);ids=np.flatnonzero(roots==f)
            vertices=xyz[np.array(s.points[indices[0]]['nodes_a'])[:3]];x=m.natural(vertices,coords[ids])
            mat[np.ix_(ids,indices)]=interpolate(order,None,x)
        matrices.append(mat)
    folders=[out/f'h0-g{g}/independent']+[out/f'contact-q{q}-g{g}/independent' for q in (6,8)]
    stats=[Stats(),Stats()];byorder={};count=0
    for states in zip(*(c.old.states(p) for p in folders),strict=True):
        arrays=[]
        for mat,(row,fields) in zip(matrices,states):
            f={k:mat@fields[k] for k in m.CONTACT};f.update({k:fields[k] for k in m.MATERIAL+m.NODAL});f['open']=(f['gap']>0).astype(float);arrays.append(f)
        for j in (0,1):stats[j].add(arrays[j],arrays[j+1],states[0][0]['stage'])
        for k in arrays[0]:
            a=arrays[0][k]-arrays[1][k];b=arrays[1][k]-arrays[2][k];v=byorder.setdefault(k,[0.,0.,0.]);v[0]+=float((a*a).sum());v[1]+=float((b*b).sum());v[2]+=float((a*b).sum())
        count+=1
    report=dict(contact_orders=[4,6,8],common_states=count,pairs=[s.report() for s in stats],difference_directions={k:float(d/np.sqrt(a*b)) if a*b>0 else None for k,(a,b,d) in byorder.items()},
                point_positions_m=coords,point_root_face=roots,production_changed=False,independent_input_only=True,accepted=False,
                interpretation='Changes GP placement and quadrature together on the same mesh. Does not uniquely separate those two mechanisms. No observed h order from non-geometric quadrature sequence.')
    c.write(target/'contact-quadrature-diagnostic.json',report)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--mode',choices=['solve','measure','water'],required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--target',type=Path);p.add_argument('--order',type=int,choices=[6,8]);p.add_argument('--steps',type=int,default=256);a=p.parse_args()
    if a.mode=='solve':solve(a.output.resolve(),a.order,a.steps)
    elif a.mode=='measure':measure(a.output.resolve(),a.target.resolve(),a.steps)
    else:water(a.output.resolve(),a.target.resolve(),m.Maps([c.read(a.output/f'mesh-h{h}.json') for h in range(3)],c.read(a.output/'stages.json')))
