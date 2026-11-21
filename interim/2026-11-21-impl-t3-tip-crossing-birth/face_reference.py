"""Dry elastic bulk and frictional face reference for the two archived construction inputs."""
from types import SimpleNamespace
import numpy as np
from scipy.optimize import root
import fem_reference as fem
from joint_ext_reference import elastic_operators, constraint_space
from joint_reference import System
from verify_high_order import physical
from verify_joint_coupling import contact
from wet_joint_cycle_reference import contact_work


def solve(mesh, ground, case):
    geo = ground['initialization']['geostatic']
    assert ground['initialization']['water_level_m'] is None
    joint = case['joint_system']['joints'][0]
    system = System(mesh,case['joint_system'])
    ops = elastic_operators(mesh,ground)
    cells = fem.cells(mesh)
    bodies = {r['id']:r['body'] for r in mesh['regions']}
    names = [bodies[i] for i in mesh['region_ids']]
    active = np.array([b not in case.get('initially_inactive_bodies',[]) for b in names])
    seed = np.array([geo['seed']['states'][b]['effective_stress_compression_positive_kpa'] for b in names])
    nd = len(mesh['nodes'])*3
    assert nd == 243 and len(ops) == 24 and len(system.points) == 32
    u = np.zeros(nd)
    origin = np.zeros((len(ops),nd))
    q = np.zeros_like(system.normal)
    jump = q.copy()
    rest = q.copy()
    work = np.zeros(len(q))
    pd = work.copy()
    numerical = work.copy()
    slip = q.copy()
    generation = np.zeros(len(q),dtype=int)
    born = np.zeros(len(q),dtype=bool)
    owner_a = np.array([p['owner_a'] for p in system.points])
    owner_b = np.array([p['owner_b'] for p in system.points])
    factor = geo['initial_boundary_factor']
    pressure = ground['initialization']['surface_surcharge_kpa']
    sample = [fem.volume(mesh,e) for e in range(len(ops))]
    states = []
    for stage in [dict(kind='initial',name='initial_equilibrated')]+case['stages']:
        kind = stage['kind']
        previous_active = active.copy()
        if kind == 'birth':
            assert stage['placement']=='cad_new_nodes'
            used = set(cells[active].ravel())
            placed = u.reshape(-1,3).copy()
            for n in set(range(len(placed)))-used:
                placed[n] = 0.
            for p in system.points:
                for a,b in zip(p['nodes_a'],p['nodes_b']):
                    if a in used and b not in used:
                        placed[b] = placed[a]
                    elif b in used and a not in used:
                        placed[a] = placed[b]
            u = placed.ravel()
            for e,b in enumerate(names):
                if b in stage['bodies']:
                    active[e] = True
                    seed[e] = stage['states'][b]['effective_stress_compression_positive_kpa']
                    origin[e] = u
        elif kind == 'excavate':
            for b in stage['bodies']:
                active[mesh['excavation_sets'][b]] = False
            pressure = stage['remaining_surface_pressure_kpa']
        elif kind == 'joint_death':
            born[:] = False
        live = active[owner_a] & active[owner_b]
        k = np.zeros((nd,nd))
        inside = np.zeros(nd)
        for e in np.flatnonzero(active):
            op,_,ids,ke = ops[e]
            k[np.ix_(ids,ids)] += ke
            inside[ids] -= np.einsum('qai,a,q->i',op['b'],seed[e],op['weight'])+ke@origin[e,ids]
        load,_ = physical(mesh,ground,active,pressure,case['pressure_on_original_surface_only'],joint_system=system)
        load = (load+factor*fem.conditions_load(mesh,case['boundary'],active)).ravel()
        if kind == 'joint_birth':
            new = live & ~born
            jump[new] = system.c[new]@u
            rest[new] = jump[new]
            q[new] = 0.
            work[new] = pd[new] = numerical[new] = 0.
            slip[new] = 0.
            generation[new] += 1
            if joint['birth_policy']=='transfer_force':
                # A-02 S1 has no pre-existing active contact at this event.
                assert not born.any()
                free,_ = constraint_space(mesh,case,active)
                pairs = sorted({(a,b) for i,p in enumerate(system.points) if new[i]
                                for a,b in zip(p['nodes_a'],p['nodes_b']) if a!=b})
                pm = np.zeros((nd,3*len(pairs)))
                for j,(a,b) in enumerate(pairs):
                    pm[3*a:3*a+3,3*j:3*j+3] = -np.eye(3)
                    pm[3*b:3*b+3,3*j:3*j+3] = np.eye(3)
                rhs = load-k@u-inside
                f = pm @ np.linalg.lstsq(free.T@pm,free.T@rhs,rcond=1e-12)[0]
                mat = np.einsum('qag,q->gqa',system.c[new],np.sqrt(system.weight[new])).reshape(nd,-1)
                x = np.linalg.lstsq(mat,f,rcond=1e-12)[0]
                assert np.max(np.abs(free.T@(f-rhs)))<=1e-10 and np.max(np.abs(mat@x-f))<=1e-10
                trial = x.reshape(-1,3)/np.sqrt(system.weight[new])[:,None]
                pn = -np.einsum('qa,qa->q',trial,system.normal[new])
                excess = np.linalg.norm(trial+pn[:,None]*system.normal[new],axis=1)-system.mu[new]*pn
                assert np.any((pn<0)|(excess>0))
                return states, dict(rejected=True,max_friction_excess_kpa=float(excess.max()))
            born[new] = True
        on = live & born
        t,offset = constraint_space(mesh,case,active,SimpleNamespace(points=[
            p for i,p in enumerate(system.points) if not born[i]]),True)
        fixed = factor*offset
        c = system.c[on]
        weights = system.weight[on]
        sub = SimpleNamespace(**{key:getattr(system,key)[on] for key in ('normal','kn','kt','mu')})
        cr = c@t
        h = k@u+inside-load+np.einsum('qag,qa,q->g',c,q[on],weights)
        if kind == 'excavate':
            # The archived S2 input removes its contact before excavating soil.
            assert not on.any()
            gone = np.zeros(nd)
            for e in np.flatnonzero(previous_active & ~active):
                op,_,ids,ke = ops[e]
                gone[ids] += ke@(u[ids]-origin[e,ids])-np.einsum('qai,a,q->i',op['b'],seed[e],op['weight'])
            h = previous_load-load-gone
        h.reshape(-1,3)[np.setdiff1d(np.arange(nd//3),cells[active].ravel())] = 0.
        previous_load = load
        reduced = t.T@k@t
        nstep = 1 if kind=='initial' else round(1/stage.get('initial_increment',case['settings']['initial_increment']))
        for i in range(1,nstep+1):
            external = load+(1-i/nstep)*h
            base = t.T@(k@fixed+inside-external)
            old = jump[on]-rest[on]
            def gp(v):
                return contact(sub,old,q[on],v,1)[0][:,0]
            def args(x):
                return c@fixed+cr@x-rest[on]
            def residual(x):
                return reduced@x+base+np.einsum('qai,qa,q->i',cr,gp(args(x)),weights)
            def jacobian(x):
                v = args(x)
                dq = np.stack([(gp(v+h)-gp(v-h))/(4e-9) for h in 2e-9*np.eye(3)],axis=2)
                return reduced+np.einsum('qai,qab,qbj,q->ij',cr,dq,cr,weights,optimize=True)
            if on.any():
                answer = root(residual,t.T@(u-fixed),jac=jacobian,tol=1e-11)
                x = answer.x
                assert np.linalg.norm(residual(x))/max(np.linalg.norm(external),1.)<=1e-10
                target = args(x)
                assert np.max(np.linalg.norm(target-old,axis=1))<=case['joint_system']['maximum_increment_m']
                plastic, num, dw = contact_work(sub,old,q[on],target,1)
                work[on] += dw
                pd[on] += plastic
                numerical[on] += num
                slip[on] += contact(sub,old,q[on],target,1)[1][:,0]
                q[on] = gp(target)
            else:
                x = np.linalg.solve(reduced,-base)
            u = fixed+t@x
            jump[on] = c@u
            gp_stress = np.concatenate([seed[e]-np.einsum('ab,qbi,i->qa',ops[e][1],sample[e]['b'],
                                                       u[ops[e][2]]-origin[e,ops[e][2]])
                                        for e in np.flatnonzero(active)])
            states.append(dict(stress=gp_stress,stage=stage['name'],progress=i/nstep,u=u.reshape(-1,3).copy(),
                               q=q.copy(),work=work.copy(),pd=pd.copy(),numerical=numerical.copy(),
                               slip=slip.copy(),reference=rest.copy(),generation=generation.copy(),active=active.copy(),live=on.copy(),
                               load=external.copy(),internal=k@u+inside+np.einsum('qag,qa,q->g',c,q[on],weights)))
    return states, dict(rejected=False)
