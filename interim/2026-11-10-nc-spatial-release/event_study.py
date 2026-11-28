"""Bounded spatial-group release of the EXISTING cut holding load.

No new physical load, mesh, constitutive law or acceptance threshold. These are
diagnostic load procedures, not production job commits or a stability verdict.
Independent fields use input-only RK4 history and independently assembled 3D B.
"""
import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [str(ROOT/'tools'), str(HERE.parent/'2026-11-09-nc-pseudo-arclength')]
from native_arc import Native
from independent_branches import StrainRate as OriginalStrainRate, ClosedStickSystem, expected_inputs, physical
from spatial_release_control import partition, fractions as group_fractions


from mcc_rate_event_floor import advance_rk4 as diagnostic_advance


class StrainRate(OriginalStrainRate):
    def integrate(self, state, strain, count):
        return diagnostic_advance(self, state, np.zeros(1), strain, count, minimum_event_step=1e-14)[0]


def read(p):
    return json.loads(p.read_text())


def write(p, x):
    p.write_text(json.dumps(x, indent=2, allow_nan=False)+'\n')


def append(p, x):
    with p.open('a') as f:
        f.write(json.dumps(x, allow_nan=False)+'\n')


def digest(x):
    return hashlib.sha256(json.dumps(x, sort_keys=True, allow_nan=False).encode()).hexdigest()


def fractions(s, order):
    return group_fractions(s)


class Model:
    def __init__(self, args, out):
        self.kind = args.map
        if self.kind == 'native':
            self.m = m = Native(args.evidence, args.binary, out/'material')
            run = read(args.evidence/'prior/fresh-head-midpoint/job/recovery-checkpoint-0.json')['run']
            unload = next(f['state'] for f in run['frames'] if f['name'] == 'unload')
            self.old = unload['material_states']
            self.u = np.array(unload['displacement_m'])
            self.contacts = unload['joints']
            self.R, self.B, self.w, self.pi = m.R, m.rb, m.w, m.pi
            self.active = np.array(m.active)
            self.final = m.load(1.)
            self.holding = m.holding.copy()
            m.cp = dict(m.cp, material_states=self.old)
            self.source = dict(run_checkpoint_sha256=hashlib.sha256((args.evidence/'prior/fresh-head-midpoint/job/recovery-checkpoint-0.json').read_bytes()).hexdigest(),
                               binary_sha256=hashlib.sha256(args.binary.read_bytes()).hexdigest())
        else:
            mesh, ground, case = expected_inputs()
            unload = read(args.reference/'unload.json')
            self.old = np.array(unload['material_state'])
            self.u = np.array(unload['displacement_m'])
            self.contacts = unload['contacts']
            self.active = np.array([r != next(x['id'] for x in mesh['regions'] if x['body'] == 'UpperClay') for r in mesh['region_ids']])
            self.m = m = ClosedStickSystem(mesh, ground, self.active, case)
            if m.live.any() or np.any(m.kjoint):
                raise ValueError('component study requires all contacts inactive after cut')
            self.R = np.zeros((m.xyz.size, len(m.free)))
            for j, g in enumerate(m.free):
                v = np.zeros(m.xyz.size); v[g] = 1.
                self.R[:, j] = m.lift(v)
            self.B = np.einsum('eqai,eij->eqaj', m.b, self.R[m.dofs]).reshape(-1, 6, len(m.free))
            self.w = m.weight.ravel()
            self.pi = np.concatenate([np.arange(i*27, (i+1)*27) for i in np.flatnonzero(self.active)])
            self.material = StrainRate(m.parameters)
            self.material_steps = args.material_steps
            self.final = physical(mesh, ground, self.active, 0., True, joint_system=m.joint)[0].reshape(-1)
            before = m.total_force(self.old[self.active].reshape(-1, 14), self.u.ravel())
            self.holding = before-self.final
            self.source = dict(inputs=[mesh, ground, case], unload_sha256=hashlib.sha256((args.reference/'unload.json').read_bytes()).hexdigest())
        self.initial_old = self.serial(self.old)
        nodes = m.mesh['nodes'] if self.kind == 'native' else mesh['nodes']
        self.horizontal, self.vertical, self.first_nodes = partition(nodes, self.holding, args.order)
        self.order = args.order

    def serial(self, x):
        return x if self.kind == 'native' else np.asarray(x).tolist()

    def state(self, x):
        return x if self.kind == 'native' else np.array(x)

    def load(self, s):
        h, v = fractions(s, self.order)
        return self.final+(1-h)*self.horizontal+(1-v)*self.vertical

    def evaluate(self, du, old, label):
        try:
            return self._evaluate(du, old, label)
        finally:
            if self.kind == 'native':
                # Trial IO is transient. Persist hashes and the actual assembled
                # linear systems/accepted states, without filling disk with IO.
                for p in sorted(self.m.out.glob('call-*.json')):
                    raw=p.read_bytes()
                    append(self.m.out/'io-sha256.jsonl', dict(name=p.name, sha256=hashlib.sha256(raw).hexdigest(), size=len(raw)))
                    p.unlink()

    def _evaluate(self, du, old, label):
        de = -np.einsum('qai,i->qa', self.B, du)
        if np.abs(de).sum(axis=1).max() > .001:
            raise ValueError('existing total increment strain bound')
        if self.kind == 'native':
            return self.m.evaluate(du, old, label)
        state = old[self.active].reshape(-1, 14)
        new = self.material.integrate(state, de, self.material_steps)
        D = np.zeros((len(state), 6, 6))
        for j in range(6):
            plus = de.copy(); minus = de.copy()
            plus[:, j] += 1e-8; minus[:, j] -= 1e-8
            D[:, :, j] = (self.material.integrate(state, plus, self.material_steps)[:, :6]-self.material.integrate(state, minus, self.material_steps)[:, :6])/(2e-8)
        full = old.copy(); full[self.active] = new.reshape(old[self.active].shape)
        K = np.einsum('qai,qab,qbj,q->ij', self.B, D, self.B, self.w, optimize=True)
        f = -np.einsum('qai,qa,q->i', self.B, new[:, :6], self.w)
        return f, K, full

    def stress(self, old):
        if self.kind == 'native':
            return np.array([old[i]['effective_stress_compression_positive_kpa'] for i in self.pi])
        return old[self.active].reshape(-1, 14)[:, :6]

    def frozen(self, old):
        if self.kind == 'native':
            return all(old[i] == self.initial_old[i] for i in range(len(old)) if i not in self.pi)
        return bool(np.array_equal(old[~self.active], np.array(self.initial_old)[~self.active]))


def solve(m, old, s, out, label, budget_end, initial=None):
    du = np.zeros(m.R.shape[1]) if initial is None else initial.copy()
    events = []
    try:
        for it in range(21):
            if time.monotonic() >= budget_end:
                raise ValueError('wall time budget; inconclusive')
            f, K, new = m.evaluate(du, old, label)
            rhs = m.R.T@m.load(s)-f
            rel = float(np.linalg.norm(rhs)/max(np.linalg.norm(m.load(s)), 1.))
            dx = np.linalg.solve(K, rhs)
            lr = float(np.linalg.norm(K@dx-rhs)/max(np.linalg.norm(rhs), 1e-30))
            append(out/'linear-systems.jsonl', dict(label=label, iteration=it, matrix=K.tolist(), rhs=rhs.tolist(), solution=dx.tolist(), true_relative_residual=lr))
            events.append(dict(iteration=it, relative_force_residual=rel, correction_norm_m=float(np.linalg.norm(dx)), true_linear_residual=lr, reduced_increment_m=du.tolist()))
            if lr > 1e-10:
                raise ValueError('existing true linear residual gate; inconclusive')
            if it >= 1 and rel < 1e-8 and np.linalg.norm(dx) < 1e-11:
                return new, du, events, None
            if it == 20:
                raise ValueError('bounded Newton limit; inconclusive')
            accepted = False
            for alpha in [1., .5, .25, .125, .0625, .03125, .015625]:
                candidate = du+alpha*dx
                if np.abs(np.einsum('qai,i->qa', m.B, candidate)).sum(axis=1).max() > .001:
                    continue
                ff, _, _ = m.evaluate(candidate, old, label+'-line-search')
                if np.linalg.norm(m.R.T@m.load(s)-ff) < np.linalg.norm(rhs) or np.linalg.norm(rhs) == 0:
                    du = candidate; accepted = True; break
            if not accepted:
                raise ValueError('bounded line search exhausted; inconclusive')
    except (ValueError, np.linalg.LinAlgError, subprocess.SubprocessError) as exc:
        return None, None, events, str(exc)


def run(args):
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic(); end = started+args.budget_seconds
    m = Model(args, out)
    sources = [ROOT/'tools/mcc_rate_event_floor.py', ROOT/'tools/spatial_release_control.py', Path(__file__), ROOT/'tools/mcc_rate_fem.py', ROOT/'tools/mcc_rate_rk4.py', ROOT/'tools/low_preload_nc_reference.py', HERE.parent/'2026-11-09-nc-pseudo-arclength/native_arc.py', HERE.parent/'2026-11-09-nc-pseudo-arclength/independent_branches.py']
    binding = dict(rk4_minimum_event_step=1e-14, map=args.map, order=args.order, initial_step=args.step_size, material_steps=args.material_steps, minimum_step=args.minimum_step,
                   source=digest(m.source), driver_sources={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    write(out/'controls.json', dict(binding=binding, force_tolerance=1e-8, true_linear_tolerance=1e-10,
          correction_check_m=1e-11, minimum_step=args.minimum_step, maximum_increment_strain_l1=.001,
          maximum_newton_iterations=20, budget_seconds=args.budget_seconds, input_only=args.map=='independent',
          final_load_kn=m.final.tolist(), holding_kn=m.holding.tolist(), horizontal_holding_kn=m.horizontal.tolist(),
          vertical_holding_kn=m.vertical.tolist(), first_group_nodes=m.first_nodes, grouping='spatial nodes, all XYZ components', final_load_exact=bool(np.array_equal(m.load(1.), m.final)),
          active_elements=m.active.tolist(), transform=m.R.tolist(), B=m.B.tolist(), weights=m.w.tolist(), active_gp=m.pi.tolist()))
    old, u, s, ds, index = m.old, m.u.copy(), 0., args.step_size, 0
    if args.resume:
        saved = read(args.resume)
        if saved['binding'] != binding or saved['payload_sha256'] != digest(saved['payload']):
            raise ValueError('checkpoint binding/payload mismatch')
        p = saved['payload']; old=m.state(p['states']); u=np.array(p['displacement_m'])
        s, ds, index = p['coordinate'], p['next_step'], p['step']+1
        if p['contacts'] != m.contacts:
            raise ValueError('contact restart mismatch')
    else:
        # A real perturbed lambda=0 Newton solve, holding load unchanged.
        direction = np.sin(np.arange(m.R.shape[1])+1.)
        initial = 1e-4*direction/np.linalg.norm(direction)
        old, du, ev, fail = solve(m, old, 0., out, 'lambda-zero', end, initial)
        append(out/'attempts.jsonl', dict(label='lambda-zero', target=0., accepted=fail is None, events=ev, stop_reason=fail))
        if fail:
            write(out/'report.json', dict(passed=False, stop_reason=fail, full_removal_reached=False, native_committed_increments=0)); return
        u += (m.R@du).reshape(u.shape)
    def save(label, du=None):
        p = dict(step=index-1, coordinate=s, fractions=list(fractions(s,args.order)), next_step=ds,
                 states=m.serial(old), displacement_m=u.tolist(), contacts=m.contacts,
                 removed_material_frozen=m.frozen(old), contact_history_frozen=True)
        stress=m.stress(old)
        f=-np.einsum('qai,qa,q->i',m.B,stress,m.w)
        p['reassembled_relative_force_residual']=float(np.linalg.norm(m.R.T@m.load(s)-f)/max(np.linalg.norm(m.load(s)),1.))
        if du is not None: p['reduced_increment_m']=du.tolist()
        append(out/'samples.jsonl',p)
        write(out/'restart.json',dict(binding=binding,payload=p,payload_sha256=digest(p)))
        if label: write(out/(label+'.json'),dict(binding=binding,payload=p,payload_sha256=digest(p)))
    save('initial')
    attempts=0; accepted=0; stop=None
    while s < 1 and accepted < args.count:
        if time.monotonic() >= end: stop='wall time budget; inconclusive';break
        corner=.5 if s < .5 and args.order!='proportional' else 1.
        target=min(s+ds,corner); before=digest(m.serial(old))
        new,du,events,fail=solve(m,old,target,out,'step-'+str(index),end)
        attempts+=1
        append(out/'attempts.jsonl',dict(step=index,start=s,target=target,step_size=ds,accepted=fail is None,
               events=events,stop_reason=fail,rejected_committed_state_unchanged=before==digest(m.serial(old))))
        if fail:
            if 'wall time' in fail: stop=fail;break
            if ds/2 < args.minimum_step: stop='minimum increment exhausted; inconclusive';break
            ds/=2;continue
        old=new;u+=(m.R@du).reshape(u.shape);s=target;index+=1;accepted+=1
        ds=min(args.step_size,ds*2)
        save('midpoint' if s==.5 else None,du)
    if s==1:
        ff,_,frozen=m.evaluate(np.zeros(m.R.shape[1]),old,'full-release-freeze')
        write(out/'full-release.json',dict(coordinate=s,material_zero_increment_exact=m.serial(frozen)==m.serial(old),
              contacts=m.contacts,removed_material_frozen=m.frozen(old),relative_force_residual=float(np.linalg.norm(m.R.T@m.final-ff)/max(np.linalg.norm(m.final),1.))))
        stop='full release diagnostic reached; acceptance evidence still required'
    elif stop is None: stop='bounded point count; inconclusive'
    if args.map=='native':write(out/'material-processes.json',m.m.processes)
    report=dict(passed=False,history_converged=False,static_limit_confirmed=False,production_adopted=False,
                native_committed_increments=0,map=args.map,order=args.order,accepted_increments=accepted,attempts=attempts,
                final_coordinate=s,final_fractions=list(fractions(s,args.order)),full_removal_reached=s==1,
                full_release_acceptance=False,restart_to_completion_verified=False,stop_reason=stop,elapsed_seconds=time.monotonic()-started)
    write(out/'report.json',report);print(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--map',choices=['independent'],required=True)
    p.add_argument('--order',choices=['left-first','center-first'],required=True)
    for n in ['evidence','binary','reference','resume']:p.add_argument('--'+n,type=Path)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--step-size',type=float,choices=[1/128,1/256],default=1/128)
    p.add_argument('--material-steps',type=int,choices=[4,8],default=4)
    p.add_argument('--minimum-step',type=float,choices=[2**-22,2**-28],default=2**-22)
    p.add_argument('--count',type=int,default=2048)
    p.add_argument('--budget-seconds',type=float,default=300)
    a=p.parse_args()
    if not 0<a.budget_seconds<=1800 or not 0<a.count<=4096:p.error('bounded controls required')
    run(a)
