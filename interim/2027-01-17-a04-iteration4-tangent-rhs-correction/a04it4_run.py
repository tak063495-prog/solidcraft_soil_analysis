"""Observe the unmodified stopping paths and replay the native fourth state."""
import argparse
import inspect
import os
from functools import partial
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / '2027-01-12-a04-fixed-order-portable-reference'))
from a04ref_dry_report import joint, mapping
from a04follow_build import replace
from a04follow_common import np, read, rows, sha, write


def replay(folder, observed, model, evaluate, projected, t, matrix):
    it = next(r for r in observed if r['kind'] == 'iterate' and r['i'] == 4)
    u = np.asarray(it['u']).ravel()
    layout = read(folder/'native/run.json')['joint_integration_points']
    order, match = mapping(layout, model.system.points)
    old = [it['old'][i] for i in order]
    history = dict(rest=np.array([p['birth']['elastic_offset_m'] for p in old]),
                   jump=np.array([p['jump_m'] for p in old]), traction=np.array([p['traction_kpa'] for p in old]))
    counts = np.array(it['counts'])[order]
    project = partial(projected, joint=history, counts=counts)
    assert np.all(np.isin(t, [-1.,0.,1.])) and np.all((t != 0).sum(axis=0) == 1)
    x = t.T @ u
    fixed = u - t @ x
    assert np.array_equal(fixed + t @ x, u)
    constant = t.T @ (model.k @ fixed + model.seed)
    constant -= t.T @ it['target']
    evaluate_at = partial(evaluate, fixed=fixed, constant=constant, projected=project)
    residual, jac = evaluate_at(x, True)
    force, candidate, target, traction, *_ = evaluate_at(x)
    assert np.array_equal(force,residual) and np.array_equal(candidate,u)
    gap = ((target-history['rest'])*model.system.normal).sum(axis=1)
    assert np.all(gap > 0) and not np.any(traction)
    dx = np.linalg.solve(jac,-residual)
    write(folder/'independent-fixed.json', dict(matrix=jac, rhs=-residual, dx=dx, basis=t,
        u=u, full_step=fixed+t@(x+dx), contact_tangent=jac-matrix, gap=gap, traction=traction,
        point_order=order, mapping=match, native_status=[it['points'][i]['status'] for i in order]))


def instrument(folder, ref, model, native):
    # Only the isolated function receives optional frozen inputs for read-only replay.
    src = inspect.getsource(ref.run)
    src = replace(src, 'def projected(target):','def projected_for(target, *, joint, counts):')
    mark = '            cr=s.c@soil_t'
    src = replace(src, mark, '            def projected(target):\n                return projected_for(target,joint=joint,counts=counts)\n'+mark)
    src = replace(src, 'def evaluate(x,tangent=False):',
                  'def evaluate(x,tangent=False, *, fixed=fixed, constant=constant, projected=projected):')
    mark = '                if relative<1e-10:'
    src = replace(src, mark, '''                selected = phase == 2 and progress == .25
                if selected:
                    if iteration == 0:
                        replay(folder,native,m,evaluate,projected_for,soil_t,matrix)
                    path.append(dict(i=iteration, u=candidate.copy(), rhs=-residual,
                                     relative=float(relative), counts=counts.copy()))
''' + mark)
    mark = '                for cut in range(18):'
    src = replace(src, mark, '''                if selected:
                    path[-1].update(matrix=jac.copy(),dx=correction.copy())
''' + mark)
    path = []
    namespace = dict(vars(ref), path=path, native=native, folder=folder, replay=replay)
    exec(compile(src, '<a04it4-independent-observer>', 'exec'),namespace)
    result = namespace['run'](model,4,lambda *_:None)
    assert result['complete']
    write(folder/'independent-path.json',path)


def run(out, probe):
    out.mkdir(parents=True,exist_ok=False)
    study, ref, _ = joint.binding('dry')
    for name in ('000','010'):
        folder = out/name
        folder.mkdir()
        mesh, ground, case = study.fixture('stress_free',4,True)
        if name == '010':
            mesh['face_sets'] = {k:list(reversed(v)) for k,v in mesh['face_sets'].items()}
        paths = [folder/(key+'.json') for key in ('mesh','ground','case')]
        for path,value in zip(paths,(mesh,ground,case),strict=True):
            write(path,value)
        env = dict(os.environ,A04IT4_TRACE=str(folder/'native.jsonl'))
        env.pop('A04FOLLOW_Q',None)
        subprocess.run([probe,*paths,folder/'native','--whole'], env=env, check=True)
        model = ref.Model(*[read(p) for p in paths])
        assert len(mesh['tetrahedra10'])==24 and len(model.system.points)==32 and not (model.frame or model.tendon)
        observed = list(rows(folder/'native.jsonl'))
        instrument(folder,ref,model,observed)
        write(folder/'binding.json',dict(input_sha256=[sha(p) for p in paths],probe_sha256=sha(probe)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--probe',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args = parser.parse_args()
    run(args.output.resolve(),args.probe.resolve())
