"""Finish the independent Q control through the dedicated hydrostatic ray rule."""
import argparse
from contextlib import contextmanager
import json
from pathlib import Path
import sys
from a03follow_binding import HERE, study, reference
from a03follow_controls import controls, Sink, layout
from a03follow_run import command, write, sha
import hydrostatic_reference as hydro


@contextmanager
def hydro_order(reverse):
    original = hydro.rule
    def rule(n):
        x, w = original(n)
        return x[::-1].copy(), w[::-1].copy()
    if reverse:
        hydro.rule = rule
    try:
        yield
    finally:
        hydro.rule = original


def worker(folder):
    spec = study.read(folder/'spec.json')
    with controls(5, True, 1), hydro_order(True):
        model = reference.Model(*[study.read(p) for p in study.paths(folder)])
        write(folder/'independent-layout.json', layout(model))
        with (folder/'independent-fields.jsonl').open('x', encoding='utf-8') as f, (folder/'independent-prepared.jsonl').open('x', encoding='utf-8') as p:
            result = reference.solve(model, spec['n'], Sink(f), Sink(p, work=False))
    result.update(dedicated_hydrostatic_rule_reversed=True, new_physical_input=False)
    write(folder/'independent-result.json', result)


def run(formal, out):
    out = out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    specs = [s for s in study.read(formal/'plan.json')['specs'] if 'order' in s['roles'] and s['reverse']]
    binding = {p.name: sha(p) for p in sorted(HERE.glob('*.py'))}
    write(out/'source-binding.json', binding)
    source = out/'source-snapshot'
    source.mkdir()
    import shutil
    for p in HERE.glob('*.py'):
        shutil.copy2(p, source/p.name)
    paths = []
    for spec in specs:
        folder = out/spec['name']
        study.inputs(folder, [study.read(p) for p in study.paths(formal/spec['name'])])
        write(folder/'spec.json', spec)
        process = command([sys.executable, HERE/'a03follow_surface.py', '--worker', folder], folder, 'independent')
        paths.append(dict(path=spec['name'], engine='independent', states=4*spec['n'], process=process,
                          reversed_contact_and_GL_faces=True, reversed_hydrostatic_inner_and_outer_ray_points=True))
    assert binding == {p.name: sha(p) for p in sorted(HERE.glob('*.py'))}
    write(out/'driver-result.json', dict(paths=paths, fresh_paths=len(paths), fresh_states=sum(p['states'] for p in paths),
                                       numerical_acceptance=False, partial_Q_control_retained_in_primary_matrix=True))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--formal', type=Path)
    p.add_argument('--output', type=Path)
    p.add_argument('--worker', type=Path)
    a = p.parse_args()
    if a.worker:
        worker(a.worker)
    else:
        run(a.formal, a.output)
