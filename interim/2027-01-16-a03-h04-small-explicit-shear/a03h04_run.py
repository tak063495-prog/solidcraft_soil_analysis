"""Compute the single h=0.4 cross-control with the frozen 01-07 observers."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SHEAR = HERE.parent / '2027-01-07-a03-explicit-shear-initial-stress'
THICK = HERE.parent / '2027-01-06-a03-unequal-thickness-anisotropic-k'
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
saved_path = sys.path[:]
try:
    sys.path[:0] = [str(SHEAR), str(THICK), str(ROOT / 'tools')]
    import a03aniso_fixture as thickness
    import a03shear_fixture as shear
    import a03shear_run as engine
    from a03shear_binding import controls, load, reference, study, compare
finally:
    sys.path[:] = saved_path

CASES = ('zero', 'xz1m', 'xz1p')
BASE = '0048954815fd4f53310417e0c5a903c6ef99b913'


def plan():
    return [dict(name=f'h4-{s}-n{n}-gl{gl}', height=.4, stress=s, n=n, volume=gl,
                 minimum=1, reverse=False, element=False, face=False,
                 engines=['native', 'independent'] if gl == 5 else ['independent'])
            for s in CASES for n, gl in ((4, 5), (8, 5), (16, 5), (8, 7))]


def fixture(spec, seal):
    if spec not in plan():
        raise ValueError('outside h=0.4 / single xz pair')
    n, s = spec['n'], spec['stress']
    original, ground, case = shear.fixture(n, .5, s)
    prior = f'formal/h5-{s}-n{n}-gl5/'
    for name, data in zip(('mesh', 'ground', 'case'), (original, ground, case)):
        encoded = (json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False)+'\n').encode()
        if hashlib.sha256(encoded).hexdigest() != seal[prior+name+'.json']['sha256']:
            raise ValueError('01-07 input hash mismatch: '+prior+name)
    audit = shear.audit(original, ground, case)
    mesh, _, _ = thickness.fixture(n, .4, 'k55')
    geometry = thickness.audit(mesh, ground, case)
    if any(mesh[k] != original[k] for k in original if k not in ('nodes', 'regions')):
        raise ValueError('non-geometric mesh change')
    audit.update(height_m=.4, geometry=geometry, inherited_input_prefix=prior)
    return (mesh, ground, case), audit


def run(out):
    out = out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    specs = plan()
    engine.write(out/'plan.json', specs)
    seal = study.read(SHEAR/'archive-index.json')['entries']
    suffix = '.exe' if os.name == 'nt' else ''
    probe = ROOT/'target/release/examples'/('wet_birth_a03_probe'+suffix)
    solver = ROOT/'target/release'/('soil-fem-3d'+suffix)
    engine.write(out/'binaries.json', {p.name: engine.sha(p) for p in (probe, solver)})
    engine.write(out/'source-binding.json', {
        p.relative_to(ROOT).as_posix(): engine.sha(p) for p in sorted([
            *HERE.glob('*.py'), *SHEAR.glob('*.py'), *THICK.glob('*.py'),
            ROOT/'examples/wet_birth_a03_probe.rs', HERE/'SOURCES.ja.md'])})
    for spec in specs:
        folder = out/spec['name']
        inputs, audit = fixture(spec, seal)
        study.inputs(folder, inputs)
        engine.write(folder/'spec.json', spec)
        with controls(spec['volume'], False, 1):
            model = reference.Model(*inputs)
            audit['volume_points'] = [len(q['weight']) for q, *_ in model.ops]
            audit['output_points'] = [len(q['weight']) for q in model.sample]
            assert audit['volume_points'] == [spec['volume']**3]*12
            assert audit['output_points'] == [27]*12
            engine.write(folder/'pre-solve-seed-forces.json', shear.pre_solve_forces(inputs))
        engine.write(folder/'input-audit.json', audit)
    paths = []
    for spec in specs:
        folder = out/spec['name']
        for implementation in spec['engines']:
            if implementation == 'native':
                value = engine.native(folder, spec, probe, solver)
            else:
                proc = engine.command([sys.executable, SHEAR/'a03shear_run.py', '--independent', folder],
                                      folder, 'independent')
                value = dict(process=proc, states=4*spec['n'])
            paths.append(dict(path=spec['name'], engine=implementation, **value))
        if len(spec['engines']) == 2:
            data = [load(out, spec, e) for e in spec['engines']]
            engine.write(folder/'implementation-comparison.json', compare(*data))
        print(spec['name'], flush=True)
    engine.write(out/'driver-result.json', dict(passed=False, paths=paths, formal_paths=len(paths),
                                                formal_states=sum(x['states'] for x in paths), failures=[]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args().output)
