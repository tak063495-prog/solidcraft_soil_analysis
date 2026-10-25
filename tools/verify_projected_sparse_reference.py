"""Compare sparse three-field rates to independent dense elimination.

Audits full vectors/GP rates and both projection operators on actual CAD cells,
including domain restriction, material-split shared nodes and plastic tangents.
This is a linear-rate equivalence check, not a spatial-convergence acceptance.
"""
import argparse
import copy
import json
from pathlib import Path
import mcc_projected_volume_diagnostic as dense
from mcc_projected_sparse_reference import SparseProjectedSystem
import mcc_projected_sparse_reference as sparse
import freeze_sources

np = dense.np


def relative(left, right):
    return float(np.linalg.norm(left-right)/max(np.linalg.norm(right), 1e-30))


def verify(study, output):
    mesh, ground, case = [dense.read(study/'inputs'/f'{name}.json') for name in ('mesh', 'ground', 'case')]
    initial_file = study/'initial-400.json'
    initial = np.asarray(dense.read(initial_file)['material_state'])
    inputs = [study/'inputs'/f'{name}.json' for name in ('mesh', 'ground', 'case')]+[initial_file]
    bindings = {str(p): dense.sha(p) for p in inputs}
    sources = freeze_sources.manifest()
    rng = np.random.default_rng(20261025)
    reports = []
    for split in (False, True):
        current_mesh, current_ground = copy.deepcopy(mesh), copy.deepcopy(ground)
        if split:
            material = copy.deepcopy(ground['materials'][0])
            material['id'] = max(m['id'] for m in ground['materials'])+1
            material['mcc']['lambda'] *= 1.2
            current_ground['materials'].append(material)
            xyz = np.asarray(mesh['nodes'])
            centers = xyz[np.asarray(mesh['tetrahedra'])].mean(axis=1)
            current_mesh['material_ids'] = [material['id'] if x[0] > np.median(centers[:, 0]) else old
                                            for x, old in zip(centers, mesh['material_ids'])]
        for restricted in (False, True):
            active = np.ones(len(mesh['tetrahedra']), dtype=bool)
            if restricted:
                stage = next(s for s in case['stages'] if s['kind'] == 'excavate')
                for body in stage['bodies']:
                    active[mesh['excavation_sets'][body]] = False
            old = dense.ProjectedSystem(current_mesh, current_ground, active, case)
            new = SparseProjectedSystem(current_mesh, current_ground, active, case)
            assert old.keys == new.keys
            stress = rng.normal(size=(new.weight.size, 14))
            displacement = rng.normal(size=new.xyz.size)
            operators = dict(mass=relative(new.mass.toarray(), old.mass),
                             divergence=relative(new.divergence.toarray(), old.divergence),
                             strain=relative(new.strain(displacement), old.strain(displacement)),
                             pressure=relative(new.pressure(stress)[0], old.pressure(stress)[0]),
                             force=relative(new.force(stress), old.force(stress)))
            assert max(operators.values()) < 2e-12, operators
            for path in ('original_initial_elastic', 'normally_consolidated_plastic', 'mixed_active_set'):
                if path == 'original_initial_elastic':
                    state = initial[active].reshape(-1, 14).copy()
                    plastic = np.zeros(len(state), dtype=bool)
                else:
                    state = np.zeros((new.weight.size, 14))
                    state[:, :3] = 100.
                    state[:, 6] = np.log(100.)
                    state[:, 7] = np.log(2.7)
                    plastic = np.ones(len(state), dtype=bool)
                    if path == 'mixed_active_set':
                        plastic[::3] = False
                force = rng.normal(size=new.xyz.size)
                expected = old.rate(state, force, plastic)
                actual = new.rate(state, force, plastic)
                rates = {key: relative(a, b) for key, a, b in zip(
                    ('displacement_rate', 'material_state_rate', 'loading_potential', 'yield_value'), actual, expected)}
                assert max(rates.values()) < 1e-10, rates
                reports.append(dict(material_split=split, restricted_domain=restricted, path=path,
                                    operator_comparison=operators, rate_comparison=rates,
                                    operator_checks=new.checks, solver=new.solver_diagnostics()))
                print(json.dumps({k: reports[-1][k] for k in ('material_split', 'restricted_domain', 'path', 'rate_comparison')}), flush=True)
    assert sources == freeze_sources.manifest(), 'source changed during verification'
    assert bindings == {str(p): dense.sha(p) for p in inputs}
    report = dict(schema='soil-fem-projected-sparse-rate-equivalence/1', passed=True,
                  goal_complete=False, p6_complete=False, spatial_convergence_accepted=False,
                  cases=reports, input_sha256=bindings, source_sha256=sources,
                  limitations='Rate/operator equivalence only; no nonlinear cycle or spatial convergence claim')
    with (output/'rate-equivalence.json').open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
    print('Sparse/dense projected rate equivalence passed', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('study', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    verify(args.study.resolve(), args.output.resolve())
