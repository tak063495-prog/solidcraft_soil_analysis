"""Audit sparse compatible MCC against independently accumulated dense solves."""
import argparse
import copy
import json
from pathlib import Path
import mcc_rate_fem as rate
import mcc_compatible_sparse_reference as sparse
import freeze_sources

np = sparse.np


class DenseCompatibleSystem(rate.System):
    def linear_solve(self, local, rhs):
        lookup = np.full(self.xyz.size, -1)
        lookup[self.free] = np.arange(len(self.free))
        indices = lookup[self.dofs]
        matrix = np.zeros((len(self.free), len(self.free)))
        for cell, stiffness in zip(indices, local):
            used = cell >= 0
            matrix[np.ix_(cell[used], cell[used])] += stiffness[np.ix_(used, used)]
        value = np.linalg.solve(matrix, rhs)
        assert np.linalg.norm(matrix@value-rhs)/max(np.linalg.norm(rhs), 1e-30) < 1e-10
        return value


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def verify(study, output):
    mesh, ground, case = [read(study/'inputs'/f'{name}.json') for name in ('mesh', 'ground', 'case')]
    original = np.asarray(read(study/'initial-400.json')['material_state'])
    source = freeze_sources.manifest()
    rng = np.random.default_rng(20261026)
    reports = []
    for directional in (False, True):
        current = copy.deepcopy(case)
        if directional:
            current['settings']['kinematics'] = 'explicit'
            current['boundary']['constraints'] = [dict(node_set=name, direction=direction, value_m=0.)
                for name, direction in (('xmin', [1., .2, .1]), ('xmax', [1., .2, .1]),
                                        ('ymin', [0., 1., .1]), ('ymax', [0., 1., .1]), ('bottom', [0., 0., 1.]))]
        for restricted in (False, True):
            active = np.ones(len(mesh['tetrahedra']), dtype=bool)
            if restricted:
                active[mesh['excavation_sets']['Pit']] = False
            expected = DenseCompatibleSystem(mesh, ground, active, current)
            actual = sparse.SparseCompatibleSystem(mesh, ground, active, current)
            for path in ('initial_elastic', 'nc_plastic', 'mixed_active_set'):
                if path == 'initial_elastic':
                    state = original[active].reshape(-1, 14).copy()
                    plastic = np.zeros(len(state), dtype=bool)
                else:
                    state = np.zeros((actual.weight.size, 14))
                    state[:, :3] = 100.
                    state[:, 6] = np.log(100.)
                    state[:, 7] = np.log(2.7)
                    plastic = np.ones(len(state), dtype=bool)
                    if path == 'mixed_active_set':
                        plastic[::3] = False
                force = rng.normal(size=actual.xyz.size)
                a = actual.rate(state, force, plastic)
                b = expected.rate(state, force, plastic)
                errors = {key: float(np.linalg.norm(x-y)/max(np.linalg.norm(y), 1e-30))
                          for key, x, y in zip(('displacement_rate', 'material_state_rate', 'loading_potential', 'yield_value'), a, b)}
                assert max(errors.values()) < 1e-10, errors
                reports.append(dict(directional_constraints=directional, restricted_domain=restricted,
                                    path=path, relative_errors=errors, solver=actual.solver_diagnostics()))
                print(json.dumps({k: reports[-1][k] for k in ('directional_constraints', 'restricted_domain', 'path', 'relative_errors')}), flush=True)
    assert source == freeze_sources.manifest(), 'source changed during verification'
    with (output/'rate-equivalence.json').open('x', encoding='utf-8') as stream:
        json.dump(dict(schema='soil-fem-compatible-sparse-rate-equivalence/1', passed=True,
                       goal_complete=False, p6_complete=False, spatial_convergence_accepted=False,
                       cases=reports, source_sha256=source), stream, indent=2, allow_nan=False)
    print('Sparse compatible rate verification passed', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('study', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    verify(args.study.resolve(), args.output.resolve())
