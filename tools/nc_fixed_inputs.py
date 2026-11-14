"""Semantic comparison for the two original proportional NC diagnostic inputs.

Only documented native serde defaults are supplied, at schema-specific paths.
Unknown fields, nonzero gradients and all physical/control values remain in the
comparison. No float tolerance, projection, or saved-state rewriting is used.
"""
from copy import deepcopy
import math


def canonical_inputs(mesh, ground, case):
    mesh, ground, case = deepcopy((mesh, ground, case))
    mesh.setdefault('excavation_sets', {})
    case['settings'].setdefault('linear_solver', {'method': 'dense_lu'})

    def loads(value, constraints=False):
        value.setdefault('nodal_loads', [])
        value.setdefault('surface_loads', [])
        if constraints:
            for constraint in value['constraints']:
                constraint.setdefault('gradient', [0., 0., 0.])
                constraint.setdefault('reference_point_m', [0., 0., 0.])
        for load in value['surface_loads']:
            if load['kind'] in ('traction', 'pressure'):
                zero = [[0., 0., 0.] for _ in range(3)] if load['kind'] == 'traction' else [0., 0., 0.]
                load.setdefault('gradient_kpa_per_m', zero)
                load.setdefault('reference_point_m', [0., 0., 0.])

    if 'boundary' in case:
        loads(case['boundary'], constraints=True)
    geo = ground['initialization'].get('geostatic', {})
    if 'base_loads' in geo:
        loads(geo['base_loads'])
    return mesh, ground, case


def equivalent(a, b):
    """JSON equality with exact finite numbers, without bool/number aliasing."""
    if isinstance(a, bool) or isinstance(b, bool):
        return type(a) is type(b) and a == b
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return math.isfinite(a) and math.isfinite(b) and a == b
    if type(a) is not type(b):
        return False
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(equivalent(a[k], b[k]) for k in a)
    if isinstance(a, (list, tuple)):
        return len(a) == len(b) and all(equivalent(x, y) for x, y in zip(a, b))
    return a == b


def fixed_inputs_match(actual, expected):
    return equivalent(canonical_inputs(*actual), canonical_inputs(*expected))


def differences(actual, expected, path=''):
    """Report exact schema/value differences for an audit, without mutating data."""
    if equivalent(actual, expected):
        return []
    if isinstance(actual, dict) and isinstance(expected, dict):
        result = []
        for key in sorted(actual.keys() | expected.keys()):
            here = path+'/'+key
            if key not in actual or key not in expected:
                result.append(dict(path=here, actual_present=key in actual, expected_present=key in expected,
                                   actual=actual.get(key), expected=expected.get(key)))
            else:
                result.extend(differences(actual[key], expected[key], here))
        return result
    if isinstance(actual, (tuple, list)) and isinstance(expected, (tuple, list)) and len(actual) == len(expected):
        return [row for i, (a, b) in enumerate(zip(actual, expected))
                for row in differences(a, b, path+'/'+str(i))]
    return [dict(path=path, actual=actual, expected=expected)]


def audit_history(folder):
    import json
    from probe_nc_midpoint import inputs
    mesh = json.loads((folder/'mesh.json').read_text(encoding='utf-8'))
    run = json.loads((folder/'run.json').read_text(encoding='utf-8'))
    actual = (mesh, run['ground'], run['case'])
    increment = actual[2]['stages'][2]['initial_increment']
    expected = inputs('midpoint')
    expected[2]['stages'][2]['initial_increment'] = increment
    allowed = increment in (1/128, 1/256)
    return dict(schema='nc-fixed-input-audit/1', allowed_initial_increment=allowed, initial_increment=increment,
                raw_equal=equivalent(actual, expected),
                fixed_inputs_match=allowed and fixed_inputs_match(actual, expected),
                raw_differences=differences(actual, expected),
                canonical_differences=differences(canonical_inputs(*actual), canonical_inputs(*expected)),
                physical_float_tolerance=0., saved_inputs_modified=False, numerical_acceptance=False)


if __name__ == '__main__':
    import argparse
    import json
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--history', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = audit_history(args.history)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(fixed_inputs_match=report['fixed_inputs_match'], raw_equal=report['raw_equal'],
                         canonical_difference_count=len(report['canonical_differences']))))
    raise SystemExit(0 if report['fixed_inputs_match'] else 1)
