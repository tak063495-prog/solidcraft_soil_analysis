"""Exact saved point-set audit; no additional equilibrium solve."""
import argparse
import hashlib
import json
from pathlib import Path
from a03follow_binding import study
from a03follow_report import physical_fingerprint
from a03follow_run import write


def point_sets(folder):
    mesh = study.read(folder/'mesh.json')
    run = study.read(folder/'native/run.json')
    elements = [tuple(sorted(cell)) for cell in mesh['tetrahedra10']]
    material = []
    for p in run['integration_points']:
        data = dict(p)
        owner = elements[data.pop('element_index')]
        material.append(((owner, data['local_index']), data))
    contact = []
    for p in run['joint_integration_points']:
        data = dict(p)
        data.pop('face_index')
        for side in ('a', 'b'):
            data['owner_'+side] = elements[data['owner_'+side]]
        key = (data['joint_id'], tuple(data['nodes_a']), tuple(data['nodes_b']), tuple(data['shape']))
        contact.append((key, data))
    material, contact = sorted(material), sorted(contact)
    assert len({k for k, _ in material}) == len(material) == 324
    assert len({k for k, _ in contact}) == len(contact) == 32
    return dict(material=material, contact=contact)


def audit(run, output):
    specs = [s for s in study.read(run/'plan.json')['specs'] if 'order' in s['roles']]
    baseline = next(s for s in specs if not any(s[k] for k in ('element', 'face', 'reverse')))
    original = point_sets(run/baseline['name'])
    rows = []
    for s in specs:
        actual = point_sets(run/s['name'])
        equal = {key: actual[key] == original[key] for key in original}
        assert all(equal.values())
        packed = json.dumps(actual, sort_keys=True, separators=(',', ':')).encode()
        rows.append(dict(path=s['name'], exact_physical_point_data_equal=equal,
                         canonical_point_data_sha256=hashlib.sha256(packed).hexdigest(),
                         physical_input_sha256=physical_fingerprint(run/s['name'])))
    assert len({r['canonical_point_data_sha256'] for r in rows}) == 1
    assert len({r['physical_input_sha256'] for r in rows}) == 1
    write(output, dict(diagnostic_only=True, new_equilibrium_solves=0, observations=rows,
                      changed_quadrature_point_sets_or_weights=False, physical_mapping_mismatches=0,
                      interpretation='Exact same weighted discrete physical problem. Differences after fixed-order aggregation are arithmetic-order/linear-solve roundoff sensitivity in this fixture; no claim about general plasticity or Windows bit identity.'))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    audit(a.run, a.output)
