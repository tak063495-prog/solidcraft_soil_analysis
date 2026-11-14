"""Re-evaluate two zero increments after an actually attained full release."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from study import Model, read, write


def run(args):
    args.output.mkdir(parents=True, exist_ok=False)
    saved = read(args.completed/'restart.json')
    c = read(args.completed/'controls.json')
    p = saved['payload']
    if p['coordinate'] != 1. or p['fractions'] != [1., 1.]:
        raise ValueError('full release was not attained')
    args.map = c['binding']['map']; args.order = c['binding']['order']
    args.material_steps = c['binding']['material_steps']
    m = Model(args, args.output)
    old = m.state(p['states'])
    records = []
    for i in range(2):
        force, tangent, new = m.evaluate(np.zeros(m.R.shape[1]), old, 'post-full-release-'+str(i))
        residual = float(np.linalg.norm(m.R.T@m.final-force)/max(np.linalg.norm(m.final), 1.))
        record = dict(coordinate=1., displacement_m=p['displacement_m'], states=m.serial(new),
                      contacts=m.contacts, zero_increment_exact=m.serial(new) == m.serial(old),
                      removed_material_frozen=m.frozen(new), contacts_exact=m.contacts == p['contacts'],
                      holding_load_kn=(m.load(1.)-m.final).tolist(), relative_force_residual=residual,
                      material_tangent=tangent.tolist())
        if not all(record[k] for k in ('zero_increment_exact', 'removed_material_frozen', 'contacts_exact')) or residual >= 1e-8:
            raise ValueError('post full-release field or equilibrium changed')
        records.append(record); old = new
    result = dict(passed=False, evaluator_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), model_source=m.source, source_checkpoint_sha256=hashlib.sha256((args.completed/'restart.json').read_bytes()).hexdigest(),
                  diagnostic_only=True, full_release_zero_increment_checks=records)
    write(args.output/'post-full-release.json', result)
    if args.map=='native':write(args.output/'material-processes.json',m.m.processes)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('completed', 'evidence', 'binary', 'reference', 'output'):
        p.add_argument('--'+name, type=Path, required=True)
    run(p.parse_args())
