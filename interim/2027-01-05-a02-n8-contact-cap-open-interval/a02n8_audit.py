"""Integrity of new fields and GP correspondence; not another solve or cause claim."""
import argparse
import importlib
from pathlib import Path
import numpy as np
from a02n8_binding import ROOT, fields, reference, norm, study
from a02n8_common import write, sha


def audit(run, output):
    paths = study.read(run/'paths.json')
    geometry = importlib.import_module('a02follow_binding').geometry
    results, solve_bindings = {}, []
    names = ('a02n8_run.py', 'a02n8_reference.py', 'a02n8_binding.py',
             'a02follow_policies.py', 'a02cause_binding.py', 'reference.py', 'fields.py', 'study.py')
    for p in paths:
        folder = run/p['path']
        assert p['complete'] and p['exit_code'] == 0 and not p['timed_out']
        binding = study.read(folder/'binding.json')
        source = {name: digest for name, digest in binding['source_sha256'].items() if Path(name).name in names}
        assert all(sha(ROOT/name) == digest for name, digest in source.items())
        solve_bindings.append(source)
        assert all(sha(folder/name) == digest for name, digest in binding['input_sha256'].items())
        model = reference.Model(*(study.read(folder/(name+'.json')) for name in ('mesh', 'ground', 'case')))
        maximum = 0.
        states = 0
        for row, f in fields.load_fields(folder/(p['engine']+'-fields.jsonl')):
            assert np.all(f['material_return'] == 0), 'elastic bulk restriction'
            assert len(f['contact_status']) == 32
            basis = model.spaces[not bool(np.all(f['joint_active']))][0]
            relative = norm(basis.T@(f['internal']-f['applied_load']).ravel())/max(norm(f['applied_load']), 1.)
            maximum = max(maximum, relative)
            states += 1
        results[p['path']] = dict(accepted_states=states, postsolve_reassembled_force_residual=maximum,
                                  gp_correspondence=geometry.audit(folder) if p['engine'] == 'native' else None)
    assert all(v == solve_bindings[0] for v in solve_bindings)
    write(output, dict(paths=results, solve_source_sha256=solve_bindings[0],
        solver_sources_and_prescribed_inputs_unchanged=True,
        analysis_bindings_are_versioned_separately=True,
        postsolve_integrity_audit_not_additional_fresh_solve=True, a02_passed=False))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    audit(a.run.resolve(), a.output.resolve())
