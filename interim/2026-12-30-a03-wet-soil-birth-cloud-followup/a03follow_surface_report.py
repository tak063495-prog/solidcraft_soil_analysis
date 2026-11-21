"""Append the complete independent surface-rule control without replacing raw paths."""
import argparse
from pathlib import Path
import numpy as np
from a03follow_binding import study
from a03follow_report import load, values, raw_observations
from a03follow_controls import compare
from a03follow_run import write, sha


def assess(formal, surface, primary, output, vectors_output):
    result = study.read(primary)
    specs = [s for s in study.read(formal/'plan.json')['specs'] if 'order' in s['roles']]
    baseline = next(s for s in specs if not any(s[k] for k in ('element', 'face', 'reverse')))
    def rows(spec):
        return load(surface if spec['reverse'] else formal, spec, 'independent')
    vectors = {}
    selected = ('u', 'stress', 'pore_pressure', 'reaction', 'traction', 'endpoint_work', 'numerical_dissipation', 'load_work_cumulative')
    def diff(name, a, b):
        vectors[name] = {k: (values(a, k)-values(b, k)).tolist() for k in selected}
        return compare(a, b)
    full = dict(primary_matrix_Q_scope='contact Tri6 and GL25 face rule; dedicated hydrostatic rule retained',
                completed_Q_scope='all the above plus independent dedicated hydrostatic inner/outer ray rules',
                native_Q_scope='common face_points_order, including fully submerged hydrostatic faces',
                all_combinations={}, independent_effects={}, dedicated_hydrostatic_additional_effect={}, implementation={})
    base = rows(baseline)
    for spec in specs:
        name = spec['name']
        r = rows(spec)
        full['all_combinations'][name] = diff('all/'+name, r, base)
        full['implementation'][name] = diff('implementation/'+name, load(formal, spec, 'native'), r)
        if spec['reverse']:
            full['dedicated_hydrostatic_additional_effect'][name] = diff('hydro-only/'+name, r, load(formal, spec, 'independent'))
            result['raw'][name+'/independent-full-surface'] = raw_observations(surface, spec, 'independent')
    for factor in ('element', 'face', 'reverse'):
        pairs = []
        for spec in specs:
            if spec[factor]:
                other = next(s for s in specs if all(s[k] == (False if k == factor else spec[k]) for k in ('element', 'face', 'reverse')))
                pairs.append(dict(changed=spec['name'], baseline=other['name'],
                                  effect=diff(factor+'/'+spec['name'], rows(spec), rows(other))))
        full['independent_effects'][factor] = pairs
    driver = study.read(surface/'driver-result.json')
    result['ordering_full_surface'] = full
    result['additional_full_surface_paths'] = driver['paths']
    result['formal_paths_including_full_surface'] = result['formal_paths']+driver['fresh_paths']
    result['formal_states_including_full_surface'] = result['formal_states']+driver['fresh_states']
    result['total_fresh_global_paths_including_compatibility'] = result['formal_paths_including_full_surface']+result['compatibility_paths']+result['standard_cli_paths']
    result['total_fresh_global_states_including_compatibility'] = result['formal_states_including_full_surface']+result['compatibility_states']+result['standard_cli_states']
    result['primary_matrix_report_sha256'] = sha(primary)
    write(vectors_output, vectors)
    write(output, result)
    print(result['formal_paths_including_full_surface'], result['total_fresh_global_paths_including_compatibility'])


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--formal', type=Path, required=True)
    p.add_argument('--surface', type=Path, required=True)
    p.add_argument('--primary', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--vectors', type=Path, required=True)
    a = p.parse_args()
    assess(a.formal, a.surface, a.primary, a.output, a.vectors)
