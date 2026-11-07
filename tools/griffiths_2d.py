"""Pinned, measured 661-node displacement comparison; never runs a slope solve.

Archived controls are remeasured separately, never fitted as an in-plane series.
New series require the matching reduced material, completed production results,
and the existing independent mechanics/VTU audit with the producing binary.
"""
import argparse
import contextlib
import hashlib
import json
from pathlib import Path
import tempfile
import zipfile
import numpy as np
from plastic_spatial_convergence import ROOT, read, write, sha
from assess_slope_h_series import contains_input

REFERENCE_DIGEST = '29822140bcfb32da9954d75db23062b65ea28f127630971475ce5c2ab0cb9715'
GEOMETRY_DIGEST = '263ac09491784831796a5e49de04ceae434ed23e51f800f617708c54f87a73bf'
DEPTHS = ('0.0', '0.5', '1.0')
LIMIT = .05  # Original displacement_sample_relative_L2; no CLI override.
ARCHIVES = {
    'interim/2026-10-30-Root635-CI-and-refinement-controls/slope/xy1-ny1/all661-nodes-and200-prisms.json':
        'c898a7305933820e866e1dbcdce9b0adced30b0a762ed537a1b45db032c6491f',
    'interim/2026-10-30-Root635-CI-and-refinement-controls/slope/xy1-ny2/all661-nodes-and200-prisms.json':
        '9642b8876e9eb35bf513f63bc4437f0319a90078d9281ab46ece3446608724fb',
    'interim/2026-10-30-Root635-CI-and-refinement-controls/slope/xy1-ny1-global64/all661-nodes-and200-prisms.json':
        '57028b96fca596626bd2c30d5a62e98854927ecbd156a564b72a5508b5b9e927',
    'interim/2026-10-30-async-UI-and-independent-SRM-bound/slope/material-half1200/all661-nodes-and200-prisms.json':
        'd55c2e715e6bca1a83306cf2d600d6c838f043b29374a9a4132e0cd39d786f61',
}


def digest(value):
    """Pin JSON content, allowing checkout CRLF and insignificant key ordering."""
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


def load_reference(path):
    if path is None:
        raise ValueError('original 661-node reference file required; no 5% verdict without data')
    value = read(path)
    if digest(value) != REFERENCE_DIGEST:
        raise ValueError('original 661-node reference changed or unverified')
    return value


def displacement_metrics(reference, samples):
    """Unweighted Frobenius L2 of all 661 x/z pairs, as in the original protocol.

    This numerical kernel also accepts manufactured fields for tests. Production
    entrypoints separately require the pinned actual reference and provenance.
    """
    target = np.asarray(reference, float)
    if target.shape != (661, 2) or not np.isfinite(target).all():
        raise ValueError('661 finite x/z reference pairs required')
    norm = float(np.linalg.norm(target))
    maximum = float(np.linalg.norm(target, axis=1).max())
    if norm <= 0 or not np.isfinite(norm):
        raise ValueError('nonzero finite reference norm required')
    if set(samples) != set(DEPTHS):
        raise ValueError('all three front/middle/back samples required')
    fields = {key: np.asarray(samples[key], float) for key in DEPTHS}
    if any(v.shape != (661, 2) or not np.isfinite(v).all() for v in fields.values()):
        raise ValueError('661 finite x/z pairs at each depth required')
    relative = {key: float(np.linalg.norm(v-target)/norm) for key, v in fields.items()}
    max_error = max(float(np.linalg.norm(v-target, axis=1).max()/maximum) for v in fields.values())
    variation = max(float(np.linalg.norm(a-b, axis=1).max()/maximum)
                    for a in fields.values() for b in fields.values())
    if not np.isfinite([*relative.values(), max_error, variation]).all():
        raise ValueError('nonfinite comparison metrics')
    return dict(displacement_sample_relative_L2_per_depth=relative,
                displacement_sample_relative_L2_limit=LIMIT,
                five_percent_verdict='pass' if max(relative.values()) <= LIMIT else 'fail',
                maximum_displacement_error_over_reference_max=max_error,
                depth_displacement_variation_over_max=variation,
                maximum_displacement_error_passed=max_error <= .1,
                depth_variation_passed=variation <= .01)


def sample_nodes(mesh, reference):
    xyz = np.asarray(mesh['nodes'], float)
    lookup = {tuple(np.round(row, 12)): i for i, row in enumerate(xyz)}
    if len(lookup) != len(xyz):
        raise ValueError('coincident distinct mesh nodes')
    maps = {}
    for depth in DEPTHS:
        indices = []
        for _, x, z, _, _ in reference['nodes']:
            point = np.array([x, float(depth), z])
            index = lookup.get(tuple(np.round(point, 12)))
            if index is None or np.max(abs(xyz[index]-point)) > 1e-10:
                raise ValueError('missing original 2D node at a required depth; no interpolation')
            indices.append(index)
        maps[depth] = indices
    return maps


def base_report(reference_path):
    return dict(reference_file_sha256=sha(reference_path), reference_content_sha256=REFERENCE_DIGEST,
                comparison_measured=True, new_slope_solve_performed=False,
                p6_complete=False, general_3d_spatial_convergence_demonstrated=False,
                griffiths_slope_five_percent_closed=False,
                fs=dict(status='not_attempted', value=None, failure_upper_bound=None))


def compare_archives(reference_path, output):
    reference = load_reference(reference_path)
    cases = []
    for name, expected in ARCHIVES.items():
        data = read(ROOT/name)
        if digest(data) != expected:
            raise ValueError('archived displacement data changed: '+name)
        if data['source2d_u_xz_m'] != reference['displacement_xz_m']:
            raise ValueError('archived source differs from pinned original reference')
        cases.append(dict(source=name, source_content_sha256=expected, in_plane_h_series=False,
                          **displacement_metrics(reference['displacement_xz_m'], data['native3d_u_xz_m'])))
    report = dict(schema='soil-fem-griffiths-archived-2d-comparison/1', **base_report(reference_path),
                  cases=cases, continuum_verdict='indeterminate',
                  note='Remeasured retained actual fields only; no new h-series solve or convergence claim.')
    write(output, report)
    return report


def verify_source(reference_path, fixtures_path, archive_path, output):
    """Reproduce the actual-data provenance check without executing Excel/FEM."""
    reference = load_reference(reference_path); source = reference['source']
    if sha(fixtures_path) != source['fixtures_sha256'] or sha(archive_path) != source['archive_sha256']:
        raise ValueError('pinned original fixture or archive bytes changed')
    fixtures = read(fixtures_path)
    fixture = next(v for v in fixtures['fixtures'] if v['name'] == 'griffiths_1999_coarse')
    if [row[:5] for row in fixture['nodes']] != reference['nodes']:
        raise ValueError('original node IDs, coordinates or supports differ')
    with zipfile.ZipFile(archive_path) as archive:
        if len(archive.namelist()) != len(set(archive.namelist())) or archive.testzip() is not None:
            raise ValueError('invalid original archive inventory or CRC')
        raw = archive.read(source['member'])
        if hashlib.sha256(raw).hexdigest() != source['result_sha256']:
            raise ValueError('original result bytes changed')
        original = json.loads(raw)
        displacement = np.asarray(original['displacements_flat']).reshape(-1, 3)
        if (not np.array_equal(displacement[:, 0], np.arange(1, 662))
                or displacement[:, 1:].tolist() != reference['displacement_xz_m']
                or dict(zip(original['metric_names'], original['metrics'], strict=True))['fss']
                != reference['material_reduction_factor']):
            raise ValueError('original 661 displacements or fixed factor differ')
    report = dict(schema='soil-fem-griffiths-reference-source-verification/1', passed=True,
                  original_661_nodes_and_displacements_exactly_equal=True,
                  source=source, reference_content_sha256=REFERENCE_DIGEST,
                  new_excel_run=False, new_slope_solve_performed=False,
                  comparison_measured=False, griffiths_slope_five_percent_closed=False, p6_complete=False)
    write(output, report)
    return report


def validate_run(run, mesh, ground, case, solver_digest):
    if (run.get('complete') is not True or run.get('geostatic', {}).get('complete') is not True
            or run['geostatic'].get('progress') != 1. or len(run.get('frames', [])) != 2):
        raise ValueError('incomplete slope solve; not convergence, collapse, or an Fs upper bound')
    if any(not contains_input(run.get(name), value) for name, value in [('ground', ground), ('case', case)]):
        raise ValueError('result ground or case differs from prepared inputs')
    if run.get('provenance', {}).get('solver_binary_sha256') != solver_digest:
        raise ValueError('producing solver binary changed or missing')
    geo = run['geostatic']; frame = run['frames'][0]
    # The production solver resets stage_progress and release-load bookkeeping
    # after initialization. Compare the total displacement, not that metadata.
    if geo['last_committed']['displacement_m'] != frame['state']['displacement_m']:
        raise ValueError('geostatic total displacement differs from initial frame')
    u = np.asarray(geo['last_committed']['displacement_m'], float)
    if u.shape != (len(mesh['nodes']), 3) or not np.isfinite(u).all() or np.any(u[:, 1] != 0.):
        raise ValueError('finite full displacement with every transverse DOF fixed required')
    if frame['state']['active_elements'] != [True]*len(mesh['tetrahedra']):
        raise ValueError('original full slope domain required')
    applied = np.asarray(frame['applied_force_kn'], float)
    reaction = np.asarray(frame['support_reaction_kn'], float)
    errors = [frame['relative_force_residual'], geo['final_relative_physical_force_residual'],
              float(np.linalg.norm(applied+reaction)/4400.),
              float(np.linalg.norm(applied-[0., 0., -4400.])/4400.),
              frame['maximum_normalized_yield'], frame['maximum_absolute_engineering_strain_from_seed'],
              geo['maximum_absolute_engineering_strain_from_seed']]
    if (not np.isfinite(errors).all() or any(v < 0 for v in errors[:4])
            or any(v > 1e-8 for v in errors[:4]) or errors[4] > 1e-7
            or any(v < 0 or v > .05 for v in errors[5:])):
        raise ValueError('original force/yield/strain physical audit gates failed')
    if any(a.get('committed_state_unchanged_on_rejection') is not True
           for a in geo['attempts'] if not a['accepted']):
        raise ValueError('rejected history was not preserved')
    return u


def independent_audit(sub, solver):
    """Run the existing verifier; no optional bypass or alternative acceptance."""
    from export_stages import export
    from verify_stage_package import verify
    with tempfile.TemporaryDirectory(prefix='griffiths-audit-') as temp:
        exports = Path(temp)/'vtu'
        with (sub/'griffiths-independent-audit.log').open('x', encoding='utf-8') as log:
            with contextlib.redirect_stdout(log):
                export(sub/'mesh.json', sub/'result.json', exports)
                verify(sub/'mesh.json', sub/'ground.json', sub/'result.json', exports, solver=solver)


def compare_series(folder, reference_path, output, solver):
    from prepare_slope_h_series import prepare
    reference = load_reference(reference_path)
    folder = Path(folder); solver = Path(solver); manifest = read(folder/'preparation.json')
    if (manifest.get('material_reduction_factor') != reference['material_reduction_factor']
            or manifest.get('reference_2d_content_sha256') != REFERENCE_DIGEST):
        raise ValueError('matching reduced 2D reference preparation required; F=1 cannot be compared')
    # Regenerate all inputs: rehashing a modified manifest must not authorize a
    # different physical problem, level order, depth, or 661-node sample map.
    with tempfile.TemporaryDirectory(prefix='griffiths-inputs-') as temp:
        expected = prepare(Path(temp)/'expected', manifest['levels'][0]['depth_divisions'], reference_path)
        if manifest != expected:
            raise ValueError('preparation manifest changed')
        for record in expected['levels']:
            for name, expected_sha in record['files_sha256'].items():
                if sha(folder/record['folder']/name) != expected_sha:
                    raise ValueError('prepared input changed: '+name)
    solver_sha = sha(solver); levels = []; files = {}; inputs = []
    # Preflight every level before exporting or auditing any one of them.
    for record in manifest['levels']:
        sub = folder/record['folder']
        mesh, ground, case, run = [read(sub/f'{name}.json') for name in ('mesh', 'ground', 'case', 'result')]
        u = validate_run(run, mesh, ground, case, solver_sha)
        mapping = read(sub/'reference-map.json')
        if mapping != sample_nodes(mesh, reference):
            raise ValueError('original node sample mapping changed')
        samples = {depth: u[indices][:, [0, 2]] for depth, indices in mapping.items()}
        levels.append(dict(folder=record['folder'], elements=record['elements'],
                           **displacement_metrics(reference['displacement_xz_m'], samples)))
        inputs.append(sub)
        for name in (*record['files_sha256'], 'result.json'):
            files[str(sub/name)] = sha(sub/name)
    for sub in inputs:
        independent_audit(sub, solver)
    if sha(solver) != solver_sha or any(sha(Path(p)) != h for p, h in files.items()):
        raise ValueError('solver, inputs or results changed during comparison')
    load_reference(reference_path)
    report = dict(schema='soil-fem-griffiths-h-series-2d-comparison/1', **base_report(reference_path),
                  levels=levels, files_sha256=files, solver_binary_sha256=solver_sha,
                  completed_user_pc_h_series_measured=True, independent_audit_passed=True,
                  all_levels_five_percent_passed=all(v['five_percent_verdict'] == 'pass' for v in levels),
                  continuum_verdict='indeterminate',
                  note='Discrete 661-node comparison only. No continuum field extrapolation, Fs or general 3D acceptance.')
    write(output, report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('archives', 'series', 'verify-source'))
    parser.add_argument('output', type=Path)
    parser.add_argument('--reference-2d', required=True, type=Path)
    parser.add_argument('--folder', type=Path)
    parser.add_argument('--solver', type=Path)
    parser.add_argument('--fixtures', type=Path)
    parser.add_argument('--archive', type=Path)
    args = parser.parse_args()
    if args.mode == 'series' and (args.folder is None or args.solver is None):
        parser.error('series requires --folder and --solver')
    if args.mode == 'verify-source' and (args.fixtures is None or args.archive is None):
        parser.error('verify-source requires --fixtures and --archive')
    if args.output.exists():
        parser.error('output already exists')
    try:
        if args.mode == 'verify-source':
            verify_source(args.reference_2d, args.fixtures, args.archive, args.output)
            return
        report = (compare_archives(args.reference_2d, args.output) if args.mode == 'archives' else
                  compare_series(args.folder, args.reference_2d, args.output, args.solver))
    except (ValueError, KeyError, OSError, AssertionError, zipfile.BadZipFile) as error:
        parser.exit(2, f'comparison refused: {error}\n')
    passed = (all(v['five_percent_verdict'] == 'pass' for v in report['cases']) if args.mode == 'archives'
              else report['all_levels_five_percent_passed'])
    raise SystemExit(0 if passed else 1)


if __name__ == '__main__':
    main()
