"""Remeasure saved observations and prior evidence, without running a solver."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from b01cause_geometry import dot, frame, metric, norm, rigid, shape

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PRIOR = ROOT/'interim/2026-12-03-tri6-curved-b01'


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def write(path, obj):
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False)+'\n',
                          encoding='utf-8', newline='\n')


def prior_assessment():
    saved = read(PRIOR/'curved-comparison.json')
    roof = [{k: r[k] for k in ('n', 'native_signed_probe', 'draft_signed_probe',
                               'native_minus_draft', 'literature_signed_probe',
                               'native_relative_literature_error')} for r in saved['roof']]
    attached = [dict(n=r['n'], uz_m=r['crown_y0_uz_m']) for r in saved['attached']]
    d1 = attached[1]['uz_m']-attached[0]['uz_m']
    d2 = attached[2]['uz_m']-attached[1]['uz_m']
    return dict(saved_reanalysis_only=True, fresh_benchmark_solves=0,
                old_source_sha256=hashlib.sha256((PRIOR/'curved-comparison.json').read_bytes()).hexdigest(),
                roof=roof, roof_observed_order=None, roof_levels_insufficient_for_order=True,
                native_minus_q4_sign_reversal=roof[0]['native_minus_draft']*roof[1]['native_minus_draft'] < 0,
                attached=attached, attached_signed_differences_m=[d1, d2],
                attached_diagnostic_order=math.log2(abs(d1/d2)),
                attached_scope='circumferential refinement only, fixed soil and one axial division, faceted planar elements',
                reference_status='finite Q4 analytic-director discretization; not exact and not the same element as faceted Tri6',
                legacy_references_unchanged=dict(roof=-.3024, cylinder=-1.8248e-5, open_top_hemisphere=.094),
                source_confirmed_variant_values=dict(
                    Ko2017_roof=-.3024, Ko2017_cylinder=-1.8248e-5,
                    Ko2017_open_top_hemisphere=.093, Ko2017_full_hemisphere=.0924,
                    Abaqus2025_LE3_target_mm=185., Abaqus2025_cylinder=1.825e-5,
                    Abaqus2025_roof_target_mm=91.2),
                reference_definitions_interchangeable=False,
                source_warning='KC2023 sections 2/4 distinguish KL/MR and point-load singularities; no reference/definition/tolerance replacement adopted',
                original_numerical_acceptance=None, curved_benchmarks_implemented=False)


def summarize(formal):
    report = read(formal/'internal-observations.json')
    requests = {r['id']: r for r in read(formal/'requests.json')}
    raw = {r['id']: r for r in read(formal/'probe.stdout.json')}
    observations = {r['id']: r for r in report['observations']}
    main = []
    for name in ('plane', 'inplane-bow', 'outplane-bulge', 'cylinder', 'sphere'):
        a = observations[name]; k = np.array(raw[name]['stiffness']); x = np.array(requests[name]['points'])
        eigen, vectors = np.linalg.eigh((k+k.T)*.5)
        z = vectors[:, eigen <= 1e-10*eigen[-1]]; r = rigid(x)
        outside = [norm(r[:, i]-z@(z.T@r[:, i]))/norm(r[:, i]) for i in range(6)]
        energies = []
        for field in a['fields']:
            value = field['recovery']
            energies.append(dict(name=field['name'], **{key: value[key] for key in (
                'membrane_energy_kn_m', 'bending_energy_kn_m',
                'transverse_shear_energy_kn_m', 'drilling_penalty_energy_kn_m')}))
        # Pull back the OLD constant-plane membrane strain to the actual
        # quadratic tangents. It lacks the normal-displacement contribution.
        isotropic = next(f['recovery'] for f in a['fields'] if f['name'] == 'isotropic-membrane')
        pullback_errors = []
        axes = np.array(raw[name]['axes'])
        for sample in isotropic['samples']:
            l = sample['barycentric']; f = frame(x, l[1], l[2])
            ex, ey, gamma = sample['membrane_strain']
            old = np.array([[ex, gamma*.5], [gamma*.5, ey]])
            projected = f['tangent']@axes[:2].T
            pullback_errors.append((projected@old@projected.T-.001*f['metric']).tolist())
        main.append(dict(id=name, production_guard=a['production_guard'],
                         inherited_zero_count=a['spectrum']['modes_below_inherited_threshold'],
                         physical_rigid_residual=a['spectrum']['physical_six_rigid_normalized_residual'],
                         physical_rigid_fraction_outside_algebraic_nullspace=outside,
                         exact_quadratic_rigid_membrane_roundoff=a['geometry']['exact_quadratic_rigid_covariant_membrane_max'],
                         geometry_normal_derivative=a['geometry']['maximum_normal_derivative'],
                         geometry_second_form=a['geometry']['maximum_second_form'],
                         internal_K_exactly_equal_to_plane_same_corners=np.array_equal(k, np.array(raw['plane']['stiffness'])) if name in ('plane','inplane-bow','outplane-bulge') else None,
                         isotropic_membrane_pullback_error_l2=norm(pullback_errors), energy_partitions=energies))
    order = [a['stiffness_after_physical_id_permutation']['relative_l2']
             for a in report['observations'] if 'stiffness_after_physical_id_permutation' in a]
    boundary = [dict(id=a['id'], accepted=a['production_guard']['accepted'],
                     cutoff_m=a['geometry']['guard_cutoff_m'],
                     actual_maximum_midpoint_distance_m=max(a['geometry']['midpoint_distances_m']))
                for a in report['observations'] if a['id'].startswith('boundary-')]
    sensitivity = []
    for a in report['observations']:
        if a['id'].startswith('bulge-t'):
            sensitivity.append(dict(id=a['id'], zero_count=a['spectrum']['modes_below_inherited_threshold'],
                                    rotation_x_energies=a['fields'][3]['recovery']))
    mesh_summary = []
    for a in report['meshes']:
        mesh_summary.append({key: a[key] for key in ('kind','n','nodes','elements','dofs','spectrum',
                           'element_stiffness_exactly_equals_same_corner_faceted_control')})
        if 'global_equilibrium_solves' in a:
            mesh_summary[-1]['global_solves'] = [{k: v for k, v in s.items()
                if k not in ('solution','prescribed_rigid','residual')} for s in a['global_equilibrium_solves']]
    return dict(schema='b01cause-measurements/1', passed=False, b01_complete=False, diagnostic_only=True,
                curved_element_acceptance=None, production_scope_expanded=False,
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                formal_observation_sha256=hashlib.sha256((formal/'internal-observations.json').read_bytes()).hexdigest(),
                internal_requests=len(requests), internal_meshes=6, global_equilibrium_solves=9,
                maximum_new_elements=8, maximum_new_structure_dofs=150,
                production_accepted_requests=sum(a['production_guard']['accepted'] for a in report['observations']),
                all_accepted_transplants_exactly_equal_production=all(
                    a['production_stiffness_exact_when_accepted'] is True for a in report['observations'] if a['production_guard']['accepted']),
                main=main, meshes=mesh_summary, boundary=boundary,
                maximum_node_order_K_relative_difference=max(order), thickness_drilling_sensitivity=sensitivity,
                saved_reference_assessment=prior_assessment(),
                conclusions=dict(guard='full 3D distance to arithmetic edge midpoint; not just a coplanarity test',
                    implementation='midpoint coordinates do not enter stiffness after guard; fixed affine geometry cannot represent the quadratic surface',
                    rigid_modes='six algebraic near-zero eigenvalues do not certify actual-coordinate rigid rotations',
                    global_equilibrium='zero-load rigid boundary tests solve assembled systems; internal force balance can coexist with erroneous artificial deformation/reactions',
                    drilling='spurious membrane/shear energy dominates selected examples; drilling is not the sole cause; no curved drilling law was tested',
                    reference='variant/model/point-load target definitions differ; finite saved series do not establish a curved Tri6 limit',
                    unique_cause_fully_identified=False),
                cloud_diagnostics_remaining=True,
                cloud_remaining=['small distorted-patch and frame-rotation cross controls',
                                 'compatible curved manufactured-load specification and independent element (requires design first)',
                                 'separate director-interpolation versus curvature versus MITC interventions after specification'],
                specification_decisions_needed=True, pc_executed=False, windows_ci_after_push='未確定')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--formal', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); write(a.output, summarize(a.formal))
