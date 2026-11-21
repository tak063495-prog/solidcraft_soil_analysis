"""Additional diagnostics of the recorded indeterminate dry-MC in-plane series.

The original report and its half-increment control are read-only evidence. No
threshold, convergence decision, or failure/collapse interpretation is changed.
"""
import argparse
from pathlib import Path
import numpy as np
from plastic_spatial_convergence import read, write, sha, richardson


def diagnose(series_path, recorded_path, half_path, output):
    series, recorded, half = [read(p) for p in (series_path, recorded_path, half_path)]
    for report in (series, recorded, half):
        if (not report['complete'] or report['finite_width']
                or report.get('refinement_axis', 'in_plane') != 'in_plane'
                or [v['divisions'] for v in report['levels']] != [2,4,8]
                or any(v['depth_divisions'] != 1 for v in report['levels'])):
            raise ValueError('completed original one-layer extruded 2/4/8 in-plane series required')
    for a,b,h in zip(series['levels'],recorded['levels'],half['levels'],strict=True):
        for name in ('mesh.json','ground.json','case.json'):
            if a['files_sha256'][name] != b['files_sha256'][name]:
                raise ValueError('remeasurement input differs from original recorded series')
        for name in ('mesh.json','ground.json'):
            if a['files_sha256'][name] != h['files_sha256'][name]:
                raise ValueError('half control uses a different mesh or ground')
    # Verify the remeasurement rather than replacing archived values with it.
    for new,old in zip(series['pairs'],recorded['pairs'],strict=True):
        for a,b in zip(new['stages'],old['stages'],strict=True):
            if a['stage'] != b['stage']:
                raise ValueError('stage mismatch')
            for field,key in [('displacement','rms_difference_m'),('construction_displacement','rms_difference_m'),
                              ('stress','rms_difference_kpa')]:
                if not np.isclose(a[field][key],b[field][key],rtol=1e-10,atol=1e-14):
                    raise ValueError('remeasurement does not reproduce archived spatial differences')
    stress = [p['stages'][-1] for p in series['pairs']]
    original_stress = [p['stages'][-1]['stress']['rms_difference_kpa'] for p in recorded['pairs']]
    component_rms = np.asarray([p['stress_diagnostic']['tensor_weighted_component_rms_kpa'] for p in stress])
    initial, total, construction = [np.asarray([v['frames'][i][key] for v in recorded['levels']])
        for i,key in [(0,'monitor_uz_m'),(-1,'monitor_uz_m'),(-1,'construction_monitor_uz_m')]]
    if not np.allclose(initial+construction,total,rtol=0,atol=1e-15):
        raise ValueError('monitor decomposition does not close')
    control = np.asarray([v['frames'][-1]['construction_monitor_uz_m'] for v in half['levels']])-construction
    changes = np.diff(construction)
    report = dict(schema='soil-fem-plastic-spatial-diagnostic/1',
        sources_sha256={str(p):sha(p) for p in (series_path,recorded_path,half_path)},
        archived_values_reproduced=True, original_records_modified=False,
        stress=dict(recorded_difference_norms_kpa=original_stress,
            recorded_ratio=original_stress[1]/original_stress[0], guard_ratio=0.7,
            order=None, estimate=None, status='indeterminate',
            component_order=stress[0]['stress_diagnostic']['component_order'],
            tensor_weighted_component_rms_kpa=component_rms.tolist(),
            component_difference_ratios=[float(b/a) if a > 1e-12 else None for a,b in component_rms.T],
            pair_decompositions=[p['stress_diagnostic'] for p in stress]),
        monitor=dict(position_m=[.5,0.,.5], initial_uz_m=initial.tolist(), total_uz_m=total.tolist(),
            construction_uz_m=construction.tolist(),
            initial_signed_changes_m=np.diff(initial).tolist(), total_signed_changes_m=np.diff(total).tolist(),
            construction_signed_changes_m=changes.tolist(),
            half_control_changes_m=control.tolist(),
            half_control_pair_change_bounds_m=(abs(control[:-1])+abs(control[1:])).tolist(),
            original_richardson=richardson(construction),
            half_control_richardson=richardson(construction+control),
            status='indeterminate',
            interpretation='Signed changes of total and initial equilibrium subtract to give construction changes; cancellation is an observation, not a localization or uniqueness diagnosis.'),
        p6_complete=False, general_3d_spatial_convergence_demonstrated=False,
        fs=dict(status='not_attempted', nonconvergence_is_upper_bound=False),
        limitations=['Component/phase partitions diagnose the same full tensor norm and cannot replace its original gate.',
                     'Half control changes load and material increments together; it does not eliminate all numerical error.',
                     'No localization, uniqueness, collapse, or continuum-convergence conclusion follows.'])
    write(output,report)
    return report


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('series','recorded','half','output'):
        parser.add_argument(key,type=Path)
    args=parser.parse_args()
    diagnose(args.series,args.recorded,args.half,args.output)
