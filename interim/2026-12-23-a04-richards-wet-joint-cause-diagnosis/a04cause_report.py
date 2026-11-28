"""Read-only causal observations from new solves, never a coupled acceptance."""
import argparse
import math
from pathlib import Path
import numpy as np
from a04cause_binding import ROOT
from a04cause_driver import read, write, sha, exact
from a04cause_metrics import compare, convergence, rows, metric
HERE = Path(__file__).resolve().parent
SOIL_KEYS = ('cell_head_m','face_head_m','outflow_m3_s','storage_offset_m3','u','stress',
             'fluid_internal','water_body','exterior','equivalent_rhs','water_m3','mass_defect_m3')


def canonical(value):
    import json
    return json.dumps(exact(value), sort_keys=True, separators=(',', ':'), allow_nan=False)


def control_rows(path):
    result = []
    allowed = {'stage', 'progress', 'hydraulic_accounting', 'stage_index', 'factor', 'relative_force_residual', 'contact_partitions'}
    for row, fields in rows(path):
        if set(row)-allowed or not {'stage', 'progress'} <= set(row):
            raise ValueError('unexpected control coordinate; no physical clock may be invented')
        result.append(({k:row[k] for k in ('stage', 'progress', 'hydraulic_accounting') if k in row}, fields))
    return result


def soil_observation(folder):
    series = rows(folder/'native-fields.jsonl')
    mesh = read(folder/'mesh.json')
    maxima = dict(psi_identity=0., pressure_identity=0., bishop_identity=0., chi_p_identity=0.,
                  mass_defect_m3=0., force_decomposition_kn=0., cell_trace_difference_m=0.)
    release = []
    negative = 0
    trace_location = None
    source_hash = sha(folder/'native.json')
    for row, f in series:
        for k in ('hydraulic', 'mechanical'):
            psi = np.asarray(f[k+'_pressure_head_m'])
            p = np.asarray(f[k+'_physical_pore_pressure_kpa'])
            chi = np.asarray(f[k+'_bishop_chi'])
            se = np.asarray(f[k+'_effective_saturation'])
            cp = np.asarray(f[k+'_effective_pore_stress_kpa'])
            maxima['pressure_identity'] = max(maxima['pressure_identity'], metric(p, 10*psi)['maximum_absolute'])
            maxima['bishop_identity'] = max(maxima['bishop_identity'], metric(chi, se)['maximum_absolute'])
            maxima['chi_p_identity'] = max(maxima['chi_p_identity'], metric(cp, chi*p)['maximum_absolute'])
        h = np.repeat(f['cell_head_m'], len(f['hydraulic_positions'])//len(f['cell_ids']))
        z = np.asarray(f['hydraulic_positions'])[:, 2]
        maxima['psi_identity'] = max(maxima['psi_identity'], metric(f['hydraulic_pressure_head_m'], h-z)['maximum_absolute'])
        negative += int(np.count_nonzero(np.asarray(f['hydraulic_pressure_head_m']) < 0))
        maxima['mass_defect_m3'] = max(maxima['mass_defect_m3'], abs(f['mass_defect_m3']))
        maxima['force_decomposition_kn'] = max(maxima['force_decomposition_kn'], metric(
            np.asarray(f['effective_internal'])+f['fluid_internal']-np.asarray(f['physical_load'])-f['holding_load'], f['physical_reaction'])['maximum_absolute'])
        lookup = {tuple(face): i for i, face in enumerate(f['face_keys'])}
        for local, e in enumerate(f['cell_ids']):
            cell = mesh['tetrahedra'][e]
            for opposite in range(4):
                key = tuple(sorted(n for i, n in enumerate(cell) if i != opposite))
                d = abs(f['cell_head_m'][local]-f['face_head_m'][lookup[key]])
                if d > maxima['cell_trace_difference_m']:
                    maxima['cell_trace_difference_m'] = d
                    trace_location = dict(row=row, element=e, face=key,
                                          cell_H_m=f['cell_head_m'][local], trace_H_m=f['face_head_m'][lookup[key]])
        if row['phase'] == 'excavation_release':
            release.append((row, f))
    hydraulic = ('cell_head_m','face_head_m','outflow_m3_s','storage_offset_m3','water_m3','removed_water_m3','net_inflow_m3')
    frozen = all(all(np.array_equal(f[k], release[0][1][k]) for k in hydraulic) for _, f in release)
    if sha(folder/'native.json') != source_hash:
        raise ValueError('observer modified native run')
    return dict(maxima=maxima, negative_pressure_point_states=negative, release_rows=[r for r, _ in release],
                release_hydraulics_frozen_exact=frozen, native_run_sha256_unchanged=source_hash,
                maximum_cell_trace_difference_location=trace_location,
                interface_mapping_created=False, cell_vs_trace_is_not_a_JOINT_head_error=True)


def report(run):
    driver = read(run/'driver-result.json')
    levels = driver['steps']
    soil = {}; dt = {}; quad = {}; reproducibility = {}; observations = {}
    for label in ('VG', 'saturated'):
        for n in levels:
            folder = run/'soil'/(label+'-n'+str(n))
            soil[folder.name] = compare(rows(folder/'native-fields.jsonl'), rows(folder/'independent-fields.jsonl'))
            observations[folder.name] = soil_observation(folder)
        for kind in ('native','independent'):
            series = [rows(run/'soil'/(label+'-n'+str(n))/(kind+'-fields.jsonl')) for n in levels]
            grid = [(r,f) for r,f in series[0] if r['phase'] == 'transient']
            dt[label+'/'+kind] = convergence(series, grid, SOIL_KEYS)
        base = run/'soil'/(label+'-n4')
        if (run/(label+'-repeat')).exists():
            native = read(base/'native.json')
            reproducibility[label] = {variant: canonical(read(run/(label+'-'+variant)/'native.json')) == canonical(native)
                                      for variant in ('repeat','reverse-keys')}
            reproducibility[label]['independent_reverse_keys'] = compare(rows(base/'independent-fields.jsonl'),
                rows(run/(label+'-reverse-keys')/'independent-fields.jsonl'))
            series = [rows(run/'soil'/(label+'-n8')/'independent-fields.jsonl')]+[
                rows(run/(label+'-quadrature'+str(n))/'independent-fields.jsonl') for n in (5,7)]
            quad[label] = dict(orders=[3,5,7], same_cell_face_nodal_quantities_only=True,
                              GP_arrays_not_remapped=True, full_equilibrium_solved_from_inputs=True,
                              pairs=[compare(series[i], series[i+1], SOIL_KEYS) for i in range(2)])
    controls = {}
    for kind in ('dry','hydrostatic'):
        if (run/kind).exists():
            data = read(run/kind/'measurements.json')
            controls[kind] = dict(fixed_order_comparison=compare(control_rows(run/kind/'case/native-fields.jsonl'), control_rows(run/kind/'case/independent-fields.jsonl')),
                legacy_comparison={k:v for k,v in data['comparison'].items() if k != 'all_points'},
                hydraulic_flow_m3_s=data['hydraulic_flow_m3_s'], hydraulic_storage_m3=data['hydraulic_storage_m3'],
                mass_balance=data['mass_balance'], physical_time_s=data['physical_time_s'],
                restarts=data['restarts'], tampering=data['tampering'], state_spoofing=data['state_spoofing'])
    gap = None
    if (run/'boundary-gap').exists():
        native = read(run/'boundary-gap/native.json')
        if 'declared reservoir and hydraulic head boundary are inconsistent' not in (run/'boundary-gap-independent.log').read_text(encoding='utf-8'):
            raise ValueError('boundary reference failed for an unexpected reason')
        gap = dict(expected_refused=True, actual_production_refused=False, actual_production_complete=native['complete'],
                   excluded_from_limited_acceptance=True, numerical_or_scope_acceptance=False,
                   reason='production validates exterior schedule names, not same-face H/level equality; independent constructor refuses before solve',
                   independent_process=read(run/'boundary-gap-independent.process.json'),
                   pairs=[dict(period='p2', head_m=1.8, reservoir_m=2.2),dict(period='p3', head_m=2.1, reservoir_m=2.2)])
    primary = {}
    for folder in ('2026-11-29-richards-wet-joint-a04','2026-11-22-wet-joint-transient-water',
                   '2026-12-19-c04-transient-water-cause-diagnosis','2026-12-20-a01-joint-combo-cause-diagnosis',
                   '2026-12-21-a02-joint-birth-cause-diagnosis','2026-11-26-wet-soil-birth-darcy-a03',
                   '2026-12-22-a03-wet-soil-birth-cause-diagnosis'):
        primary[folder] = {p.name:sha(p) for p in sorted((ROOT/'interim'/folder).glob('*.json'))}
    main = read(run/'soil/measurements.json')
    checkpoints = [r for v in main['cases'].values() for r in v.get('restarts',[])]+[r for v in controls.values() for r in v['restarts']]
    mutations = [r for v in main['cases'].values() for r in v.get('tampering',[])]+[r for v in controls.values() for r in v['tampering']]+[r for v in controls.values() for r in v['state_spoofing'] if r['outer_hash_recomputed']]
    return dict(schema='a04-cause-measurements/1', diagnostic_only=True, passed=False, a04_complete=False,
                all_field_numerical_acceptance=False, source_binding_sha256=sha(run/'source-binding.json'),
                postprocessing_source_sha256={p.name:sha(p) for p in (HERE/'a04cause_report.py', HERE/'a04cause_metrics.py')},
                soil=soil, soil_definitions_clock_force=observations, physical_dt=dt,
                independent_quadrature=quad, reproducibility=reproducibility, controls=controls,
                original_numpy_soil_comparisons_retained_at='formal/soil/measurements.json in sealed archive',
                scope_refusals=read(run/'refusals/refusals.json') if (run/'refusals').exists() else None,
                same_surface_boundary_gap=gap,
                checkpoints=dict(real_files=len(checkpoints), suffix_states_exact=sum(r['suffix_states_exact'] for r in checkpoints),
                                 rehashed_content_refusals=len(mutations), not_fresh_solves=True),
                primary_record_sha256=primary,
                interface=dict(head_mapping=None, hydraulic_local_law=None, physical_time_generation=None,
                               storage_m3=None, flow_m3_s=None, mass_balance=None, coupled_dt_acceptance=None,
                               hydraulic_mechanical_coupled_reassembly=None, bidirectional_controller=False),
                diagnosis='contract gaps block coupled numerical attribution; fixed soil/control agreement is limited evidence',
                single_cause_identified=False, cloud_diagnostics_remaining=True, pc_executed=False)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    write(a.output,report(a.run.resolve()))


if __name__ == '__main__':
    main()
