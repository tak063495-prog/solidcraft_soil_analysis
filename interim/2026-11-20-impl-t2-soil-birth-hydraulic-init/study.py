"""Run the small placement controls through the production CLI and audit the outputs."""
import json
import argparse
import copy
import math
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools'))
from t1_water_balance import LIMITS, audit
from verify_joint_ext_combos import wet_fixture
from importlib.util import module_from_spec, spec_from_file_location

t1_spec = spec_from_file_location('t1_controls', ROOT/'interim/2026-11-20-impl-t1-wet-joint-richards/study.py')
t1 = module_from_spec(t1_spec)
t1_spec.loader.exec_module(t1)


def prescribed(kind):
    mesh, ground, case = wet_fixture()
    case['stages'] = case['stages'][:3]
    case['verification_reject_once_stages'] = []
    case['settings']['kinematics'] = 'side_rollers'
    case.pop('boundary', None)
    birth = case['stages'][0]
    init = dict(kind=kind, placement_water=dict(material_id=2, water_source=dict(reservoir_id='carried_with_fill')))
    birth['hydraulic_initialization'] = init
    if kind == 'post_darcy':
        init.update(preceding_stage='flow', initial_head=dict(head_m=1.5, gradient=[0.,0.,.2], reference_point_m=[0.,0.,.5]))
        flow = t1.read(ROOT/'examples/stages-joint-water-flow.json')['stages'][0]
        flow['flow']['heads'] = [dict(face_set='bottom', head_m=1.6), dict(face_set='Seam:a',head_m=1.5)]
        flow['exterior_water_level_m'] = 1.5
        case['stages'].insert(0, flow)
    elif kind == 'two_region':
        init.update(level_m=.75, mixture=dict(porosity=.7/1.7, upper_saturation=.25, gravity_m_s2=10.,
                    particle_density_kg_m3=2700., water_density_kg_m3=1000.))
        ground['initialization']['water_level_m'] = .75
    else:
        template = t1.read(ROOT/'examples/ground-joint-water-mcc.json')
        seed = template['initialization']['geostatic']['seed']['states']['UpperClay']
        e = seed['void_ratio']
        ground['materials'][1] = template['materials'][1]
        ground['initialization']['geostatic']['seed']['states']['UpperClay'] = seed
        birth['states']['UpperClay'] = seed
        init.update(level_m=1.5, formation=dict(history='Synthetic drained isotropic calibration, 200 kPa to 100 kPa',
                    hydraulic_porosity=e/(1.+e), pc_over_mean_effective=2.))
    return mesh, ground, case


def construction():
    case = dict(t1.cases())['lifecycle']
    ground = t1.read(ROOT/'examples/ground-joint-water-mc.json')
    ground['initialization']['water_level_m'] = None
    ground['initialization']['geostatic']['hydrostatic_face_sets'] = []
    material = copy.deepcopy(ground['materials'][1])
    material.update(id=3, name='Replacement soil')
    ground['materials'].append(material)
    ground['construction_material_ids'] = [3]
    hydraulic = copy.deepcopy(case['one_way']['flow']['materials'][1])
    hydraulic['id'] = 3
    case['one_way']['flow']['materials'].append(hydraulic)
    event = copy.deepcopy(case['transitions'][1])
    case['transitions'][1] = event
    event['action']['material_id'] = 3
    event['action']['initial_head']['head_m'] = 2.25
    event['action']['states']['UpperClay']['effective_stress_compression_positive_kpa'] = [1.,1.,1.,0.,0.,0.]
    return t1.read(ROOT/'examples/joint-water-column-mesh.json'), ground, case


def check_restarts(solver, paths, folder):
    def read_numerical_report(path):
        def omit_wall_time(fields):
            fields.pop('elapsed_seconds', None)
            return fields
        return json.loads(path.read_bytes(), object_hook=omit_wall_time)

    reference = read_numerical_report(folder / 'result.json')
    results = []
    for count in (6, 8, 12, 14):
        checkpoint = folder / f'checkpoint-{count}.json'
        resumed = folder / f'resumed-{count}.json'
        for command, source, output in (
            ('checkpoint-one-way-construction', str(count), checkpoint),
            ('resume-one-way-construction', checkpoint, resumed),
        ):
            proc = subprocess.run([solver, command, *paths, source, output], capture_output=True,
                                  text=True, timeout=3600)
            (folder / f'{command}-{count}.log').write_text(proc.stdout + proc.stderr, encoding='utf-8')
            proc.check_returncode()
        restored = read_numerical_report(resumed)
        identical = restored == reference
        if not restored['complete'] or not identical:
            raise ValueError(f'restart disagrees with uninterrupted placement at frame {count}')
        results.append(dict(accepted_frames=count, identical=identical))
    return results


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('output', type=Path)
    p.add_argument('--solver', type=Path, default=ROOT/'target/release'/('soil-fem-3d.exe' if sys.platform=='win32' else 'soil-fem-3d'))
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    summary = {}
    for kind in ('post_darcy','two_region','saturated_mcc','lifecycle'):
        folder = args.output/kind
        folder.mkdir()
        mesh, ground, case = construction() if kind=='lifecycle' else prescribed(kind)
        paths = [folder/(n+'.json') for n in ('mesh','ground','case')]
        for path, value in zip(paths, (mesh,ground,case), strict=True):
            t1.write(path, value)
        output = folder/'result.json'
        command = 'solve-one-way-construction' if kind=='lifecycle' else 'solve-stages'
        proc = subprocess.run([args.solver, command, *paths, output], capture_output=True, text=True, timeout=3600)
        (folder/'solve.log').write_text(proc.stdout+proc.stderr, encoding='utf-8')
        item = summary[kind] = dict(exit_code=proc.returncode, complete=False)
        if output.exists():
            run = t1.read(output)
            item.update(complete=run['complete'], failure=run.get('failure_reason'))
            if kind=='lifecycle':
                budget = {}
                t1.record(kind, mesh, run, folder, budget)
                item['budget'] = budget[kind]
                item['checkpoint_restarts'] = check_restarts(args.solver, paths, folder)
                item['events'] = [dict(name=r['name'], time_s=r['time_s'], generations=r['placed_generations'],
                                      added_m3=r['water']['water_added_m3']) for r in run['transitions']]
                item['excavations'] = [dict(name=r['name'], time_s=r['time_s'], generations=r['removed_generations'],
                                      removed_m3=r['removal']['water_removed_m3']) for r in run['excavations']]
            else:
                item['budgets'] = {profile:audit(mesh, run, profile) for profile in LIMITS}
                placement = next((r for r in run['transition_audits'] if r.get('placed_water')), None)
                if placement:
                    item['physical_quantities'] = {key:math.fsum(c[key] for c in placement['placed_water'])
                        for key in ('upper_volume_m3','lower_volume_m3','upper_water_m3','lower_water_m3','upper_weight_kn','lower_weight_kn')}
                    item['zero_release_relative_force_residual'] = placement['zero_release_relative_force_residual']
        else:
            item['failure'] = proc.stderr.strip()
    t1.write(args.output/'results.json', summary)


if __name__ == '__main__':
    main()
