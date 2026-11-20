"""Run the small T1 controls and audit the accepted hydraulic states in m³."""
import argparse
import copy
import gzip
import hashlib
import json
import math
import subprocess
import sys
import zipfile
from pathlib import Path
from tempfile import TemporaryFile

ROOT = Path(__file__).resolve().parents[2]
EX = ROOT / 'examples'
sys.path.insert(0, str(ROOT / 'tools'))
from t1_water_balance import Account, LIMITS, audit, split


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def write(path, value):
    with path.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write('\n')


def cases():
    base = read(EX / 'one-way-mc.json')
    joint = read(EX / 'stages-joint-water-cycle.json')
    for key in ('settings', 'boundary', 'joint_system', 'monitors', 'monitor_sets',
                'pressure_on_original_surface_only'):
        base['mechanics'][key] = joint[key]
    flow = base['flow']
    flow['initial_head'] = {'head_m': 2.}
    flow['materials'][0]['volume_source_per_s'] = 1e-6
    bottom = {'face_set': 'bottom', 'condition': {'kind': 'head', 'value': {'head_m': 2.1}}}
    flow['analysis']['periods'] = [dict(name='step', end_time_s=.25, boundaries=[bottom])]
    flow['settings'].update(initial_time_step_s=.125, maximum_time_step_s=.125)
    base['exterior_periods'] = [dict(name='step', reservoirs=[])]
    base['mechanics']['joint_system']['joints'][0]['water'] = {'kind': 'continuous'}
    yield 'continuous', base
    for name in ('soil', 'prescribed', 'table', 'staggered'):
        case = copy.deepcopy(base)
        system = case['mechanics']['joint_system']
        if name == 'staggered':
            case['schema'] = 'soil-fem-staggered-transfer/1'
        elif name == 'soil':
            system['structure_bodies'] = ['UpperClay']
            system['joints'][0]['water'] = dict(kind='soil', soil_side='a')
        elif name == 'prescribed':
            system['joints'][0]['water'] = dict(kind='prescribed', level_m=2.)
        else:
            system['joints'][0]['water'] = dict(kind='prescribed')
            case['joint_reservoir_heads'] = {'Seam': [[0., 2.], [.0625, 2.25], [.25, 1.875]]}
        yield name, case
    one = copy.deepcopy(base)
    one['mechanics']['joint_system']['joints'][0]['birth_policy'] = 'stress_free'
    periods = []
    names = ['hold', 'place', 'remove', 'replace']
    for i, name in enumerate(names):
        faces = ['xmin', 'xmax', 'ymin', 'ymax', 'top' if i % 2 else 'Seam:a']
        bcs = [dict(face_set='bottom', condition=dict(kind='head', value={'head_m': 2.}))]
        bcs.extend(dict(face_set=f, condition=dict(kind='outward_flux', flux_m_s=0.)) for f in faces)
        periods.append(dict(name=name, end_time_s=(i+1)*.125, boundaries=bcs))
    one['flow']['analysis']['periods'] = periods
    one['exterior_periods'] = [dict(name=n, reservoirs=[]) for n in names]
    birth = read(EX / 'one-way-lifecycle.json')['transitions'][0]
    birth['action']['initial_head'] = {'head_m': 2.}
    birth['water_source'] = {'reservoir_id': 'fill_water'}
    exc = dict(name='remove', before_period='remove', bodies=['UpperClay'],
               water_boundary_policy='explicit_period_boundaries',
               cut_traction_policy='explicit_period_reservoirs',
               water_destination={'reservoir_id': 'haul_water'})
    yield 'lifecycle', dict(schema='soil-fem-one-way-construction/1', one_way=one,
                            initially_inactive_bodies=['UpperClay'], excavations=[exc],
                            transitions=[dict(birth, name=n, before_period=n) for n in ('place', 'replace')])


def fixed_domain(mesh, report, profile):
    """A C-04 fixed mesh has constant reference storage, canceled in all differences."""
    rows = []
    accounts = {}
    sources = {m['id']: m.get('volume_source_per_s', 0.) for m in report['case']['materials']}
    for prev, frame in zip(report['frames'], report['frames'][1:]):
        state = frame['state']
        time = state['time_s']
        dt = time - prev['state']['time_s']
        before, after = (f['state']['cell_storage_offset_m3'] for f in (prev, frame))
        q = state['cell_outflow_m3_s']
        def add(kind, key, m0, m1, p, n):
            label = f'{kind}/{key}'
            if label not in accounts:
                accounts[label] = Account(m0)
            rows.append(dict(kind=kind, id=label, time_s=time,
                             **accounts[label].record(m0, m1, p, n, profile)))
        vols = [0.] * len(before)
        for point in frame['points']:
            vols[point['element_index']] += point['weight_m3']
        flows = []
        for i, (v, mat) in enumerate(zip(vols, mesh['material_ids'], strict=True)):
            source = -dt * v * sources[mat]
            flows.append(source)
            add('cell', i, before[i], after[i], *split([dt*x for x in q[i]]+[source]))
        for k, face in enumerate(report['faces']):
            exchange = [dt*q[i][j] for i,j in face['owners']]
            if len(exchange) == 2:
                add('face_pair', k, 0., 0., *split(exchange))
            else:
                flows.extend(exchange)
        add('global', 'domain', math.fsum(before), math.fsum(after), *split(flows))
    return dict(rows=rows, complete=report['complete'],
                passed=bool(report['complete'] and rows and all(r['passed'] for r in rows)))


def record(name, mesh, report, out, summary, fixed=False):
    primary = 'practical' if name in ('table', 'lifecycle') else 'verification'
    item = summary[name] = dict(primary=primary, frames=len(report['frames']), profiles={})
    for profile in LIMITS:
        result = (fixed_domain if fixed else audit)(mesh, report, profile)
        rows = result['rows']
        # Each profile has its own step and cumulative upper bounds.
        with gzip.open(out / f'{name}-{profile}.json.gz', 'wt', encoding='utf-8') as stream:
            json.dump(rows, stream, allow_nan=False)
        kinds = {}
        for kind in sorted({r['kind'] for r in rows if 'e' in r}):
            vs = [r for r in rows if r['kind'] == kind]
            worst = max(vs, key=lambda r: abs(r['e']))
            kinds[kind] = dict(rows=len(vs), failures=sum(not r['passed'] for r in vs),
                               worst_step=worst, max_abs_E=max(abs(r['E']) for r in vs),
                               max_U=max(r['U'] for r in vs))
        item['profiles'][profile] = dict(passed=result['passed'], kinds=kinds)
    if name == 'staggered':
        item['iterations'] = [f['staggered_iterations'] for f in report['frames'] if f['phase']=='transient']


def legacy_refusals(out, solver):
    old = ROOT / 'interim/2026-11-29-richards-wet-joint-a04'
    manifest = read(old / 'evidence-manifest.json')
    results = {}
    with TemporaryFile() as tmp:
        for part in manifest['parts']:
            data = (old / part['name']).read_bytes()
            if hashlib.sha256(data).hexdigest() != part['sha256']:
                raise ValueError(part['name'])
            tmp.write(data)
        with zipfile.ZipFile(tmp) as archive:
            for p in archive.namelist():
                if not p.startswith('runs/formal-refusals-01/') or not p.endswith('/case.json'):
                    continue
                folder = Path(p).parent
                d = out / folder.name
                d.mkdir()
                for file in ('mesh.json', 'case.json', 'ground.json'):
                    (d / file).write_bytes(archive.read((folder / file).as_posix()))
                cmd = 'solve-one-way-construction' if folder.name.startswith('construction-') else 'solve-one-way'
                run = subprocess.run([solver, cmd, d/'mesh.json', d/'ground.json', d/'case.json', d/'result.json'],
                                     capture_output=True, text=True)
                results[folder.name] = dict(code=run.returncode, diagnostic=(run.stdout+run.stderr).strip(),
                                            produced=(d/'result.json').exists())
            # Replay the earlier soil-only floor inconsistency from its archived input.
            d = out / 'floor-conflict'
            d.mkdir()
            for file in ('mesh.json', 'case.json', 'ground.json'):
                (d / file).write_bytes(archive.read('runs/formal-soil-02/saturated-n4/' + file))
            proc = subprocess.run([solver, 'solve-one-way-construction', d/'mesh.json', d/'ground.json',
                                   d/'case.json', d/'result.json'], capture_output=True, text=True)
            write(d/'observation.json', dict(expected_refused=True, code=proc.returncode,
                 diagnostic=(proc.stdout+proc.stderr).strip(), produced=(d/'result.json').exists()))
    if len(results) != 16:
        raise ValueError('incomplete A-04 refusal matrix')
    write(out/'refusal-results.json', results)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=['controls','a04','c04','refusals'])
    p.add_argument('output', type=Path)
    p.add_argument('--run', type=Path)
    p.add_argument('--solver', type=Path, default=ROOT/'target/release'/('soil-fem-3d.exe' if sys.platform=='win32' else 'soil-fem-3d'))
    a = p.parse_args()
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    summary = {}
    if a.mode == 'controls':
        mesh = EX/'joint-water-column-mesh.json'
        ground = read(EX/'ground-joint-water-mc.json')
        ground['initialization']['water_level_m'] = None
        ground['initialization']['geostatic']['hydrostatic_face_sets'] = []
        write(out/'ground.json', ground)
        for name, case in cases():
            case_path, report = out/f'{name}.json', out/f'{name}-report.json'
            write(case_path, case)
            command = 'solve-one-way-construction' if name=='lifecycle' else 'solve-one-way'
            subprocess.run([a.solver, command, mesh, out/'ground.json', case_path, report], check=True)
            record(name, read(mesh), read(report), out, summary)
    elif a.mode == 'refusals':
        legacy_refusals(out, a.solver)
        return
    else:
        patterns = ['VG-n*','saturated-n*'] if a.mode=='a04' else ['hydraulic-h*','saturated-controls/h*']
        for pattern in patterns:
            for folder in sorted(a.run.glob(pattern)):
                if folder.is_dir():
                    fixed = a.mode=='c04'
                    mesh = a.run/'hydraulic-mesh.json' if fixed else folder/'mesh.json'
                    record(folder.name, read(mesh), read(folder/'native.json'), out, summary, fixed)
        if len(summary) != 6:
            raise ValueError(f'incomplete {a.mode} matrix: {list(summary)}')
    write(out/'budget-summary.json', summary)


if __name__ == '__main__':
    main()
