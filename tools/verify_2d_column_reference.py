"""Actual pinned Excel MC versus Rust 3D and an independent closed-form column.

The two programs use different meshes/quadrature and displacement origins.
Compare every node and GP, with a piecewise quadratic/linear 2D field sampled
at the actual 3D coordinates. No constitutive parameter is fitted to results.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tempfile
import numpy as np


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))


def maximum(a, b):
    return float(np.max(np.abs(np.asarray(a)-np.asarray(b))))


def serialized_identity(path, key):
    # Hash the actual serde-rendered object, removing formatting whitespace
    # outside strings only. Python's float spelling differs from serde's.
    text = path.read_text(encoding='utf-8-sig')
    decoder = json.JSONDecoder()
    index = 1
    while True:
        while text[index].isspace() or text[index] == ',':
            index += 1
        name, index = decoder.raw_decode(text, index)
        while text[index].isspace() or text[index] == ':':
            index += 1
        start = index
        _, index = decoder.raw_decode(text, index)
        if name == key:
            compact = re.sub(r'("(?:[^"\\]|\\.)*")|\s+', lambda m: m[1] or '', text[start:index])
            return hashlib.sha256(compact.encode('utf-8')).hexdigest()
        assert text[index:].lstrip()[0] != '}', key


def input_subset(expected, actual):
    # Deserialization adds declared defaults; every supplied field must still
    # match, including every list length and entry.
    if isinstance(expected, dict):
        for key, value in expected.items():
            input_subset(value, actual[key])
    elif isinstance(expected, list):
        assert len(expected) == len(actual)
        for a, b in zip(expected, actual):
            input_subset(a, b)
    else:
        assert expected == actual, (expected, actual)


def theory(z, stage):
    """Compression-positive stress/plastic strain; z-up displacement.

    With psi=0 and equal constrained transverse strains, plastic volume is
    zero. Thus ev=p/K even after yielding. The MC compression edge gives
    sv=A*p+B, A=3(1+sin(phi))/(3-sin(phi)), B=4c*cos(phi)/(3-sin(phi)).
    Gravity is elastic; unloading by P is elastic with modulus M.
    """
    z = np.asarray(z)
    e, nu, gamma, height, pressure, c = 1000., .25, 20., .2, 12., 2.
    sine, cosine = np.sin(np.deg2rad(10.)), np.cos(np.deg2rad(10.))
    bulk = e/(3*(1-2*nu))
    constrained = e*(1-nu)/((1+nu)*(1-2*nu))
    a, b = 3*(1+sine)/(3-sine), 4*c*cosine/(3-sine)
    sv0 = gamma*(height-z)
    ug = -gamma*(height*z-z*z/2)/constrained
    ul = -((pressure+gamma*height-b)*z-gamma*z*z/2)/(a*bulk)
    sv = sv0 if stage == 0 else sv0+pressure
    sh = nu/(1-nu)*sv0 if stage == 0 else (3*(sv-b)/a-sv)/2
    plastic = np.zeros((len(z), 6))
    if stage:
        epv = (sv-b)/(a*bulk)-(sv-2*nu*sh)/e
        plastic[:, 0] = plastic[:, 1] = -.5*epv
        plastic[:, 2] = epv
    if stage == 2:
        sv = sv-pressure
        sh = sh-nu*pressure/(1-nu)
    stress = np.column_stack((sh, sh, sv, np.zeros((len(z), 3))))
    u = (ug, ul, ul+pressure*z/constrained)[stage]
    return u, stress, plastic


def excel(base, ny, bound):
    directory = base/f'ny{ny}/excel-v3'
    data = read(directory/'result.json')
    fixture = read(base/f'ny{ny}/fixture-2d.json')
    assert data['status'] == 'PASS' and data['application_code_unchanged']
    assert not read(directory/'code-comparison.json')['program_changes']
    root = Path(__file__).resolve().parents[1]
    assert sha(root/'scripts/Invoke-2DReference.ps1') == data['runner_sha256'].lower()
    before, after = [read(directory/f'code-{name}.json') for name in ('before', 'after')]
    for name in set(before) | set(after):
        if name in before and name in after:
            assert before[name] == after[name], name
        else:
            entry = (before if name in before else after)[name]
            assert entry['lines'] == 0 and entry['component_type'] == 100, name
    assert sha(base/f'ny{ny}/fixture-2d.json') == data['fixture_sha256'].lower()
    assert sha(base/'source/2DSoilFEM_ui_20261024.xlsm') == data['source_workbook_sha256'].lower()
    assert sha(base/'source/soilfem_native_x64.dll') == data['source_dll_sha256'].lower()
    for table, name in (('節点データ', 'nodes'), ('要素データ', 'elements'), ('材料データ', 'materials'), ('載荷', 'loads'), ('ステージ', 'stages')):
        assert data['tables'][table][1:] == fixture[name], table
    actual_settings = {row[4]: row[2] for row in data['tables']['設定'] if row[4]}
    for key, value in fixture['settings'].items():
        assert actual_settings[key] == value or (value == '' and actual_settings[key] is None), key
    out = []
    for stage, name in enumerate(('stg01_GRAVITY', 'stg02_APPLY', 'stg03_UNLOAD')):
        folder = directory/'case_out'/name
        meta = dict(line.split('=', 1) for line in (folder/'meta.txt').read_text(encoding='utf-8-sig').splitlines() if '=' in line)
        assert meta['STATUS'] == 'SUCCESS' and meta['GRAVITY'] == '1'
        node, gauss, cumulative = [rows(folder/f'{n}.csv') for n in ('node', 'gauss', 'plastic_cumulative')]
        assert len(node) == len(fixture['nodes']) and len(gauss) == len(cumulative) == 4*ny
        assert {int(r['node']) for r in node} == set(range(1, len(node)+1))
        for r in node:
            expected = fixture['nodes'][int(r['node'])-1]
            assert maximum([float(r['x']), float(r['y'])], expected[1:3]) < 6e-11
            assert [int(r['condx']), int(r['condy'])] == expected[3:5]
        ids = {(int(r['elem']), int(r['gp'])) for r in gauss}
        assert ids == {(e, g) for e in range(1, ny+1) for g in range(4)}
        assert ids == {(int(r['elem']), int(r['gp'])) for r in cumulative}
        cumulative_map = {(int(r['elem']), int(r['gp'])): r for r in cumulative}
        z = np.array([float(r['y']) for r in node])
        u = np.array([float(r['uy']) for r in node])
        gp_z = np.array([((int(r['elem'])-1)+(1+(-1 if int(r['gp']) < 2 else 1)*.5773502692)/2)*.2/ny for r in gauss])
        stress = np.array([[float(r[k]) for k in ('sx', 'sz', 'sy')]+[0., 0., float(r['txy'])] for r in gauss])
        plastic = []
        for r in gauss:
            p = cumulative_map[(int(r['elem']), int(r['gp']))]
            assert p['coverage'] == 'FROM_ZERO_STATE'
            plastic.append([float(p[k]) for k in ('ex_pl_cumulative', 'ez_pl_cumulative', 'ey_pl_cumulative')]+[0., 0., float(p['gxy_pl_cumulative'])])
        plastic = np.asarray(plastic)
        ut, _, _ = theory(z, stage)
        _, st, pt = theory(gp_z, stage)
        error = dict(displacement_maximum_m=maximum(u, ut), stress_maximum_kpa=maximum(stress, st),
                     plastic_strain_maximum=maximum(plastic, pt))
        assert error['displacement_maximum_m'] < bound['excel_u_m']
        assert error['stress_maximum_kpa'] < bound['excel_stress_kpa']
        assert error['plastic_strain_maximum'] < bound['excel_plastic_strain']
        assert max(abs(float(r['ux'])) for r in node) < 1e-11
        # Reaction totals are actual workbook output, independent of the CSV
        # field reconstruction. The retained selfweight is 0.8 kN.
        summary = data['tables']['ステージ結果'][stage+1]
        assert summary[3] == 'PASS' and abs(summary[7]-(.8 if stage != 1 else 3.2)) < 1e-8
        full_nodes = [r for r in data['tables']['結果節点'][1:] if int(r[0]) == stage+1]
        assert len(full_nodes) == len(node)
        assert abs(sum(float(r[9]) for r in full_nodes)-summary[7]) < 1e-8
        node_map = {int(r[1]): r for r in full_nodes}
        for r in node:
            source = node_map[int(r['node'])]
            assert maximum(source[2:6], [float(r[k]) for k in ('x', 'y', 'ux', 'uy')]) < 6e-13
        full_elements = {int(r[1]): r for r in data['tables']['結果要素'][1:] if int(r[0]) == stage+1}
        assert len(full_elements) == ny
        for r in gauss:
            source = full_elements[int(r['elem'])]
            offset = 4+17*int(r['gp'])
            assert maximum([source[offset+i] for i in (0, 1, 2, 16)],
                           [float(r[k]) for k in ('sx', 'sy', 'txy', 'sz')]) < 6e-12
        top = u[z == .2]
        out.append(dict(z=z, u=u, gp_z=gp_z, stress=stress, plastic=plastic,
                        metrics=dict(stage=name, **error, reaction_y_kn=summary[7], top_absolute_m=float(top.mean()))))
    assert np.array_equal(out[1]['plastic'], out[2]['plastic']), 'unloading changed cumulative plastic history'
    return out


def sample_2d(field, z, ny, gp=False):
    """Interpolate unsmoothed values inside each source Q8 element."""
    coordinate = field['gp_z' if gp else 'z']
    values = np.column_stack((field['stress'], field['plastic'])) if gp else field['u']
    samples = []
    for target in z:
        element = min(int(target/(.2/ny)), ny-1)
        lo, hi = element*.2/ny, (element+1)*.2/ny
        mask = (coordinate >= lo-1e-12) & (coordinate <= hi+1e-12)
        abscissa = np.unique(coordinate[mask])
        assert len(abscissa) == (2 if gp else 3)
        ordinates = np.asarray([values[np.isclose(coordinate, t, atol=1e-12, rtol=0)].mean(axis=0) for t in abscissa])
        sample = np.zeros_like(ordinates[0])
        for i, x in enumerate(abscissa):
            weight = np.prod([(target-y)/(x-y) for j, y in enumerate(abscissa) if j != i])
            sample += weight*ordinates[i]
        samples.append(sample)
    return np.asarray(samples)


def native(base, folder, excel_fields, bound):
    where = base/folder if folder else base
    names = ('mesh.json', 'stages.json', 'result.json') if folder else ('mesh-3d.json', 'stages-3d.json', 'result-3d.json')
    mesh, case, result = [read(where/n) for n in names]
    provenance = result['provenance']
    solver = Path(__file__).resolve().parents[1]/'target/release/soil-fem-3d.exe'
    assert provenance['solver_binary_sha256'] == sha(solver)
    for key, data in (('case', case), ('ground', read(base/'ground-3d.json'))):
        input_subset(data, result[key])
        assert serialized_identity(where/names[2], key) == provenance[key+'_sha256'], key
    with tempfile.TemporaryDirectory(prefix='soil-fem-column-identity-') as temporary:
        identity_path = Path(temporary)/'identity.json'
        subprocess.run([str(solver), 'mesh-identity', str(where/names[0]), str(identity_path)], check=True, capture_output=True)
        identity = read(identity_path)
    assert identity['mesh_sha256'] == provenance['mesh_sha256']
    assert identity['integration_points'] == result['integration_points']
    assert result['complete'] and not result['failure_reason'] and len(result['frames']) == 3
    assert all(a['accepted'] for a in result['attempts'])
    xyz = np.asarray(mesh['nodes'])
    gp_z = np.asarray([r['position_m'][2] for r in result['integration_points']])
    assert len(gp_z) == 27*len(mesh['tetrahedra10'])
    assert {(g['element_index'], g['local_index']) for g in result['integration_points']} == {(e, i) for e in range(len(mesh['tetrahedra10'])) for i in range(27)}
    report = []
    states = []
    for stage, frame in enumerate(result['frames']):
        state = frame['state']
        assert all(state['active_elements'])
        u = np.asarray(state['displacement_m'])
        stress = np.array([s['effective_stress_compression_positive_kpa'] for s in state['material_states']])
        plastic = np.array([s['plastic_strain_compression_positive'] for s in state['material_states']])
        assert len(stress) == len(gp_z)
        ut, _, _ = theory(xyz[:, 2], stage)
        ug, _, _ = theory(xyz[:, 2], 0)
        _, st, pt = theory(gp_z, stage)
        error = dict(displacement_maximum_m=maximum(u[:, 2], ut-ug), stress_maximum_kpa=maximum(stress, st),
                     plastic_strain_maximum=maximum(plastic, pt))
        gates = dict(displacement=error['displacement_maximum_m'] < bound['native_u_m'],
            stress=error['stress_maximum_kpa'] < bound['native_stress_kpa'],
            plastic_strain=error['plastic_strain_maximum'] < bound['native_plastic_strain'])
        assert maximum(u[:, :2], 0.) < 1e-12
        assert frame['relative_force_residual'] < 1e-8
        assert max(abs(v) for v in frame['force_balance_error_kn']) < 1e-7
        expected_reaction = .8 if stage != 1 else 3.2
        assert abs(frame['support_reaction_kn'][2]-expected_reaction) < 1e-7
        comparisons = []
        for ny, fields in excel_fields.items():
            field = fields[stage]
            du = sample_2d(field, xyz[:, 2], ny)-sample_2d(fields[0], xyz[:, 2], ny)
            gp = sample_2d(field, gp_z, ny, gp=True)
            delta = dict(ny=ny, displacement_maximum_m=maximum(u[:, 2], du),
                stress_maximum_kpa=maximum(stress, gp[:, :6]), plastic_strain_maximum=maximum(plastic, gp[:, 6:]))
            delta['passed'] = (delta['displacement_maximum_m'] < bound['native_u_m']+bound['excel_u_m']
                and delta['stress_maximum_kpa'] < bound['native_stress_kpa']+bound['excel_stress_kpa']
                and delta['plastic_strain_maximum'] < bound['native_plastic_strain']+bound['excel_plastic_strain'])
            comparisons.append(delta)
        report.append(dict(stage=frame['name'], **error, gates=gates, comparisons_to_actual_2d=comparisons,
            top_increment_from_gravity_m=frame['monitors']['top']['mean_displacement_z_m'],
            relative_force_residual=frame['relative_force_residual']))
        states.append(plastic)
    assert np.array_equal(states[1], states[2]), 'native unloading changed plastic history'
    return dict(folder=folder, nodes=len(xyz), elements=len(mesh['tetrahedra10']), gauss_points=len(gp_z),
        passed=all(all(f['gates'].values()) and all(c['passed'] for c in f['comparisons_to_actual_2d']) for f in report),
        material_substep=case['settings']['material_maximum_substep'], construction_step=case['settings']['initial_increment'], frames=report)


def verify(base, output):
    # Bounds are absolute engineering quantities. They are deliberately kept
    # separate from refinement ratios and from the programs' solver tolerances.
    bounds = dict(excel_u_m=2e-12, excel_stress_kpa=2e-8, excel_plastic_strain=2e-10,
                  native_u_m=5e-8, native_stress_kpa=5e-4, native_plastic_strain=1e-6)
    provenance = read(base/'provenance.json')
    assert provenance['revision'] == '440d246d9e32f5892c79539d4a7b797ad29442ee'
    for item in provenance['objects']:
        data = (base/'source'/Path(item['path']).name).read_bytes()
        assert hashlib.sha256(data).hexdigest() == item['sha256']
        assert hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest() == item['git_blob']
    fields = {ny: excel(base, ny, bounds) for ny in (2, 4)}
    results = [native(base, folder, fields, bounds) for folder in ('', '3d-n2-fine', '3d-n4-fine')]
    refinement = []
    for stage in (1, 2):
        coarse, fine, spatial = [r['frames'][stage] for r in results]
        ratios = {k: spatial[k]/fine[k] for k in ('displacement_maximum_m', 'stress_maximum_kpa', 'plastic_strain_maximum')}
        refinement.append(dict(stage=spatial['stage'], spatial_error_ratios=ratios,
            all_spatial_errors_reduced=all(r < .8 for r in ratios.values()),
            n2_material_and_construction_error_ratios={k: fine[k]/coarse[k] for k in ratios}))
    for stage in range(3):
        common = sorted(set(fields[2][stage]['z']) & set(fields[4][stage]['z']))
        assert maximum(sample_2d(fields[2][stage], common, 2), sample_2d(fields[4][stage], common, 4)) < 2e-12
    passed = results[-1]['passed'] and all(r['all_spatial_errors_reduced'] for r in refinement)
    report = dict(schema='soil-fem-actual-2d-3d-column-comparison/1', passed=passed, goal_complete=False, p6_complete=False,
        scope='Pinned Excel Q8 and Rust Tet10: dry confined MC gravity/preload/unload column, phi=10 deg, psi=0',
        reference_origin='Subtract actual Excel gravity displacement; Rust K0 starts at zero displacement with identical initial stress',
        interpolation='Q8 vertical quadratic displacement; elementwise linear GP stress and cumulative plastic strain; no nodal smoothing',
        bounds=bounds, accepted_native_resolution='3d-n4-fine' if passed else None,
        coarse_and_intermediate_are_resolution_diagnostics=True, refinement=refinement,
        historical_assessment='assessment-v1.log stops on the coarse resolution; its unchanged bounds are retained here',
        bounds_unchanged_from_first_assessment=True,
        excel=[dict(vertical_elements=ny, frames=[f['metrics'] for f in fields[ny]]) for ny in (2, 4)], native=results,
        limitations=['One confined MC loading path; not arbitrary 2D/3D equivalence',
                     'No excavation, free lateral boundary, SRM, structures, JOINT or water in this comparison',
                     'Plastic strain is cumulative CSV history; gauss.csv only contains the final increment plastic strain',
                     'Closed-form solution is exactly representable by Q8/Tet10; general spatial convergence remains required'],
        input_and_evidence_sha256={p.relative_to(base).as_posix(): sha(p) for p in base.rglob('*') if p.is_file()},
        verifier_sha256=sha(Path(__file__)))
    with output.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2)
    print(json.dumps({k: report[k] for k in ('passed', 'scope', 'refinement', 'native')}), flush=True)
    assert passed, 'Finest resolution or mesh-refinement gate failed; inspect the written report'


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('base', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    verify(args.base.resolve(), args.output.resolve())
