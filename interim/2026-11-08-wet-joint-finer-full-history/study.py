"""Existing small wet JOINT only; save accepted fields, never certify history."""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[_key] = '1'
import argparse
import hashlib
import io
import json
from pathlib import Path
import sys
import time
import zipfile
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from diagnose_cut_volume_history import FIELDS, inputs, validate_path, summary, write
from diagnose_wet_joint_simultaneous import recorded
from joint_water_reference import Water
from wet_joint_cycle_reference import differences, serializable
from reference_history import cut_volume_reference

BASE = 'c4de9fdafa4d7147f4f2bd8c8c877ca969b269b6'
PRIOR = ROOT/'interim/2026-11-08-wet-joint-simultaneous-force'
ROUTES = ((2048, 64), (4096, 64))


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


ENCODING = 'float64-le-integer-linear-predictor-xor-byte-shuffle-npz'
LEGACY_ENCODING = 'float64-le-xor-time-byte-shuffle-npz'


def encode(array):
    """Lossless prediction on uint64 bits modulo 2**64; no float arithmetic."""
    raw = np.ascontiguousarray(array, dtype='<f8').view('<u8')
    delta = raw.copy()
    delta[1:2] = raw[1:2] ^ raw[:1]
    delta[2:] = raw[2:] ^ (raw[1:-1]+(raw[1:-1]-raw[:-2]))
    return np.ascontiguousarray(np.moveaxis(delta.view('u1').reshape(*raw.shape, 8), -1, 0))


def decode(array, encoding=ENCODING):
    assert array.dtype == np.uint8 and array.shape[0] == 8
    shape = array.shape[1:]
    delta = np.ascontiguousarray(np.moveaxis(array, 0, -1)).view('<u8').reshape(shape)
    if encoding == LEGACY_ENCODING:
        return np.bitwise_xor.accumulate(delta, axis=0).view('<f8')
    assert encoding == ENCODING
    raw = delta.copy()
    if len(raw) > 1:
        raw[1] ^= raw[0]
    for i in range(2, len(raw)):
        raw[i] ^= raw[i-1]+(raw[i-1]-raw[i-2])
    return raw.view('<f8')


class HistoryWriter:
    def __init__(self, folder, chunk_size=128):
        self.folder = folder
        folder.mkdir(parents=True, exist_ok=False)
        self.chunk_size = chunk_size
        self.rows = []; self.fields = []; self.chunks = []; self.count = 0

    def __call__(self, row, fields):
        if self.rows and self.rows[-1]['stage'] != row['stage']:
            self.flush()
        self.rows.append(row)
        self.fields.append({key: np.asarray(fields[key]).copy() for key in FIELDS})
        if len(self.rows) == self.chunk_size:
            self.flush()

    def flush(self):
        if not self.rows:
            return
        name = f'{len(self.chunks):04d}.npz'
        arrays = {key: encode(np.stack([f[key] for f in self.fields])) for key in FIELDS}
        with (self.folder/name).open('xb') as stream:
            np.savez_compressed(stream, **arrays)
        raw = (self.folder/name).read_bytes()
        self.chunks.append(dict(file=name, sha256=sha(raw), start=self.count,
                                count=len(self.rows), rows=self.rows, bytes=len(raw)))
        self.count += len(self.rows)
        print('saved', self.folder.parent.name, self.rows[-1]['stage'],
              self.rows[-1]['progress'], 'accepted', self.count, flush=True)
        self.rows = []; self.fields = []

    def close(self):
        self.flush()
        manifest = dict(schema='wet-joint-accepted-fields/1', encoding=ENCODING,
                        fields={k: list(v) for k, v in FIELDS.items()}, count=self.count, chunks=self.chunks)
        write(self.folder/'manifest.json', manifest)
        return manifest


def check_fields(fields, system):
    """Apply existing endpoint invariants to every saved accepted state."""
    n = fields['u'].shape[0]
    assert set(fields) == set(FIELDS)
    for key, shape in FIELDS.items():
        assert fields[key].shape == (n, *shape) and np.isfinite(fields[key]).all(), key
    jump = fields['jump']; q = fields['traction']
    gap = np.einsum('tqa,qa->tq', jump, system.normal)
    qn = np.einsum('tqa,qa->tq', q, system.normal)
    tau = q-qn[:, :, None]*system.normal
    np.testing.assert_allclose(np.einsum('qag,tg->tqa', system.c, fields['u'].reshape(n, -1)), jump, atol=1e-10, rtol=0)
    np.testing.assert_allclose(fields['gap'], gap, atol=1e-10, rtol=0)
    np.testing.assert_allclose(qn, system.kn*np.minimum(gap, 0), atol=1e-7, rtol=1e-9)
    assert np.all(np.linalg.norm(tau, axis=2) <= system.mu*np.maximum(-qn, 0)+1e-7)
    energy = qn**2/(2*system.kn)+np.einsum('tqa,tqa->tq', tau, tau)/(2*system.kt)
    np.testing.assert_allclose(fields['elastic_energy'], energy, atol=1e-9, rtol=1e-8)
    np.testing.assert_allclose(fields['endpoint_work'], energy+fields['plastic_dissipation']+
                               fields['numerical_dissipation'], atol=1e-9, rtol=1e-8)
    for key in ('plastic_dissipation', 'numerical_dissipation', 'accumulated'):
        assert np.min(fields[key]) >= -1e-12
    assert np.all(fields['accumulated']+1e-10 >= np.linalg.norm(fields['plastic_slip'], axis=2))
    assert not fields['plastic_strain'].any()


def check_history(folder, result, values):
    manifest = json.loads((folder/'manifest.json').read_text())
    assert manifest['fields'] == {k: list(v) for k, v in FIELDS.items()}
    assert manifest['encoding'] in (ENCODING, LEGACY_ENCODING)
    assert manifest['count'] == len(result['path']) == result['steps']*7
    system = Water(*values).contact
    count = 0; endpoints = 0; maximum = 0.
    for chunk in manifest['chunks']:
        assert chunk['start'] == count and chunk['count'] > 0
        assert chunk['rows'] == result['path'][count:count+chunk['count']], 'accepted metadata differs'
        raw = (folder/chunk['file']).read_bytes()
        assert sha(raw) == chunk['sha256'], 'changed accepted chunk'
        with np.load(io.BytesIO(raw), allow_pickle=False) as archive:
            fields = {key: decode(archive[key], manifest['encoding']) for key in archive.files}
        check_fields(fields, system)
        assert fields['u'].shape[0] == chunk['count']
        for i, row in enumerate(chunk['rows']):
            assert np.isfinite(row['relative_force_residual']) and row['relative_force_residual'] < 1e-10
            maximum = max(maximum, row['relative_force_residual'])
            if row['progress'] == 1.:
                frame = result['frames'][endpoints]
                assert frame['name'] == ('initial_equilibrated' if row['stage'] == 'geostatic' else row['stage'])
                for key in FIELDS:
                    # Verify exact binary reconstruction, including signed zero.
                    expected = np.ascontiguousarray(frame[key], dtype='<f8').view('<u8')
                    np.testing.assert_array_equal(fields[key][i].view('<u8'), expected)
                endpoints += 1
        count += chunk['count']
    assert count == manifest['count'] and endpoints == 7
    return dict(accepted_states=count, chunks=len(manifest['chunks']), fields_per_state=len(FIELDS),
                material_points_per_state=324, contact_points_per_state=32,
                endpoint_fields_bitwise_equal=True, all_saved_field_invariants_checked=True,
                maximum_relative_force_residual=maximum, full_history_acceptance=False)


def accepted_states(folder, stride=1):
    manifest = json.loads((folder/'manifest.json').read_text())
    for chunk in manifest['chunks']:
        raw = (folder/chunk['file']).read_bytes()
        assert sha(raw) == chunk['sha256']
        with np.load(io.BytesIO(raw), allow_pickle=False) as archive:
            fields = {key: decode(archive[key], manifest['encoding']) for key in archive.files}
        for i, row in enumerate(chunk['rows']):
            if (chunk['start']+i+1) % stride == 0:
                yield row, {key: value[i] for key, value in fields.items()}


def compare_accepted(coarse, fine, ratio=2):
    """All common accepted progress points; measurements, with no pass gate."""
    from itertools import zip_longest
    sums = {key: dict(difference_squared=0., reference_squared=0., maximum_absolute=0.,
                      maximum_state_relative_l2=0., maximum_state=None) for key in FIELDS}
    count = 0
    for a, b in zip_longest(accepted_states(coarse), accepted_states(fine, ratio)):
        assert a is not None and b is not None, 'missing common state'
        ra, fa = a; rb, fb = b
        assert all(ra[k] == rb[k] for k in ('stage', 'progress', 'factor')), 'misaligned history'
        count += 1
        for key in FIELDS:
            delta = fa[key]-fb[key]
            d2 = float(np.sum(delta*delta)); r2 = float(np.sum(fb[key]*fb[key]))
            row = sums[key]
            row['difference_squared'] += d2; row['reference_squared'] += r2
            row['maximum_absolute'] = max(row['maximum_absolute'], float(np.max(np.abs(delta))))
            relative = float(np.sqrt(d2)/max(np.sqrt(r2), 1e-12))
            if row['maximum_state'] is None or relative > row['maximum_state_relative_l2']:
                row['maximum_state_relative_l2'] = relative
                row['maximum_state'] = {k: ra[k] for k in ('stage', 'progress', 'factor')}
    for row in sums.values():
        row['pooled_relative_l2'] = float(np.sqrt(row['difference_squared'])/
                                        max(np.sqrt(row['reference_squared']), 1e-12))
    return dict(common_accepted_states=count, finer_unmatched_states_compared=False,
                fields=sums, full_history_acceptance=False,
                limitation='Measurements at common accepted progress only; no interpolated fields or new acceptance thresholds.')


def measure(folder):
    values = inputs()
    for name, expected in zip(('mesh', 'ground', 'stages'), values):
        assert json.loads((folder/(name+'.json')).read_text()) == expected
    paths = {}; checks = {}; histories = {}
    for g, local in ROUTES:
        key = f'g{g}-l{local}'
        result = json.loads((folder/key/'result.json').read_text())
        checks[key] = validate_path(result, *values, g, local)
        histories[key] = check_history(folder/key/'accepted', result, values)
        paths[key] = result
    # Historical mechanical fields are first read after both input-only solves.
    old, provenance = recorded(PRIOR, ['paths/g2048-l32.json', 'paths/g1024-l32.json'])
    for name, result in old.items():
        key = Path(name).stem
        validate_path(result, *values, result['steps'], result['contact_steps'])
        paths[key] = result
    comparisons = {}
    for a, b in (('g2048-l32', 'g2048-l64'), ('g2048-l64', 'g4096-l64'),
                 ('g2048-l32', 'g4096-l64'), ('g1024-l32', 'g2048-l64')):
        rows = differences(paths[a], paths[b])
        comparisons[a+'__'+b] = dict(fields=rows, maxima=summary(rows))
    return dict(schema='soil-fem-wet-joint-finer-full-history/1', starting_commit=BASE,
                passed=False, history_converged=False, production_adopted=False,
                diagnostic_completed=True, reference_uses_native_results=False,
                native_rerun=False, nc_rerun=False, baseline=provenance,
                timings=json.loads((folder/'timings.json').read_text()),
                path_field_checks=checks, accepted_history_checks=histories, comparisons=comparisons,
                common_accepted_comparison=compare_accepted(folder/'g2048-l64/accepted', folder/'g4096-l64/accepted'),
                full_history_acceptance=dict(accepted=False, reasons=[
                    'No existing all-history error criterion certifies these measured differences.',
                    'Old baselines contain seven endpoint fields only; full accepted-field pair comparison is absent.',
                    'Production analytic volume integration and restart through that path remain unverified.']),
                windows_numerical_reproduction='unverified', p3b_p4_complete=False,
                structure_water_joint_complete=False, p6_complete=False)


def run(output):
    output.mkdir(parents=True, exist_ok=False)
    values = inputs()
    for name, value in zip(('mesh', 'ground', 'stages'), values):
        write(output/(name+'.json'), value)
    timings = []
    for steps, local in ROUTES:
        key = f'g{steps}-l{local}'
        sink = HistoryWriter(output/key/'accepted')
        start = time.monotonic()
        result = cut_volume_reference(*values, steps, local, 7, history_sink=sink)
        sink.close()
        validate_path(result, *values, steps, local)
        write(output/key/'result.json', result)
        timings.append(dict(route=key, elapsed_seconds=time.monotonic()-start, complete=True))
        print('completed', timings[-1], flush=True)
    write(output/'timings.json', timings)
    result = measure(output)
    write(output/'measurements.json', result)
    return result


def remeasure(evidence):
    manifest = json.loads((HERE/'archive-sha256.json').read_text())
    for name, digest in manifest['entries'].items():
        assert sha((evidence/name).read_bytes()) == digest, ('changed evidence', name)
    actual = serializable(measure(evidence/'paths'))
    expected = json.loads((evidence/'paths/measurements.json').read_text())
    assert actual == expected, 'measurement mismatch'
    assert not actual['passed'] and not actual['history_converged'] and not actual['production_adopted']
    return dict(software_check_passed=True, numerical_acceptance_passed=False,
                all_entries_sha256_verified=len(manifest['entries']), solver_rerun=False,
                all_measurements_exactly_reproduced=True,
                accepted_history_checks=actual['accepted_history_checks'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--remeasure', type=Path)
    args = parser.parse_args()
    if args.remeasure:
        result = remeasure(args.remeasure.resolve()); write(args.output, result)
    else:
        result = run(args.output.resolve())
    print(json.dumps({k: result[k] for k in result if k in ('passed', 'history_converged', 'diagnostic_completed', 'all_measurements_exactly_reproduced')}))
