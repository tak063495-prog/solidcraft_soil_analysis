"""Read authenticated accepted histories; compare common points, never interpolate."""
import io
import json
from pathlib import Path
import zipfile
import numpy as np
from slip_study import ROOT, OBSERVED, sha, decode, write


class OldHistory:
    def __init__(self, key):
        self.root = ROOT/'interim/2026-11-08-wet-joint-finer-full-history'
        self.manifest = json.loads((self.root/'archive-sha256.json').read_text())
        self.members = {}
        self.checked_archives = {}
        self.used_entries = {}
        for row in self.manifest['archives']:
            path = self.root/row['file']
            assert sha(path) == row['sha256']
            self.checked_archives[row['file']] = row['sha256']
            with zipfile.ZipFile(path) as archive:
                for name in archive.namelist():
                    if name.startswith(f'paths/{key}/accepted/'):
                        assert name not in self.members
                        self.members[name] = path
        self.prefix = f'paths/{key}/accepted/'
        self.meta = json.loads(self.read(self.prefix+'manifest.json'))
        assert self.meta['encoding'] == 'float64-le-integer-linear-predictor-xor-byte-shuffle-npz'

    def read(self, name):
        with zipfile.ZipFile(self.members[name]) as archive:
            raw = archive.read(name)
        import hashlib
        digest = hashlib.sha256(raw).hexdigest()
        assert digest == self.manifest['entries'][name]
        self.used_entries[name] = digest
        return raw

    def chunks(self):
        for row in self.meta['chunks']:
            raw = self.read(self.prefix+row['file'])
            with np.load(io.BytesIO(raw), allow_pickle=False) as archive:
                fields = {k: decode(archive[k]) for k in OBSERVED}
            yield row, fields


def new_chunks(folder):
    meta = json.loads((folder/'manifest.json').read_text())
    count = 0
    for row in meta['chunks']:
        assert row['start'] == count and row['count'] == len(row['rows'])
        assert sha(folder/row['file']) == row['sha256']
        with np.load(folder/row['file'], allow_pickle=False) as saved:
            fields = {k: decode(saved[k]) for k in OBSERVED}
        for key in OBSERVED:
            assert fields[key].shape == (row['count'], *({'u':(54,3),'accumulated':(32,)}.get(key,(32,3))))
            assert np.isfinite(fields[key]).all()
        assert all(np.isfinite(r['relative_force_residual']) and r['relative_force_residual'] < 1e-10 for r in row['rows'])
        count += row['count']
        yield row, fields
    assert count == meta['count']


def states(chunks, stride=1):
    for chunk, fields in chunks:
        for i, row in enumerate(chunk['rows']):
            if (chunk['start']+i+1) % stride == 0:
                yield row, {k:v[i] for k,v in fields.items()}


def compare(coarse, fine, count=None):
    """count permits a recorded prefix; it is explicitly not a complete history."""
    from itertools import zip_longest, islice
    rows = {}
    common = 0
    for a, b in zip_longest(islice(coarse, count), islice(fine, count)):
        assert a is not None and b is not None, 'missing matching accepted point'
        ra, fa = a
        rb, fb = b
        assert all(ra[k] == rb[k] for k in ('stage','progress','factor')), 'changed path'
        common += 1
        stage = rows.setdefault(ra['stage'], dict(count=0, fields={}, status_count_mismatches=0))
        stage['count'] += 1
        stage['status_count_mismatches'] += any(ra[k] != rb[k] for k in ('open_points','slipping_points'))
        for key in OBSERVED:
            delta = fa[key]-fb[key]
            diff2 = float(np.sum(delta**2))
            ref2 = float(np.sum(fb[key]**2))
            relative = float(np.sqrt(diff2)/max(np.sqrt(ref2), 1e-12))
            entry = stage['fields'].setdefault(key, dict(difference_squared=0., reference_squared=0.,
                maximum_absolute=0., maximum_state_relative_l2=-1.))
            entry['difference_squared'] += diff2
            entry['reference_squared'] += ref2
            entry['maximum_absolute'] = max(entry['maximum_absolute'], float(abs(delta).max()))
            if relative > entry['maximum_state_relative_l2']:
                entry['maximum_state_relative_l2'] = relative
                entry['maximum_at_progress'] = ra['progress']
    for stage in rows.values():
        for entry in stage['fields'].values():
            entry['pooled_relative_l2'] = float(np.sqrt(entry['difference_squared'])/max(np.sqrt(entry['reference_squared']),1e-12))
    return dict(common_accepted_states=common, stages=rows, interpolation=False,
        partial_prefix_only=count is not None, passed=False, history_converged=False)


def audit_new(folder, steps):
    meta = json.loads((folder/'manifest.json').read_text())
    last = None
    for chunk, fields in new_chunks(folder):
        for row in chunk['rows']:
            last = row
    assert last is not None
    return dict(flushed_accepted_states=meta['count'], last_flushed_accepted=last,
        complete_accepted_sequence=meta['count'] == steps*7,
        fields=list(OBSERVED), bulk_fields_at_every_accepted_state_saved=False,
        partial_state_is_resumable_checkpoint=False)


def measure(folder):
    paths = folder/'paths'
    old = OldHistory('g4096-l64')
    results = {}
    for name in ('g512-l64','g1024-l64'):
        terminal = json.loads((paths/name/'terminal.json').read_text())
        if terminal['complete']:
            steps = int(name.split('-')[0][1:])
            results[name+'__g4096-l64'] = compare(states(new_chunks(paths/name/'accepted')),
                                                 states(old.chunks(),4096//steps))
    if (paths/'g512-l128/terminal.json').exists():
        term = json.loads((paths/'g512-l128/terminal.json').read_text())
        if term['complete']:
            results['g512-l64__g512-l128'] = compare(states(new_chunks(paths/'g512-l64/accepted')),
                                                      states(new_chunks(paths/'g512-l128/accepted')))
    fine = paths/'g4096-l128'
    if (fine/'terminal.json').exists():
        meta = json.loads((fine/'accepted/manifest.json').read_text())
        results['g4096-l64__g4096-l128_flushed_prefix'] = compare(states(old.chunks()),
            states(new_chunks(fine/'accepted')), meta['count'])
    checks = {p.name:audit_new(p/'accepted',int(p.name.split('-')[0][1:])) for p in paths.glob('g*') if p.is_dir()}
    return dict(comparisons=results, saved_history_checks=checks,
        old_history_provenance=dict(archives=old.checked_archives, entries=old.used_entries),
        saved_contact_observations_only=True, all_18_fields_at_endpoints_only=True,
        all_history_acceptance=False)


if __name__ == '__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--evidence',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    write(args.output,measure(args.evidence.resolve()))
