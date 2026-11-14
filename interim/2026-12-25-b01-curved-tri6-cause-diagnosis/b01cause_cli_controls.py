"""Replay the exact saved 21 refusals and four existing accepted CLI controls."""
import argparse
import copy
import hashlib
import io
import json
from pathlib import Path
import subprocess
import time
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PRIOR = ROOT/'interim/2026-12-03-tri6-curved-b01'


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def clean_run(obj):
    def without_time(value):
        if isinstance(value, dict):
            return {k: without_time(v) for k,v in value.items() if k != 'elapsed_seconds'}
        if isinstance(value, list):
            return [without_time(v) for v in value]
        return value
    obj = without_time(copy.deepcopy(obj))
    # Different executable identities are separately authenticated. Within one
    # executable, no such exclusion is used by the historical restart tests.
    obj.get('provenance', {}).pop('solver_binary_sha256', None)
    return obj


def run(output, solver):
    output.mkdir(parents=True, exist_ok=False)
    manifest = json.loads((PRIOR/'evidence-manifest.json').read_text())
    pieces = []
    for part in manifest['parts']:
        raw = (PRIOR/part['path']).read_bytes()
        if len(raw) != part['bytes'] or digest(raw) != part['sha256']:
            raise ValueError('prior archive part mismatch')
        pieces.append(raw)
    raw = b''.join(pieces)
    if len(raw) != manifest['archive_bytes'] or digest(raw) != manifest['archive_sha256']:
        raise ValueError('prior combined archive mismatch')
    rows = []; old = json.loads((PRIOR/'scope-refusals.json').read_text())['records']
    for a in old:
        rows.append(dict(name=a['name'], prefix='probes/'+a['name'], accepted=a['expected_acceptance'],
                         expected_diagnostic=a['expected_diagnostic'], input_sha256=a['input_sha256']))
    for a in json.loads((PRIOR/'spherical-refusals.json').read_text())['records']:
        name='sphere-e'+str(a['e_m'])
        rows.append(dict(name=name, prefix='probes/spherical/'+name, accepted=False,
                         expected_diagnostic='Mindlin Tri6 requires straight edges', input_sha256=a['input_sha256']))
    results = []
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        if set(z.namelist()) != {m['path'] for m in manifest['members']}:
            raise ValueError('prior member inventory differs')
        for m in manifest['members']:
            b = z.read(m['path'])
            if len(b) != m['bytes'] or digest(b) != m['sha256']:
                raise ValueError('prior member differs: '+m['path'])
        for row in rows:
            folder = output/row['name']; folder.mkdir()
            for key in ('mesh','ground','stages'):
                b = z.read(row['prefix']+'/'+key+'.json')
                if digest(b) != row['input_sha256'][key]:
                    raise ValueError('saved input binding differs')
                (folder/(key+'.json')).write_bytes(b)
            command = [str(solver), 'solve-stages', *[str(folder/(k+'.json')) for k in ('mesh','ground','stages')], str(folder/'run.json')]
            t = time.monotonic()
            with (folder/'cli.log').open('w', encoding='utf-8') as log:
                p = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT); code = p.wait()
            text = (folder/'cli.log').read_text(encoding='utf-8')
            exists = (folder/'run.json').exists()
            value = dict(**row, command=command, pid=p.pid, exit_code=code,
                         elapsed_seconds=time.monotonic()-t, result_exists=exists,
                         diagnostic=text, log_sha256=digest((folder/'cli.log').read_bytes()))
            if row['accepted']:
                current = json.loads((folder/'run.json').read_text()) if exists else {}
                historical = json.loads(z.read(row['prefix']+'/run.json'))
                value.update(complete=current.get('complete'),
                             original_run_equal_except_authenticated_binary_and_time=json.dumps(clean_run(current), sort_keys=True)==json.dumps(clean_run(historical), sort_keys=True),
                             old_solver_binary_sha256=historical.get('provenance',{}).get('solver_binary_sha256'),
                             new_solver_binary_sha256=current.get('provenance',{}).get('solver_binary_sha256'),
                             contract_confirmed=code==0 and current.get('complete') is True)
            else:
                value['contract_confirmed'] = code==1 and not exists and row['expected_diagnostic'] in text
            results.append(value)
    report = dict(diagnostic_only=True, production_scope_expanded=False,
                  solver_sha256=digest(solver.read_bytes()), prior_archive_sha256=manifest['archive_sha256'],
                  prior_all_members_sha256_verified=True, records=results,
                  refusals=sum(not r['accepted'] for r in results), accepted_controls=sum(r['accepted'] for r in results),
                  all_original_contracts_confirmed=all(r['contract_confirmed'] for r in results))
    (output/'cli-controls-observations.json').write_text(json.dumps(report, indent=2)+'\n')
    print({k:v for k,v in report.items() if k!='records'})


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True); p.add_argument('--solver', type=Path, required=True)
    a=p.parse_args(); run(a.output.resolve(), a.solver.resolve())
