"""Actual CLI geometry/type refusals and unchanged legacy controls."""
import argparse
import copy
import sys
from pathlib import Path
from b01_cases import B03, make, read, write
from study import command, sha
sys.path.insert(0, str(B03))
from b03_cases import detached_fold, make as b03_make


def run(out, cli):
    out.mkdir(parents=True, exist_ok=False)
    rows = []
    def record(name, data, expected=None, accepted=False):
        m, g, c, _ = data
        f = out/name
        f.mkdir()
        for key, value in [('mesh', m), ('ground', g), ('stages', c)]:
            write(f/(key+'.json'), value)
        target = f/'run.json'
        proc = command([cli, 'solve-stages', *(f/(k+'.json') for k in ('mesh', 'ground', 'stages')), target], f/'cli.log', expect=None)
        log = (f/'cli.log').read_text()
        row = dict(name=name, expected_diagnostic=expected, expected_acceptance=accepted, process=proc,
                   diagnostic=log, result_exists=target.exists(),
                   complete=read(target).get('complete') if target.exists() else None,
                   input_sha256={k: sha(f/(k+'.json')) for k in ('mesh', 'ground', 'stages')})
        row['contract_passed'] = (proc['exit_code'] == 0 and row['complete']) if accepted else (proc['exit_code'] == 1 and expected in log and not target.exists())
        rows.append(row)
        write(out/'scope.json', dict(passed=False, b01_complete=False, records=rows))
    for formulation in ('mindlin_tri6', 'mindlin_tri6_offset'):
        for axis in (0, 2):
            for e in ((0.,) if formulation == 'mindlin_tri6' else (0., .05, -.05)):
                data = detached_fold(0, e, e)
                cfg = data[2]['beam_system']
                for p in cfg['plates']:
                    p['formulation'] = formulation
                    if formulation == 'mindlin_tri6':
                        p.pop('attachment_offset_m')
                node = next(v for v in cfg['nodes'] if v['id'] == cfg['plates'][0]['nodes'][3])
                node['position_m'][axis] += .001
                record(f'curved-mid-{formulation}-axis{axis}-e{e}', data, 'Mindlin Tri6 requires straight edges')
    for e in (0., .05):
        data = make('roof', 2)
        for node in data[2]['beam_system']['nodes']:
            x, y, z = node['position_m']
            radius = (x*x+z*z)**.5
            node['position_m'] = [25*x/radius, y, 25*z/radius]
        for p in data[2]['beam_system']['plates']:
            p['attachment_offset_m'] = e
        record(f'true-cylinder-patch-e{e}', data, 'Mindlin Tri6 requires straight edges')
    for key, value in [('geometry', 'planar'), ('normal', [0., 0., 1.]), ('directors', [[0., 0., 1.]]*6)]:
        data = detached_fold(0, 0, 0)
        cfg = data[2]['beam_system']
        node = next(v for v in cfg['nodes'] if v['id'] == cfg['plates'][0]['nodes'][3])
        node['position_m'][2] += .001
        cfg['plates'][0][key] = value
        record('false-planar-declaration-'+key, data, 'unknown field')
    for e in (-.051, .051):
        record(f'outside-half-thickness-{e}', detached_fold(0, e, e), 'offset')
    for e1, e2, angle, label in [(.05, -.05, 0, 'e-sign'), (.05, .05, 180, 'normal-sign'), (.05, .05, 45, 'curved-facet-nonzero-arm')]:
        record(label, detached_fold(angle, e1, e2), 'same signed normal offset')
    for formulation in ('mindlin_tri6', 'cst_dkt', 'default'):
        data = b03_make('kink-zero')
        for p in data[2]['beam_system']['plates']:
            p.pop('attachment_offset_m', None)
            if formulation == 'default':
                p.pop('formulation', None)
            else:
                p['formulation'] = formulation
            if formulation != 'mindlin_tri6':
                p['nodes'] = p['nodes'][:3]
        record('legacy-e0-'+formulation, data, accepted=True)
    record('planar-offset-positive', b03_make('planar'), accepted=True)
    write(out/'scope.json', dict(passed=False, b01_complete=False, records=rows,
                                all_expected_contracts_passed=all(r['contract_passed'] for r in rows)))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('output', type=Path)
    p.add_argument('--cli', type=Path, required=True)
    a = p.parse_args()
    run(a.output.resolve(), a.cli.resolve())
