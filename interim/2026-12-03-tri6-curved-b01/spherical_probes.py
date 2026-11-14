"""Double-curvature sphere probes, additional coverage within the fixed SPEC."""
import argparse
import copy
from pathlib import Path
import numpy as np
from b01_cases import make, read, write
from study import command, sha


def run(out, cli):
    out.mkdir(parents=True, exist_ok=False)
    rows = []
    vertices = np.array([[0., 0., 2.], [.5, 0., np.sqrt(3.75)], [0., .5, np.sqrt(3.75)]])
    mids = np.asarray([(vertices[i]+vertices[j])/2 for i, j in ((0,1), (1,2), (2,0))])
    mids *= (2/np.linalg.norm(mids, axis=1))[:, None]
    points = np.vstack((vertices, mids))
    for e in (0., .05, -.05):
        f = out/f'sphere-e{e}'
        f.mkdir()
        m, g, c, _ = make('attached', 2)
        cfg = c['beam_system']
        cfg['nodes'] = [dict(id=f's{i}', position_m=p.tolist(), match_soil_node=False) for i, p in enumerate(points)]
        plate = copy.deepcopy(cfg['plates'][0])
        plate.update(id='sphere', nodes=[f's{i}' for i in range(6)], attachment_offset_m=e)
        cfg['plates'] = [plate]
        for stage in (c['stages'][0], c['stages'][-1]):
            stage['members'] = ['sphere']
        for key, value in [('mesh', m), ('ground', g), ('stages', c)]:
            write(f/(key+'.json'), value)
        process = command([cli, 'solve-stages', *(f/(k+'.json') for k in ('mesh','ground','stages')), f/'run.json'], f/'cli.log', expect=None)
        text = (f/'cli.log').read_text()
        rows.append(dict(e_m=e, process=process, diagnostic=text,
            input_sha256={k: sha(f/(k+'.json')) for k in ('mesh','ground','stages')},
            refused=process['exit_code'] == 1 and 'Mindlin Tri6 requires straight edges' in text and not (f/'run.json').exists()))
    write(out/'spherical-refusals.json', dict(passed=False, b01_complete=False, records=rows,
                                           refusal_contract_passed=all(r['refused'] for r in rows)))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('output', type=Path)
    p.add_argument('--cli', type=Path, required=True)
    a = p.parse_args()
    run(a.output.resolve(), a.cli.resolve())
