"""Full-field arc-step sensitivity at every matched saved point, no solve."""
import argparse
import itertools
from pathlib import Path

from audit import compare, stream
from continue_path import read, write


def measure(evidence):
    result = dict(passed=False, history_converged=False,
                  scope='arc-step sensitivity over every matched saved point; not prehistory, material or spatial convergence', maps={})
    for kind in ['native', 'independent']:
        coarse = evidence/f'extended-{kind}-whole'
        fine = evidence/f'extended-{kind}-fine'
        assert read(coarse/'controls.json')['binding']['step_size_m'] == 2*read(fine/'controls.json')['binding']['step_size_m']
        fine_rows = itertools.chain(stream(fine/'samples.jsonl'), stream(evidence/f'extended-{kind}-fine-tail/samples.jsonl'))
        maxima = {}; count = 0
        for a, b in itertools.zip_longest(stream(coarse/'samples.jsonl'), itertools.islice(fine_rows, 1, None, 2)):
            if a is None or b is None or b['step'] != 2*a['step']+1:
                raise ValueError('matched history coverage differs')
            d = compare(a, b, kind == 'native')
            groups = dict(displacement=d['displacement'], **d['material'])
            for name, metrics in groups.items():
                for metric, value in metrics.items():
                    key = name+'.'+metric
                    if key not in maxima or value > maxima[key]['value']:
                        maxima[key] = dict(value=value, coarse_step=a['step'], fine_step=b['step'],
                                           coarse_progress=a['progress'], fine_progress=b['progress'])
            count += 1
        result['maps'][kind] = dict(matched_points=count, field_difference_maxima=maxima)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--evidence', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); write(a.output.resolve(), measure(a.evidence.resolve()))
