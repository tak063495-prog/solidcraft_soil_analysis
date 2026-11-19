"""Solve one 24-Tet10 path using the January 5 fixture and native/independent drivers."""
import argparse
import os
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
PREVIOUS = HERE.parent / '2027-01-05-a02-n8-contact-cap-open-interval'
sys.path.insert(0, str(PREVIOUS))
try:
    from a02n8_run import controlled, fields, native, schedule, study
    from a02n8_common import sha, write
finally:
    sys.path.remove(str(PREVIOUS))


def grid(level):
    if level not in (0, 1, 2, 3):
        raise ValueError('level must be 0 through 3')
    parts = 2 ** level
    original = schedule(8, 8)['open']
    endpoints = []
    for end in original:
        if end in (3/64, 4/64):
            endpoints.extend(end - 1/64 + k/(64*parts) for k in range(1, parts))
        endpoints.append(end)
    return {'open': endpoints}


def run(level, engine, out, build):
    endpoints = grid(level)
    values = study.fixture('stress_free', 8, True)
    assert len(values[0]['tetrahedra']) == 24 and len(values[0]['nodes']) == 81
    values[2]['joint_system']['maximum_increment_m'] = .05
    values[2]['settings']['minimum_increment'] = 1/(64*2**level)
    study.make_inputs(out, values)
    schedule_file = out/'global-input-endpoints.json'
    write(schedule_file, endpoints)
    if engine == 'native':
        os.environ.update(A02N8_P='native', A02N8_R='native', A02N8_MINIMUM='4',
                          A02N8_TRACE='1', A02N8_SCHEDULE=str(schedule_file))
        binary = 'a02n8_probe.exe' if os.name == 'nt' else 'a02n8_probe'
        count = native(out, build/'target/release'/binary, 8, endpoints)
    else:
        with (out/'independent-fields.jsonl').open('x', encoding='utf-8') as stream:
            result = controlled.run(controlled.Model(*values), 8, fields.Sink(stream),
                                    minimum_parts=4, predictor='independent',
                                    repartition='independent', stage_progress=endpoints)
        write(out/'independent-result.json', result)
        assert result['complete']
        count = result['accepted_states']
    assert count == 6*8 + len(endpoints['open'])
    write(out/'completion.json', dict(level=level, engine=engine, predictor=engine,
          repartition=engine, accepted_states=count, minimum_contact_parts=4,
          inputs={p.name: sha(p) for p in [*study.input_paths(out), schedule_file]},
          source_base='c4963221e18534ee87cb7cb2a9c88b5bfeca8755'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--level', type=int, required=True, choices=(1, 2, 3))
    parser.add_argument('--engine', choices=('native', 'independent'), required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--build', type=Path, required=True)
    args = parser.parse_args()
    run(args.level, args.engine, args.output.resolve(), args.build.resolve())
