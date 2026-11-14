"""Input-only full equilibrium controls; no production policy is changed."""
from contextlib import contextmanager
from a03cause_binding import reference, study, fields
import fem_reference as fem


@contextmanager
def volume_quadrature(order, reverse=False):
    if type(order) is not int or order not in (5, 7):
        raise ValueError('diagnostic independent volume order must be 5 or 7')
    original = fem.volume
    if order == 5 and not reverse:
        yield
        return
    def controlled(mesh, ei, independent=False, independent_order=5):
        q = original(mesh, ei, independent, order if independent else independent_order)
        if independent and reverse:
            q = {k:(v[::-1].copy() if k in ('natural','position','weight','shape','b') else v) for k,v in q.items()}
        return q
    fem.volume = controlled
    try:yield
    finally:fem.volume = original


def independent(folder, steps, order=5, reverse=False):
    # Every datum comes from the three input files. No native state is read.
    with volume_quadrature(order, reverse):
        model = reference.Model(*[study.read(p) for p in study.paths(folder)])
        with (folder/'independent-fields.jsonl').open('x', encoding='utf-8') as f, (folder/'independent-prepared.jsonl').open('x', encoding='utf-8') as p:
            result = reference.solve(model, steps, fields.Sink(f), fields.Sink(p, work=False))
    result['diagnostic_volume_order'] = order
    result['diagnostic_reversed_volume_accumulation'] = reverse
    result['observed_volume_quadrature_points_per_element'] = len(model.ops[0][0]['weight'])
    study.write(folder/'independent-result.json', result)
