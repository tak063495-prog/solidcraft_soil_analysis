"""Independent volume/face controls and stable physical correspondence."""
from contextlib import contextmanager
import hashlib
import json
import math
import numpy as np
from a03follow_binding import fields, reference, study
import fem_reference as fem

SHAPES = {'wide': [1.4, .9, 1.], 'tall': [.8, 1.1, 1.2]}
STRESSES = {'k05': dict(k0=.5, lower=90., upper=70., fresh=5.),
            'k08': dict(k0=.8, lower=110., upper=85., fresh=8.)}


def fixture(n, shape, stress, element_order=False, face_order=False):
    m, g, c = study.fixture(n, 2.)
    scale = np.asarray(SHAPES[shape])
    m['nodes'] = (np.asarray(m['nodes'])*scale).tolist()
    for region in m['regions']:
        region['bounds_m'] = (np.asarray(region['bounds_m'])*scale).tolist()
        region['volume_m3'] *= math.prod(scale)
        region['source_sha256'] = hashlib.sha256(json.dumps(dict(
            diagnostic='a03follow affine input', scale=scale.tolist(), region=region['id']), sort_keys=True).encode()).hexdigest()
    s = STRESSES[stress]
    def sigma(z):
        return [s['k0']*z, s['k0']*z, z, 0., 0., 0.]
    states = g['initialization']['geostatic']['seed']['states']
    states['LowerClay']['effective_stress_compression_positive_kpa'] = sigma(s['lower'])
    states['UpperClay']['effective_stress_compression_positive_kpa'] = sigma(s['upper'])
    c['stages'][0]['states']['UpperClay']['effective_stress_compression_positive_kpa'] = sigma(s['fresh'])
    if element_order:
        for key in ('tetrahedra', 'tetrahedra10', 'material_ids', 'region_ids'):
            m[key] = list(reversed(m[key]))
    if face_order:
        m['face_sets'] = {k: list(reversed(v)) for k, v in m['face_sets'].items()}
    return m, g, c


@contextmanager
def controls(volume_order=5, reverse_faces=False, minimum=1):
    """No global control survives the explicit diagnostic context."""
    if type(volume_order) is not int or volume_order not in (5, 7):
        raise ValueError('GL5 or GL7 required')
    if type(minimum) is not int or minimum not in (1, 2, 4):
        raise ValueError('local minimum must be 1/2/4')
    volume, face, jacobi, work = fem.volume, fem.face, fem.jacobi, reference.contact_work
    def controlled_volume(mesh, ei, independent=False, independent_order=5):
        return volume(mesh, ei, independent, volume_order if independent else independent_order)
    def controlled_face(*args, **kwargs):
        q = face(*args, **kwargs)
        return {k: v[::-1].copy() if k != 'nodes' else v for k, v in q.items()}
    def controlled_jacobi(n, alpha):
        x, w = jacobi(n, alpha)
        return (x[::-1].copy(), w[::-1].copy()) if n == 4 and alpha in (0, 1) else (x, w)
    def controlled_work(system, old, traction, target, count):
        return work(system, old, traction, target, max(minimum, count))
    fem.volume = controlled_volume
    reference.contact_work = controlled_work
    if reverse_faces:
        fem.face, fem.jacobi = controlled_face, controlled_jacobi
    try:
        yield
    finally:
        fem.volume, fem.face, fem.jacobi, reference.contact_work = volume, face, jacobi, work


def dot(a, b):
    a, b = np.asarray(a, float).ravel(order='C'), np.asarray(b, float).ravel(order='C')
    if a.shape != b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('invalid finite correspondence')
    return math.fsum(float(x)*float(y) for x, y in zip(a, b))


def norm(a):
    return math.sqrt(dot(a, a))


def metric(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if a.shape != b.shape:
        raise ValueError('shape mismatch')
    d = a-b
    raw = norm(b)
    absolute = norm(d)
    return dict(absolute_l2=absolute, reference_l2=raw, denominator=max(raw, 1e-12),
                relative_l2=absolute/max(raw, 1e-12), maximum_absolute=float(np.max(np.abs(d), initial=0.)),
                denominator_floor_active=raw < 1e-12, components=int(a.size))


class Sink:
    """Same endpoint rule; only the observation summation is fixed with fsum."""
    def __init__(self, stream, work=True):
        self.stream, self.include_work = stream, work
        self.previous = None
        self.load_increments, self.reaction_increments = [], []

    def placement(self, u):
        self.previous = np.asarray(u).copy()

    def __call__(self, row, values):
        values = {k: np.asarray(v) for k, v in values.items()}
        if self.include_work:
            old = np.zeros_like(values['u']) if self.previous is None else self.previous
            du = values['u']-old
            for kind, increments in (('load', self.load_increments), ('reaction', self.reaction_increments)):
                source = values['applied_load' if kind == 'load' else 'reaction']
                product = source*du
                values['nodal_'+kind+'_work_increment'] = product
                increments.append(math.fsum(map(float, product.ravel(order='C'))))
                values[kind+'_work_cumulative'] = np.asarray(math.fsum(increments))
            self.previous = values['u'].copy()
        self.stream.write(json.dumps(dict(row=row, fields={k: v.tolist() for k, v in values.items()}),
                                     allow_nan=False, separators=(',', ':'))+'\n')
        self.stream.flush()


MATERIAL = set('stress strain void total_stress pore_pressure head plastic_strain equivalent_plastic_strain plastic_work material_return'.split())
CONTACT = set('jump traction reference_jump initial_traction plastic_slip accumulated plastic_dissipation numerical_dissipation endpoint_work elastic_offset initial_energy contact_status contact_open opening gap normal_pressure shear_traction elastic_energy joint_water_pressure'.split())


def layout(model):
    """Keys use exact connectivity and local rule index, never nearest points."""
    elements = [tuple(sorted(cell)) for cell in model.mesh['tetrahedra10']]
    eorder = sorted(range(len(elements)), key=elements.__getitem__)
    material = [e*27+i for e in eorder for i in range(27)]
    # A rule reversal is known exactly, so undo it by the local shape vector.
    contact_keys = [(tuple(sorted(p['nodes_a'])), tuple(sorted(p['nodes_b'])), tuple(p['shape']))
                    for p in model.system.points]
    if len(set(contact_keys)) != len(contact_keys):
        raise ValueError('duplicate physical contact key')
    contact = sorted(range(len(contact_keys)), key=contact_keys.__getitem__)
    face_keys = [tuple(sorted(model.mesh['face_sets'][model.case['joint_system']['joints'][0]['side_a']][i]))
                 for i in range(len(model.mesh['face_sets'][model.case['joint_system']['joints'][0]['side_a']]))]
    return dict(element_order=eorder, material_order=material, contact_order=contact,
                joint_face_order=sorted(range(len(face_keys)), key=face_keys.__getitem__),
                element_keys=[list(elements[i]) for i in eorder],
                contact_keys=[contact_keys[i] for i in contact])


def canonical(values, mapping):
    result = {}
    for key, value in values.items():
        a = np.asarray(value)
        if key in MATERIAL:
            a = a[mapping['material_order']]
        elif key in CONTACT:
            a = a[mapping['contact_order']]
        elif key == 'active_elements':
            a = a[mapping['element_order']]
        elif key == 'joint_face_water_force':
            a = a[mapping['joint_face_order']]
        result[key] = a
    return result


def compare(a, b):
    if not a or len(a) != len(b):
        raise ValueError('state count differs')
    for (ra, fa), (rb, fb) in zip(a, b):
        if ra != rb or set(fa) != set(fb):
            raise ValueError('exact coordinates/fields differ; interpolation forbidden')
    result = {}
    for key in sorted(a[0][1]):
        x, y = (np.stack([f[key] for _, f in rows]) for rows in (a, b))
        result[key] = (dict(mismatches=int(np.count_nonzero(x != y)), components=int(x.size))
                       if key in fields.CATEGORIES else metric(x, y))
    return dict(states=len(a), fields=result, numerical_acceptance=False)


def common(rows, n):
    selected = [r for r in rows if r[0]['progress']*n == int(r[0]['progress']*n)]
    if len(selected) != 4*n:
        raise ValueError('incomplete common coordinates')
    return selected
