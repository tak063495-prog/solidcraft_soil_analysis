"""Compare raw frozen runs and record failed gates without changing them."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np

p = Path(__file__).resolve().parent
sys.path.insert(0, str(p/'runtime/tools'))
from richards_column_reference import retention
read = lambda q: json.loads(q.read_text('utf-8-sig'))
sha = lambda q: hashlib.sha256(q.read_bytes()).hexdigest()
rms = lambda x: float(np.sqrt(np.mean(np.asarray(x)**2)))


def interpolate(field, target):
    """Trilinear cell-centre interpolation, constant extension at the boundary."""
    n = field.shape[0]
    value = field
    for axis in range(3):
        coordinate = np.clip(target*n-.5, 0., n-1.)
        lower = np.minimum(np.floor(coordinate).astype(int), n-2)
        upper = lower+1
        fraction = coordinate-lower
        shape = [1,1,1]
        shape[axis] = len(target)
        value = np.take(value, lower, axis=axis)*(1-fraction.reshape(shape)) + np.take(value, upper, axis=axis)*fraction.reshape(shape)
    return value


def endpoints(run):
    return {name: next(f for f in reversed(run['frames']) if f['name']==name) for name in ('rain','pause','rain_again')}


def native(run, n, name):
    f = endpoints(run)[name]
    h = np.asarray(f['state']['relative_cell_head_m'])
    assert h.size == 6*n**3
    return h.reshape(n,n,n,6).mean(axis=3), float(sum(f['state']['cell_storage_offset_m3']))


def front(head, law):
    # Diagnostic from cube-average P0 head and cell-centre z, not a GP field.
    n = len(head)
    z = (np.arange(n)+.5)/n
    delta = retention(head-z[:,None,None],law)[0]-retention(-z[:,None,None],law)[0]
    depth = np.zeros((n,n))
    for y in range(n):
        for x in range(n):
            wet = np.flatnonzero(delta[:,y,x]>=.005)
            if not wet.size:
                continue
            i = wet[0]
            height = z[0] if i==0 else z[i-1]+(.005-delta[i-1,y,x])/(delta[i,y,x]-delta[i-1,y,x])*(z[i]-z[i-1])
            depth[y,x] = 1-height
    return depth


assert read(p/'raw-terminal.json')['passed']
scope = read(p/'scope.json')
gates = scope['predeclared_gates']
solver_sha = scope['solver_sha256']
refs = {(n,dt):read(p/f'fv-{n}-dt{dt}.json') for n,dt in [(8,100),(16,100),(32,100),(16,50),(16,25)]}
runs = {(n,dt):read(p/f'rt0-{n}-dt{dt}/result.json') for n,dt in [(4,100),(8,100),(16,100),(8,50),(8,25)]}
for (n,dt),run in runs.items():
    assert run['complete'] and run['solver_binary_sha256']==solver_sha
    mesh = read(p/f'rt0-{n}-dt{dt}/mesh.json')
    xyz = np.asarray(mesh['nodes'])
    areas = {name:sum(np.linalg.norm(np.cross(xyz[f[1]]-xyz[f[0]],xyz[f[2]]-xyz[f[0]]))/2 for f in mesh['face_sets'][name]) for name in ['rain','dry','bottom']}
    np.testing.assert_allclose([areas['rain'],areas['dry'],areas['bottom']],[.25,.75,1.],rtol=0,atol=1e-12)
    assert read(p/f'rt0-{n}-dt{dt}/audit.json')['passed']
checks, comparisons, refinements = [], [], []


def check(name, actual, limit):
    checks.append(dict(name=name, actual=float(actual), upper_bound=float(limit), passed=bool(np.isfinite(actual) and actual<limit)))


sample = (np.arange(32)+.5)/32
law = read(p/'reference-case.json')['materials'][0]['retention']
for stage in ['rain','pause','rain_again']:
    rf = next(f for f in refs[32,100]['frames'] if f['name']==stage)
    reference = np.asarray(rf['head_m'])
    for n in [4,8,16]:
        h, storage = native(runs[n,100],n,stage)
        href = interpolate(reference,(np.arange(n)+.5)/n)
        row = dict(stage=stage,native_subdivisions=n,reference_subdivisions=32,dt_s=100,
            head_rms_error_m=rms(h-href),head_maximum_error_m=float(np.max(np.abs(h-href))),
            storage_error_m3=abs(storage-rf['storage_offset_m3']),native_storage_offset_m3=storage,
            reference_storage_offset_m3=rf['storage_offset_m3'])
        comparisons.append(row)
        if n==16:
            check(stage+' finest head RMS',row['head_rms_error_m'],gates['finest_head_rms_m'])
            check(stage+' finest storage',row['storage_error_m3'],gates['finest_storage_error_m3'])
            fd = front(h,law)
            rd = front(reference,law)
            row.update(native_head_derived_front_depth_central_m=float(fd[n//2,n//2]),
                reference_head_derived_front_depth_central_m=float(rd[16,16]),
                front_diagnostic_is_not_gp_average=True)
    ns = [interpolate(native(runs[n,100],n,stage)[0],sample) for n in [4,8,16]]
    fs = [interpolate(np.asarray(next(f for f in refs[n,100]['frames'] if f['name']==stage)['head_m']),sample) for n in [8,16,32]]
    nt = [native(runs[8,dt],8,stage)[0] for dt in [100,50,25]]
    ft = [np.asarray(next(f for f in refs[16,dt]['frames'] if f['name']==stage)['head_m']) for dt in [100,50,25]]
    row = dict(stage=stage)
    for label,fields in [('native_space',ns),('fv_space',fs),('native_time',nt),('fv_time',ft)]:
        differences = [rms(fields[0]-fields[1]),rms(fields[1]-fields[2])]
        row[label+'_successive_rms_m'] = differences
        ratio = differences[1]/differences[0] if differences[0]>1e-16 else 1e300
        check(stage+' '+label+' difference ratio',ratio,gates['successive_space_and_time_difference_ratio'])
    refinements.append(row)
for (n,dt),ref in refs.items():
    check(f'FV{n} dt{dt} integrated mass',ref['maximum_integrated_mass_error_m3'],gates['independent_global_mass_error_m3'])
    assert ref['rain_patch_area_m2']==.25 and not ref['reads_native_results']
    assert max(f['transverse_flux_rms_m_s'] for f in ref['frames'])>1e-12
report = dict(schema='soil-fem-local-rain-3d-comparison/1',passed=all(c['passed'] for c in checks),
    comparisons=comparisons,refinements=refinements,checks=checks,solver_sha256=solver_sha,
    predeclared_scope_sha256=sha(p/'scope.json'),comparison_procedure_sha256=sha(Path(__file__)),
    checkpoint_controls_pending=True,mechanical_or_stability_validation=False,p6_complete=False)
with (p/'comparison.json').open('x',encoding='utf-8') as f:
    json.dump(report,f,indent=2,allow_nan=False)
print('Local-rain3D gates:',sum(c['passed'] for c in checks),'/',len(checks),flush=True)
raise SystemExit(0 if report['passed'] else 1)
