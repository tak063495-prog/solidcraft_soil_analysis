"""Reproduce the sealed mean partition's float mismatch without a new NC solve.

Read an already SHA-recorded A-01 controls entry and the original fem_reference
source. This is a small diagnostic, never a substitute for the full strict
archive replay. No Windows float-difference allowance is used.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
LEGACY = ROOT/'interim/2026-11-23-joint-combo-extension/validation/legacy-replay-failures'
CHILD = r'''
import hashlib,importlib.util,json,sys
from pathlib import Path
import numpy as np
root,control,source=map(Path,sys.argv[1:4])
c=json.loads(control.read_text(encoding='utf-8'))
spec=importlib.util.spec_from_file_location('fem_reference',source)
fem=importlib.util.module_from_spec(spec);spec.loader.exec_module(fem);sys.modules['fem_reference']=fem
sys.path.insert(0,str(root/'tools'))
from make_free_contact import make
from verify_high_order import prepare
from bounded_mean_release_control import partition
# These mathematical dependencies are byte-identical to the sealed drivers.
source_checks={}
for name in ('make_free_contact','verify_high_order','high_order_io','mean_stress_release_control','bounded_mean_release_control','spatial_release_control'):
    p=Path(sys.modules[name].__file__);key='tools/'+p.name
    actual=hashlib.sha256(p.read_bytes()).hexdigest();expected=c['binding']['driver_sources'][key]
    if actual!=expected:raise ValueError('minimal diagnostic source differs: '+key)
    source_checks[key]=actual
mesh=make(mcc=True,pc=100.,preload_kpa=5.,force=(1.,.5,0.))[0]
q,_=prepare(mesh);active=np.asarray(c['active_elements']);R=np.asarray(c['transform'])
qs=[q[ei] for ei in np.flatnonzero(active)]
dofs=np.stack([(v['nodes'][:,None]*3+np.arange(3)).ravel() for v in qs])
B=np.einsum('eqai,eij->eqaj',np.stack([v['b'] for v in qs]),R[dofs]).reshape(-1,6,R.shape[1])
w=np.concatenate([v['weight'] for ei,v in enumerate(q) if active[ei]])
first,second=partition(mesh,active,c['split_source_stress_kpa'],c['holding_kn'],c['binding']['order'])
result={}
for name,a,b in [('B',B,np.asarray(c['B'])),('weights',w,np.asarray(c['weights'])),('first',first,np.asarray(c['horizontal_holding_kn'])),('second',second,np.asarray(c['vertical_holding_kn']))]:
    result[name]=dict(exact=bool(np.array_equal(a,b)),changed_components=int(np.count_nonzero(a!=b)),maximum_absolute_difference=float(np.max(np.abs(a-b))))
    if a.ndim==1:result[name]['differences']=[dict(component=int(i),recorded=float(b[i]),remeasured=float(a[i])) for i in np.flatnonzero(a!=b)]
result['geometry_cache_equals_standalone_volume']=all(np.array_equal(q[ei]['b'],fem.volume(mesh,int(ei))['b']) and np.array_equal(q[ei]['weight'],fem.volume(mesh,int(ei))['weight']) for ei in np.flatnonzero(active))
stress=np.asarray(c['split_source_stress_kpa']);padded=np.zeros((len(stress),14));padded[:,:6]=stress
result['source_stress_layout_checks']={}
for name,values in [('C',stress),('state_stride_14',padded[:,:6]),('Fortran',np.asfortranarray(stress))]:
    a,b=partition(mesh,active,values,c['holding_kn'],c['binding']['order'])
    result['source_stress_layout_checks'][name]=dict(strides=list(values.strides),maximum_partition_difference=float(max(np.max(np.abs(a-c['horizontal_holding_kn'])),np.max(np.abs(b-c['vertical_holding_kn'])))))
result['source_sha256']=source_checks
result['jacobi_rules']=[dict(alpha=i,nodes=x.tolist(),weights=w.tolist()) for i in (2,1,0) for x,w in [fem.jacobi(3,i)]]
result['numpy']=np.__version__
result['numpy_build']=np.__config__.show(mode='dicts')
print(json.dumps(result,allow_nan=False))
'''


def run(output, fem_source=None):
    output.mkdir(parents=True, exist_ok=False)
    controls = LEGACY/'a01-mean-controls.json'
    raw = controls.read_bytes()
    control = json.loads(raw)
    manifest = json.loads((ROOT/'interim/2026-11-10-nc-mean-stress-release/archive-sha256.json').read_text())
    expected_control = manifest['entries']['bounded-independent-complement-resumed/controls.json']['sha256']
    if hashlib.sha256(raw).hexdigest() != expected_control:
        raise ValueError('original controls entry SHA differs')
    if fem_source is None:
        source = subprocess.check_output(['git', 'show', '024835c:tools/fem_reference.py'], cwd=ROOT)
    else:
        source = fem_source.read_bytes()
    source_sha = hashlib.sha256(source).hexdigest()
    if source_sha != control['binding']['driver_sources']['tools/fem_reference.py']:
        raise ValueError('original fem_reference source SHA differs')
    source_path = output/'original-fem-reference.py'
    source_path.write_bytes(source)
    cases = [{}, *({'OPENBLAS_CORETYPE': core} for core in ('Haswell', 'Sandybridge', 'Prescott', 'Nehalem', 'Zen')),
             {'NPY_DISABLE_CPU_FEATURES': 'X86_V3'},
             {'OPENBLAS_CORETYPE': 'Prescott', 'NPY_DISABLE_CPU_FEATURES': 'X86_V3'}]
    observations = []
    for extra in cases:
        env = dict(os.environ, OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', OMP_NUM_THREADS='1')
        env.pop('OPENBLAS_CORETYPE', None); env.pop('NPY_DISABLE_CPU_FEATURES', None)
        env.update(extra)
        child = subprocess.run([sys.executable, '-I', '-X', 'utf8', '-c', CHILD, str(ROOT), str(controls), str(source_path.resolve())],
                               env=env, capture_output=True, text=True, check=True)
        observations.append(dict(environment=extra, result=json.loads(child.stdout), stderr=child.stderr))
    result = dict(diagnostic_only=True, new_nc_solve=False, full_strict_archive_replay_executed=False,
                  strict_archive_test_success_claimed=False, platform_float_allowance_used=False,
                  reference_values_or_snapshots_modified=False,
                  controls_archive_entry_sha256=expected_control, original_fem_source_sha256=source_sha,
                  observations=observations)
    (output/'diagnosis.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--fem-source', type=Path, help='optional SHA-verified extracted original source')
    args = parser.parse_args()
    report = run(args.output.resolve(), args.fem_source)
    print(json.dumps([dict(environment=r['environment'], first=r['result']['first']['maximum_absolute_difference'],
                           second=r['result']['second']['maximum_absolute_difference']) for r in report['observations']], indent=2))
