"""Check physical force, the MCC ellipse and the stored binary exponents."""
import gzip
import json
from pathlib import Path
import sys
import numpy as np
from nc_restart import Section
from nc_extended import Extended

inputs,binary,extended,route,comparison=map(Path,sys.argv[1:])
base=Section(inputs,binary.resolve())
scaled=Extended(inputs,extended.resolve())
rows=[json.loads(line) for line in (route/'trials.jsonl').open()]
archive=json.loads((route/'summary.json').read_bytes())
state=json.load(gzip.open(route/'last.json.gz','rt'))
s=state['active_material_states']
mant=np.array([r['effective_stress_compression_positive_kpa'] for r in s],dtype=np.longdouble)
exp=np.array([r['stress_exponent'] for r in s])
pc=np.array([r['preconsolidation_pressure_kpa'] for r in s],dtype=np.longdouble)
M=np.array([base.laws[i]['critical_state_m'] for i in base.ids])
assert all(r['void_ratio']>0 for r in s) and (pc>0).all()
stress=mant/pc[:,None]
p=stress[:,:3].mean(1)
dev=stress.copy()
dev[:,:3]-=p[:,None]
f=1.5*((dev[:,:3]**2).sum(1)+2*(dev[:,3:]**2).sum(1))/M**2+p*(p-1)
assert np.all(p>0) and np.all(f<=1e-8)
stress=np.ldexp(mant,exp[:,None]).reshape(48,27,6)
force=np.zeros_like(base.base,dtype=np.longdouble)
for b,w,dofs,sigma in zip(base.b,base.w,base.dofs,stress):
    force[dofs]-=np.einsum('qai,qa,q->i',b,sigma,w)
load=base.load(state['load'])
res=np.linalg.norm((force-load)[base.free])/max(1.,np.linalg.norm(load))
assert res<=base.settings['relative_force_tolerance']
physical_p=np.ldexp(p*pc,exp)
assert physical_p.min()>0
fixed=np.setdiff1d(np.arange(len(force)),base.free)
np.testing.assert_array_equal(np.ravel(state['displacement_m'])[fixed],np.ravel(base.cp['displacement_m'])[fixed])
# Compare a full, nonzero incremental material assembly at the identical P3 parent.
repo=Path(__file__).resolve().parents[2]
start=json.load(gzip.open(repo/archive['source'],'rt'))['active_material_states']
secant=np.load(route.parent/'fine-secant.npy')
du=.0002*secant[:-1]
blocks=[]
for method in [base,scaled]:
    counts=np.ones(len(start),int)
    result=method.evaluate(start,du,counts)
    if result is None:
        result=method.evaluate(start,du,counts)
    blocks.append(result)
a,b=blocks
relative=[np.linalg.norm(x-y)/np.linalg.norm(x) for x,y in zip(a[:2],b[:2])]
assert max(relative)<1e-10 and a[3]==b[3]
other=[json.loads(line) for line in (comparison/'trials.jsonl').open()]
same=[]
for r,t in zip(rows,other):
    if not (r['accepted'] and t['accepted'] and r['ds_m']==t['ds_m']):
        break
    same.append(abs(r['load']-t['load']))
accepted=[r for r in rows if r['accepted']]
report=dict(last_lambda=state['load'],independent_physical_residual=float(res),min_pressure_kpa=str(physical_p.min()),
            min_pc_kpa=str(np.ldexp(pc,exp).min()),min_exponent=int(exp.min()),max_normalized_yield=float(f.max()),
            same_parent_force_tangent_relative=relative,same_parent_plastic_points=a[3],
            common_steps=len(same),max_lambda_difference=max(same,default=None),original_settings=base.settings,
            minimum_abs_projection_kn=min(abs(r['projection_kn']) for r in accepted),
            maximum_true_linear=max(max(v['linear'],v['physical_linear']) for r in accepted for v in r['iterations'] if 'linear' in v))
print(json.dumps(report,indent=2,allow_nan=False))
