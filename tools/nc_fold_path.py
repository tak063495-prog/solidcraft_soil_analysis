"""Bound the existing moving-plane diagnostic to the original common-lambda load.

This adds observation and checkpoint control binding. The equations, Newton
iterations, material maps, tolerances and step sizes are those of continue_path.
It cannot certify NC-01, physical stability, or a completed production job.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import platform
import scipy
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
PRIOR=ROOT/'interim/2026-11-09-nc-pseudo-arclength'
CONTROL='common_lambda_moving_tangent_plane/1'


def load_driver():
    sys.path.insert(0,str(PRIOR))
    p=ROOT/'interim/2026-11-09-nc-beyond-fold/continue_path.py'
    spec=importlib.util.spec_from_file_location('nc_fold_existing_path',p)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,allow_nan=False).encode()).hexdigest()


def control_binding(base):
    return dict(base,fold_control=CONTROL,arc_lambda_scale_m=1.,force_tolerance=1e-8,
                true_linear_tolerance=1e-10,arc_constraint_tolerance_m=1e-12,
                checkpoint_state_fields=['progress','tangent','states','displacement_m','contacts','arc_length_m','applied_nodal_load_kn'])


def native_loads(mesh,ground,case,origin,active):
    """Independent reconstruction of the original removal bookkeeping.

    At restrained DOFs, retained internal force is NOT the prescribed holding
    load. Preserve the existing external-minus-removed-minus-contact policy.
    """
    from verify_high_order import prepare,physical,internal
    from joint_reference import System
    import fem_reference as fem
    q,offsets=prepare(mesh);active=np.asarray(active,dtype=bool)
    joint=System(mesh,case['joint_system'])
    removed=np.asarray(origin['active_elements'])&~active
    removed_force,_=internal(mesh,ground,origin['material_states'],removed,q,offsets)
    final=physical(mesh,ground,active,0.,True,joint_system=joint)[0]
    retained=dict(origin,active_elements=active.tolist())
    holding=(np.asarray(origin['applied_nodal_load_kn']).reshape(-1,3)-final-removed_force
             -joint.audit(origin)+joint.audit(retained))
    used=np.unique(fem.cells(mesh)[active])
    holding[np.setdiff1d(np.arange(len(mesh['nodes'])),used)]=0.
    return final.ravel(),holding.ravel()


def reconstructed_loads(model):
    from verify_high_order import physical
    m=model.m
    active=np.asarray(m.active if model.kind=='independent' else m.cp['active_elements'],dtype=bool)
    if model.kind=='independent':
        final=physical(m.mesh,m.ground,active,0.,True,joint_system=m.system.joint)[0].ravel()
        unload=json.loads((Path(model.reference)/'unload.json').read_text(encoding='utf-8'))
        origin=m.system.total_force(np.asarray(unload['material_state'])[active].reshape(-1,14),
                                   np.asarray(unload['displacement_m']).ravel())
        holding=origin-final
    else:
        origin_state=m.stage_origin if hasattr(m,'stage_origin') else m.checkpoint['run']['continuation']['origin']
        final,holding=native_loads(m.mesh,m.ground,m.case,origin_state,active)
    for lam in (0.,model.lam,1.):
        from fold_stability import audit_load
        audit_load(m.load(lam),final,holding,lam)
    return final,holding


def validate_checkpoint(path,binding):
    """Accept only a checkpoint actually emitted with its unchanged evidence.

    Payload checksums alone cannot distinguish a caller re-sealing a changed
    branch vector. Match the original bytes and the actual accepted row too.
    """
    data=json.loads(path.read_text(encoding='utf-8'))
    if data['binding']!=binding or data['payload_sha256']!=digest(data['payload']):
        raise ValueError('new fold checkpoint binding or payload differs')
    if binding.get('fold_control')!=CONTROL:raise ValueError('old checkpoint has no common-lambda control binding')
    folder=path.parent;manifest=json.loads((folder/'evidence-sha256.json').read_text(encoding='utf-8'))
    if hashlib.sha256(path.read_bytes()).hexdigest()!=manifest.get('restart.json'):
        raise ValueError('checkpoint is not the original emitted checkpoint bytes')
    for name,sha in manifest.items():
        if hashlib.sha256((folder/name).read_bytes()).hexdigest()!=sha:
            raise ValueError('checkpoint sidecar evidence changed')
    matches=[json.loads(line) for line in (folder/'samples.jsonl').read_text(encoding='utf-8').splitlines()
             if json.loads(line)['step']==data['payload']['step']]
    if matches!=[data['payload']]:raise ValueError('checkpoint does not match the actual accepted full field')
    return data['payload']


def run(args):
    module=load_driver();module.validate_controls(args.step_size,args.count,args.budget_seconds)
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    if args.map=='independent':
        from independent_map_arc import read_arc_seed,seed_failure
        try:read_arc_seed(args.arc)
        except ValueError as exc:return seed_failure(out,exc)
    if getattr(args,'fresh_native',None):
        if args.map!='native':raise ValueError('fresh native history cannot seed independent reference')
        from nc_fold_fresh import Adapter
        model=Adapter(args.fresh_native,args.binary,out)
    else:
        model=module.Adapter(args.map,args.evidence,args.binary,args.reference,args.arc,out,args.material_steps)
    model.reference=args.reference
    final,holding=reconstructed_loads(model)
    source_paths=[Path(__file__),ROOT/'tools/nc_fold_fresh.py',ROOT/'tools/nc_fixed_inputs.py',ROOT/'tools/fold_stability.py',ROOT/'interim/2026-11-09-nc-beyond-fold/continue_path.py',
                  PRIOR/'native_arc.py',PRIOR/'independent_map_arc.py',PRIOR/'independent_branches.py',PRIOR/'independent_arc.py',
                  ROOT/'tools/mcc_rate_fem.py',ROOT/'tools/mcc_rate_rk4.py',ROOT/'tools/mcc_rate_midpoint.py',
                  ROOT/'tools/fem_reference.py',ROOT/'tools/verify_high_order.py',ROOT/'tools/joint_reference.py',
                  ROOT/'tools/low_preload_nc_reference.py',ROOT/'tools/probe_nc_midpoint.py',ROOT/'tools/probe_low_preload_nc.py',ROOT/'tools/make_free_contact.py']
    source_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}
    binding=control_binding(dict(map=args.map,source_sha256=source_hashes,model_sha256=digest(model.source),
                                geometry_sha256=digest(dict(transform=model.R.tolist(),reduced_b=(model.m.B if args.map=='independent' else model.m.rb).tolist(),weights=model.m.w.tolist())),
                                runtime=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,platform=sys.platform,machine=platform.machine()),
                                step_size_m=args.step_size,material_steps=args.material_steps,
                                original_initial_increment=model.m.case['stages'][2].get('initial_increment'),
                                independent_reference_steps=(json.loads((args.reference/'reference.json').read_text(encoding='utf-8'))['steps'] if args.map=='independent' else None),
                                final_load_sha256=digest(final.tolist()),holding_sha256=digest(holding.tolist())))
    if args.resume is not None:validate_checkpoint(args.resume,binding)
    module.write(out/'controls.json',dict(binding=binding,diagnostic_only=True,production_adopted=False,
                                         count=args.count,budget_seconds=args.budget_seconds))
    module.write(out/'original-load.json',dict(final_load_kn=final.tolist(),holding_load_kn=holding.tolist(),
                                               reconstruction=('original external load minus removed soil force minus lost contact force minus final physical load; inactive nodes zero' if args.map=='native' else 'own original pre-cut force minus physical final load')))
    m=model.m
    module.write(out/'geometry.json',dict(map=args.map,mesh=m.mesh,ground=m.ground,case=m.case,
        active_elements=np.asarray(m.active if args.map=='independent' else m.cp['active_elements']).tolist(),
        transform=model.R.tolist(),reduced_b_rows=(m.B if args.map=='independent' else m.rb).tolist(),
        weights_m3=(m.w).tolist()))
    original_load=m.load
    def common_load(lam):
        value=original_load(lam)
        from fold_stability import audit_load
        audit_load(value,final,holding,lam)
        return value
    m.load=common_load

    # Observe existing evaluations. No extra constitutive integration or solver call.
    original_evaluate=model.evaluate
    latest={}
    if args.map=='independent':
        original_integrate=model.m.material.integrate
        integrations=[]
        def integrate(*a,**kw):
            value=original_integrate(*a,**kw);integrations.append(value.copy());return value
        model.m.material.integrate=integrate
    def evaluate(du,old,label):
        if args.map=='independent':integrations.clear()
        result=original_evaluate(du,old,label)
        if args.map=='independent':
            if len(integrations)!=13:raise ValueError('existing independent material evaluation layout changed')
            local=np.stack([(integrations[1+2*j][:,:6]-integrations[2+2*j][:,:6])/(2e-8) for j in range(6)],axis=2)
        else:
            file=model.m.out/f'call-{model.m.call-1:05d}-result.json'
            response=json.loads(file.read_text(encoding='utf-8'))
            local=np.array([p['implicit_tangent'] for p in response['points']])
        latest.update(local_tangents=local.tolist(),reduced_increment_m=np.asarray(du).tolist())
        return result
    model.evaluate=evaluate
    original_append=module.append
    def append(path,row):
        if path.name=='samples.jsonl':row['applied_nodal_load_kn']=model.m.load(row['progress']).tolist()
        original_append(path,row)
        if path.name=='samples.jsonl':
            original_append(out/'gp-tangents.jsonl',dict(step=row['step'],**latest))
    module.append=append
    try:
        report=module.follow(model,out,args.step_size,args.count,args.budget_seconds,binding,args.resume)
    finally:
        module.append=original_append
    hashes={str(p.relative_to(out)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.rglob('*')) if p.is_file()}
    module.write(out/'evidence-sha256.json',hashes)
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--map',choices=['native','independent'],required=True)
    for name in ('evidence','binary','reference','arc','resume'):parser.add_argument('--'+name,type=Path)
    parser.add_argument('--fresh-native',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--step-size',type=float,required=True)
    parser.add_argument('--count',type=int,required=True)
    parser.add_argument('--material-steps',type=int,choices=[4,8],default=4)
    parser.add_argument('--budget-seconds',type=float,default=60.)
    result=run(parser.parse_args())
    print(json.dumps(result,indent=2))
    sys.exit(1 if result.get('diagnostic_status')=='failed' or result['equilibrium_samples']==0 else 0)
