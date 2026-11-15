"""Unmodified element observations and bounded assembled global direct solves."""
import argparse
import json
import math
from pathlib import Path
import shutil
import subprocess
import time
import numpy as np
from b01follow_common import BASE, HERE, ROOT, SECTION, dot, metric, mm, mv, norm, rigid, rotate, sha, write
from b01follow_plan import FIELDS, SHAPES, patch, rejection_controls, rotations, specs, target_fields


def query(probe, output, label, requests):
    write(output/(label+'-requests.json'), requests)
    raw = (output/(label+'-requests.json')).read_bytes(); started = time.monotonic()
    p = subprocess.Popen([str(probe)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdout,stderr = p.communicate(raw,timeout=300)
    (output/(label+'-responses.json')).write_bytes(stdout); (output/(label+'-stderr.log')).write_bytes(stderr)
    write(output/(label+'-process.json'),dict(pid=p.pid,exit_code=p.returncode,elapsed_seconds=time.monotonic()-started,
          input_sha256=sha(output/(label+'-requests.json')),binary_sha256=sha(probe)))
    if p.returncode: raise subprocess.CalledProcessError(p.returncode,[str(probe)])
    results = json.loads(stdout)
    if [r['id'] for r in requests] != [r['id'] for r in results]: raise ValueError('request/result inventory differs')
    return results


def spectrum(k,x):
    values,vectors = np.linalg.eigh((k+k.T)*.5)
    threshold = 1e-10*float(values[-1]); null = vectors[:,values <= threshold]; r = rigid(x)
    outside = r-mm(null,mm(null.T,r))
    residual = mm(k,r)
    return dict(eigenvalues=values.tolist(),inherited_threshold=threshold,
                modes_below_inherited_threshold=int(np.sum(values <= threshold)),
                negative_below_minus_threshold=int(np.sum(values < -threshold)),
                physical_rigid_basis=r.tolist(),algebraic_null_basis=null.tolist(),
                physical_rigid_residual=residual.tolist(),
                physical_rigid_normalized_residual=[norm(residual[:,i])/(norm(k)*norm(r[:,i])) for i in range(6)],
                physical_rigid_fraction_outside_nullspace=[norm(outside[:,i])/norm(r[:,i]) for i in range(6)],
                physical_rigid_translation_force_l2_kn=[norm(residual[:,i].reshape(-1,6)[:,:3]) for i in range(6)],
                physical_rigid_rotation_moment_l2_kn_m=[norm(residual[:,i].reshape(-1,6)[:,3:]) for i in range(6)],
                physical_rigid_energy_kn_m=[.5*dot(r[:,i],residual[:,i]) for i in range(6)])


def solve(k,x,boundary,targets,q):
    bd = [6*i+a for i in boundary for a in range(6)]; free = sorted(set(range(len(k)))-set(bd))
    reduced = k[np.ix_(free,free)]; rows = []
    for name in FIELDS:
        d = np.zeros(len(k)) if name == 'load' else targets[name].copy()
        f = np.zeros(len(k))
        if name == 'load':
            ids = sorted(set(i//6 for i in free))
            load = rotate(np.array([.3,-.2,1.,.01,-.02,.015]),q)/len(ids)
            for i in ids: f[6*i:6*i+6] = load
        rhs = f[free]-mv(k[np.ix_(free,bd)],d[bd]); started = time.monotonic()
        d[free] = np.linalg.solve(reduced,rhs)
        internal = mv(k,d); residual = internal-f; reaction = np.zeros(len(k)); reaction[bd] = residual[bd]
        z = residual[free].reshape(-1,6)
        row = dict(name=name,displacement_rotation=d.tolist(),applied_force_moment=f.tolist(),
                   internal_force_moment=internal.tolist(),residual=residual.tolist(),reaction=reaction.tolist(),
                   free_translation_residual_l2_kn=norm(z[:,:3]),free_rotation_residual_l2_kn_m=norm(z[:,3:]),
                   mixed_dof_algebraic_relative_residual=norm(residual[free])/max(norm(internal),norm(f),1e-12),
                   mixed_dof_residual_is_algebraic_only=True,
                   strain_energy_kn_m=.5*dot(d,internal),applied_work_kn_m=dot(d,f),reaction_work_kn_m=dot(d,reaction),
                   elapsed_seconds=time.monotonic()-started)
        if name != 'load':
            error = (d-targets[name]).reshape(-1,6)
            row.update(target=targets[name].tolist(),translation_error=metric(d.reshape(-1,6)[:,:3],targets[name].reshape(-1,6)[:,:3]),
                       rotation_error=metric(d.reshape(-1,6)[:,3:],targets[name].reshape(-1,6)[:,3:]),
                       prescribed_field_is_exact_patch=not np.any(x[:,2]) or name.startswith('rigid'))
        rows.append(row)
    return bd,free,rows,float(np.linalg.cond(reduced))


def run(output,probe):
    output.mkdir(parents=True,exist_ok=False); write(output/'plan.json',specs())
    for source in sorted(HERE.iterdir()):
        if source.suffix in ('.py','.rs') or source.name == 'SPEC.ja.md':
            target = output/'analysis-source-snapshot'/source.name
            target.parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(source,target)
    write(output/'analysis-source-binding.json',{p.name:sha(p) for p in sorted((output/'analysis-source-snapshot').iterdir())})
    sources = sorted(p for p in ROOT.glob('src/**/*.rs'))
    source_binding = {p.relative_to(ROOT).as_posix():sha(p) for p in sources}
    write(output/'source-binding-before.json',source_binding)
    requests=[]; meshes=[]
    for shape in SHAPES:
        original,elements,boundary = patch(shape)
        for rotation,q in rotations():
            name = shape+'--'+rotation; x = rotate(original,q); ids=[]
            for i,nodes in enumerate(elements):
                rid = name+'--e'+str(i); ids.append(rid)
                requests.append(dict(id=rid,points=x[nodes].tolist(),section=SECTION,fields=[]))
            meshes.append(dict(id=name,shape=shape,rotation=rotation,q=q,original=original,x=x,elements=elements,boundary=boundary,ids=ids))
    for name,x in rejection_controls():
        for rotation,q in rotations():
            requests.append(dict(id='reject--'+name+'--'+rotation,points=rotate(x,q).tolist(),section=SECTION,fields=[]))
    results=query(probe,output,'kernel',requests); byid={r['id']:r for r in results}
    recovery_requests=[]; inventory=[]
    for meta in meshes:
        x,elements,q=meta['x'],meta['elements'],meta['q']; size=6*len(x)
        terms={}; k=np.zeros((size,size))
        for rid,nodes in zip(meta['ids'],elements):
            record=byid[rid]
            if not record['accepted']: raise ValueError('permitted patch rejected: '+rid)
            ke=np.asarray(record['stiffness']); dofs=[6*i+a for i in nodes for a in range(6)]
            for a,i in enumerate(dofs):
                for b,j in enumerate(dofs): terms.setdefault((i,j),[]).append(float(ke[a,b]))
        for (i,j),v in sorted(terms.items()): k[i,j]=math.fsum(v)
        bd,free,solves,condition=solve(k,x,meta['boundary'],target_fields(meta['original'],q),q)
        rec=dict(id=meta['id'],shape=meta['shape'],rotation=meta['rotation'],rotation_matrix=q.tolist(),
                 points=x.tolist(),elements=elements,boundary_nodes=meta['boundary'],boundary_dofs=bd,free_dofs=free,
                 stiffness=k.tolist(),spectrum=spectrum(k,x),reduced_mixed_dof_condition_number=condition,
                 condition_times_epsilon=condition*np.finfo(float).eps,solves=solves)
        write(output/'assemblies'/(meta['id']+'.json'),rec)
        inventory.append(dict(id=meta['id'],shape=meta['shape'],rotation=meta['rotation']))
        for rid,nodes in zip(meta['ids'],elements):
            fields=[]; dofs=[6*i+a for i in nodes for a in range(6)]
            for row in solves:
                fields.append(dict(name=row['name'],d=np.array(row['displacement_rotation'])[dofs].tolist()))
            for i in range(6): fields.append(dict(name='physical-rigid-'+str(i),d=rigid(x)[:,i][dofs].tolist()))
            recovery_requests.append(dict(id=rid,points=x[nodes].tolist(),section=SECTION,fields=fields))
        print(meta['id'],len(solves),'solves',flush=True)
    query(probe,output,'recovery',recovery_requests)
    repeat=query(probe,output,'kernel-repeat',requests)
    write(output/'same-input-repeat.json',dict(separate_processes=True,exact_raw_bytes_equal=
          (output/'kernel-responses.json').read_bytes()==(output/'kernel-repeat-responses.json').read_bytes(),
          semantic_equal=repeat==results))
    after={p.relative_to(ROOT).as_posix():sha(p) for p in sources}; write(output/'source-binding-after.json',after)
    if source_binding!=after:raise ValueError('production sources changed during observation')
    write(output/'solve-ledger.json',dict(base=BASE,diagnostic_only=True,passed=False,b01_complete=False,
          patches=inventory,new_global_solves=192,permitted_element_requests=256,rejected_controls=16,
          same_input_repeat_is_new_physics=False,recovery_is_new_global_solve=False,production_sources_unchanged=True))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--probe',type=Path,required=True)
    a=p.parse_args();run(a.output.resolve(),a.probe.resolve())
