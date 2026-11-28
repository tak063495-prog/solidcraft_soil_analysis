"""Read saved new common-lambda diagnostics; never solve or resume a path."""
import argparse
import json
from pathlib import Path
import numpy as np
from fold_stability import point,judge,low_stress,audit_load
from analyze_nc_fold import linear_audit,material_points
from verify_high_order import prepare
from nc_fold_path import CONTROL,digest
import hashlib


def rows(path):return [json.loads(line) for line in path.read_text(encoding='utf-8').splitlines()]
def read(path):return json.loads(path.read_text(encoding='utf-8'))


def analyze(folder):
    hashes=read(folder/'evidence-sha256.json')
    for name,sha in hashes.items():
        if hashlib.sha256((folder/name).read_bytes()).hexdigest()!=sha:raise ValueError('new diagnostic evidence checksum differs')
    controls=read(folder/'controls.json');binding=controls['binding']
    if binding.get('fold_control')!=CONTROL:raise ValueError('not a new common-lambda fold diagnostic')
    g=read(folder/'geometry.json');loads=read(folder/'original-load.json');initial_data=read(folder/'initial.json')
    if initial_data['binding']!=binding or initial_data['payload_sha256']!=digest(initial_data['payload']):raise ValueError('initial checkpoint control or payload changed')
    initial=initial_data['payload']
    samples=rows(folder/'samples.jsonl') if (folder/'samples.jsonl').exists() else []
    systems=rows(folder/'linear-systems.jsonl') if (folder/'linear-systems.jsonl').exists() else []
    local=rows(folder/'gp-tangents.jsonl') if (folder/'gp-tangents.jsonl').exists() else []
    if samples:
        terminal=read(folder/'restart.json')
        if terminal['binding']!=binding or terminal['payload_sha256']!=digest(terminal['payload']) or terminal['payload']!=samples[-1]:raise ValueError('terminal actual checkpoint differs from accepted path')
    if len(local)!=len(samples):raise ValueError('missing accepted GP tangent records')
    q,offset=prepare(g['mesh']);active=np.asarray(g['active_elements']);r=np.asarray(g['transform']);b=np.asarray(g['reduced_b_rows']);w=np.asarray(g['weights_m3'])
    hfull=np.asarray(loads['holding_load_kn']);final=np.asarray(loads['final_load_kn']);h=r.T@hfull
    ids=np.concatenate([np.arange(offset[i],offset[i+1]) for i in np.flatnonzero(active)])
    points=[];stress=None;peak=max((row['progress'] for row in samples),default=None)
    native=g['map']=='native'
    for row,gp in zip(samples,local):
        if row['step']!=gp['step']:raise ValueError('GP tangent step mismatch')
        actual=[eq for eq in systems if eq['step']==row['step']]
        if not actual or actual[-1]['kind']!='tangent':raise ValueError('accepted actual tangent system missing')
        a=np.asarray(actual[-1]['matrix']);k=a[:-1,:-1];plane=a[-1]/1000.
        matrices=np.einsum('qai,qab,qbj,q->qij',b,np.asarray(gp['local_tangents']),b,w,optimize=True)
        if np.linalg.norm(matrices.sum(0)-k)/max(np.linalg.norm(k),1e-30)>1e-10:raise ValueError('GP tangent decomposition differs from actual accepted matrix')
        if not np.allclose(a[:-1,-1],h,rtol=1e-10,atol=1e-10):raise ValueError('actual augmented load column changes original common lambda')
        pieces=[dict(kind='material_gp',point=int(i),matrix=m) for i,m in zip(ids,matrices)]
        pieces.extend(dict(kind='material_gp',point=i,active=False,matrix=np.zeros_like(k)) for i in range(int(offset[-1])) if i not in ids)
        pieces.extend(dict(kind='contact',point=i,active=False,matrix=np.zeros_like(k)) for i in range(len(initial['contacts']['points']) if native else len(initial['contacts']['jump_m'])))
        d=point(k,h,plane,contributions=pieces);audits=linear_audit(actual)
        load=np.asarray(row['applied_nodal_load_kn']);load_audit=audit_load(load,final,hfull,row['progress'])
        states=row['states'];values=np.array([p['effective_stress_compression_positive_kpa'] for p in states])[ids] if native else np.asarray(states)[active].reshape(-1,14)[:,:6]
        force=-np.einsum('qai,qa,q->i',b,values,w);res=float(np.linalg.norm(r.T@load-force)/max(np.linalg.norm(load),1.))
        frozen=all(states[i]==initial['states'][i] for i in range(len(states)) if i not in ids) if native else bool(np.array_equal(np.asarray(states)[~active],np.asarray(initial['states'])[~active]))
        if not frozen or row['contacts']!=initial['contacts']:raise ValueError('removed GP/contact history changed')
        if row['progress']==peak:
            stress=low_stress(material_points(g['mesh'],g['ground'],q,offset,active,states,native,b,np.asarray(d['right_weak_mode'])))
        points.append(dict(step=row['step'],**{'lambda':row['progress']},accepted=True,diagnostic=d,
                           relative_force_residual=res,augmented_linear_residual=max(s['augmented'] for s in audits),
                           force_block_linear_residual=max(s['force_block'] for s in audits),
                           actual_system_audits=audits,full_component_load_audit=load_audit,removed_all_gp_contacts_frozen=frozen))
    return dict(verdict=judge(points),points=points,peak_lambda=peak,passed=False,nc01_complete=False,
                source_binding=binding,interpretation='local static diagnostic; not a dynamic stability or global nonexistence certificate'),stress


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--path',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False);result,stress=analyze(a.path)
    for name,value in [('fold-results.json',result),('low-stress.json',stress)]:
        (a.output/name).write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(verdict=result['verdict'],peak=result['peak_lambda']),indent=2))
