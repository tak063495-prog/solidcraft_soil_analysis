"""Same captured linear systems: SVD/range tests and direct solves, no continuation."""
import argparse
import json
from pathlib import Path
import numpy as np
from a01cause_metrics import norm,dot
from a01cause_driver import write


def matrix(event):
    if 'matrix' in event:return np.asarray(event['matrix'])
    storage=event['storage']
    if storage['kind']=='dense':return np.asarray(storage['matrix_columns']).T
    n=len(event['rhs']);a=np.zeros((n,n))
    for i in range(n):
        for j in range(storage['row_ptr'][i],storage['row_ptr'][i+1]):a[i,storage['columns'][j]]=storage['values'][j]
    return a


def assess(a,b):
    a=np.asarray(a);b=np.asarray(b);u,s,vh=np.linalg.svd(a,full_matrices=False)
    cutoff=max(a.shape)*np.finfo(float).eps*s[0];rank=int(np.sum(s>cutoff))
    result=dict(dimension=len(b),sigma_max=float(s[0]),sigma_min=float(s[-1]),
                condition_2=float(s[0]/s[-1]) if s[-1] else None,numerical_rank=rank,
                rank_classification_cutoff=float(cutoff),rank_cutoff_is_acceptance_tolerance=False,
                asymmetry_relative=norm(a-a.T)/max(norm(a),1e-300),
                rhs_outside_numerical_range_relative=norm(u[:,rank:].T@b)/max(norm(b),1e-300),
                smallest_singular_values=s[-min(8,len(s)):].tolist())
    eigen=np.linalg.eigvalsh(.5*(a+a.T));result['symmetric_part_eigen_min']=float(eigen[0]);result['symmetric_part_eigen_max']=float(eigen[-1])
    try:
        x=np.linalg.solve(a,b);result['direct_lu']=dict(true_relative_residual=norm(a@x-b)/max(norm(b),1e-300),solution_norm=norm(x),solution=x.tolist())
    except np.linalg.LinAlgError as e:result['direct_lu']=dict(error=str(e))
    least,_,_,_=np.linalg.lstsq(a,b,rcond=max(a.shape)*np.finfo(float).eps)
    result['svd_range_solve']=dict(true_relative_residual=norm(a@least-b)/max(norm(b),1e-300),solution_norm=norm(least))
    scale=np.max(np.abs(a),axis=1);scale=np.where(scale==0,1.,scale)
    scaled=a/scale[:,None];ss=np.linalg.svd(scaled,compute_uv=False)
    result['row_scaled_condition_2']=float(ss[0]/ss[-1]) if ss[-1] else None
    result['null_right_vectors']=vh[-min(6,len(s)):].tolist()
    return result


def native_events(path):
    context={};events=[];last_tangent=None
    for line in path.read_text(encoding='utf-8').splitlines():
        event=json.loads(line);kind=event['kind']
        if kind=='iterate':context={k:event[k] for k in ('stage','progress','iteration','relative_force_residual')}
        if kind=='tangent':last_tangent=event
        if kind=='a01_linear_system':events.append(dict(id=len(events),context=context.copy(),**event,original_tangent=last_tangent))
        if kind=='a01_linear_result':events[-1]['native_result']=event
    return events


def rigid_fit(mesh,transform,vectors):
    xyz=np.asarray(mesh['nodes']);pending=set(range(len(xyz)));adj={i:set() for i in pending};groups=[]
    for cell in mesh['tetrahedra']:
        for i in cell:adj[i].update(cell)
    while pending:
        stack=[min(pending)];nodes=[]
        while stack:
            i=stack.pop()
            if i not in pending:continue
            pending.remove(i);nodes.append(i);stack.extend(sorted(adj[i],reverse=True))
        groups.append(sorted(nodes))
    records=[]
    for vector in vectors:
        displacement=(transform@np.asarray(vector)).reshape(-1,3);blocks=[]
        for nodes in groups:
            x=xyz[nodes]-xyz[nodes].mean(axis=0)
            basis=np.vstack([np.c_[np.eye(3),np.cross(np.eye(3),r).T] for r in x])
            value=displacement[nodes].ravel();coeff=np.linalg.lstsq(basis,value,rcond=None)[0]
            blocks.append(dict(z_min=float(xyz[nodes,2].min()),z_max=float(xyz[nodes,2].max()),
                motion_norm=norm(value),rigid_fit_relative_residual=norm(basis@coeff-value)/max(norm(value),1e-30),translation_rotation=coeff.tolist()))
        records.append(blocks)
    return records


def main():
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False)
    events=native_events(a.input/'native/trace.jsonl');native=[];selected=[]
    mesh=json.loads((a.input/'mesh.json').read_text())
    for e in events:
        assessment=assess(matrix(e),e['rhs'])
        original=e['original_tangent'];a0=np.asarray(original['matrix_columns']).T
        # These are distinct if existing cold/gauge predictor modifies A,b.
        assessment['pre_predictor_tangent']=assess(a0,original['reduced_residual'])
        transform=np.zeros((len(original['dof_rows']),len(e['rhs'])))
        for i,row in enumerate(original['dof_rows']):
            for j,coefficient in row:transform[i,j]=coefficient
        assessment['actual_small_mode_rigid_fits']=rigid_fit(mesh,transform,assessment['null_right_vectors'][-2:])
        assessment['matrix_change_from_pre_predictor_max']=float(np.max(np.abs(matrix(e)-a0)))
        native.append(dict(id=e['id'],context=e['context'],native_result=e['native_result'],assessment=assessment))
        if 'error' in e['native_result'] or e['id']>=len(events)-3:selected.append(e)
    write(a.output/'native-spectra.json',native);write(a.output/'replay-input.json',selected)
    independent=[]
    for line in (a.input/'independent/trace.jsonl').read_text(encoding='utf-8').splitlines():
        e=json.loads(line)
        if e['kind']=='linear':
            independent.append(dict(stage=e['stage'],progress=e['progress'],iteration=e['iteration'],relative=e['relative'],
                gap=e['gap'],traction=e['traction'],contact_partitions=e['contact_partitions'],
                assessment=assess(e['matrix'],e['rhs'])))
    write(a.output/'independent-spectra.json',independent)
    print('native',len(native),'independent',len(independent))
    print('native last',native[-1]['context'],{k:v for k,v in native[-1]['assessment'].items() if k not in ('null_right_vectors','pre_predictor_tangent','direct_lu')},native[-1]['assessment']['direct_lu'].get('true_relative_residual'))
    print('reference last', {k:v for k,v in independent[-1]['assessment'].items() if k not in ('null_right_vectors','direct_lu')},independent[-1]['assessment']['direct_lu'].get('true_relative_residual'))

if __name__=='__main__':main()
