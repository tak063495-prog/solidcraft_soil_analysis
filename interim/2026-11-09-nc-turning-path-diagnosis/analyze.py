"""Reassemble derivative diagnostics from saved states, never certify a solution."""
import argparse
import json
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tools'))
from verify_high_order import prepare

def read(p): return json.loads(p.read_text())
def write(p,v): p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def spectral(k, rhs):
    left,s,right = np.linalg.svd(k)
    return dict(singular_values_kn_m=s.tolist(),condition_number=float(s[0]/s[-1]),
        symmetric_eigenvalues_kn_m=np.linalg.eigvalsh((k+k.T)/2).tolist(),
        asymmetry_relative=float(np.linalg.norm(k-k.T)/np.linalg.norm(k)),
        weak_left_residual_projection_kn=float(left[:,-1]@rhs),
        weak_left_residual_fraction=float(abs(left[:,-1]@rhs)/max(np.linalg.norm(rhs),1e-30)))
def measure(out):
    mesh=read(out/'prior/fresh-head-midpoint/input/mesh.json')
    q,offsets=prepare(mesh)
    records=read(out/'field-records.json')
    points=read(out/'material-results-final.json')['points']
    by={}
    for p in points:
        group=by.setdefault(p['label'],{})
        if p['point'] in group: raise ValueError('duplicate material point')
        group[p['point']]=p
    reports=[]; matrices={}
    for record in records:
        label=record['label'];e=record['event'];ps=by[label]
        expected={int(offsets[ei]+qi) for ei in np.flatnonzero(e['active']) for qi in range(len(q[ei]['weight']))}
        if set(ps)!=expected: raise ValueError('missing or extra active material point')
        errors=[p for p in ps.values() if 'error' in p]
        if errors: raise ValueError(f'{label}: {len(errors)} material errors')
        if not all(p['ad_primal_exact'] and p['substep_replay_exact'] for p in ps.values()):
            raise ValueError('AD/primal or substep replay differs')
        if not all(p['updated_state']==e['states'][i] for i,p in ps.items()):
            raise ValueError('current-head replay differs from saved Newton field')
        rows=e.get('dof_rows') or record['tangent']['dof_rows']
        nd=max(r for row in rows for r,c in row)+1
        transform=np.zeros((len(mesh['nodes'])*3,nd))
        for g,row in enumerate(rows):
            for r,c in row:transform[g,r]=c
        rb={}
        for ei in np.flatnonzero(e['active']):
            eq=q[ei];gd=(eq['nodes'][:,None]*3+np.arange(3)).ravel()
            rb[int(ei)]=eq['b']@transform[gd]
        def assemble(get):
            k=np.zeros((nd,nd))
            for ei,b in rb.items():
                ds=np.array([get(ps[int(offsets[ei]+i)]) for i in range(len(b))])
                k+=np.einsum('qai,qab,qbj,q->ij',b,ds,b,q[ei]['weight'],optimize=True)
            return k
        ad=assemble(lambda p:p['implicit_tangent'])
        rhs=np.array(e['reduced_residual'])
        u,s,vt=np.linalg.svd(ad);v=vt[-1]
        native=None if record['tangent'] is None else np.array(record['tangent']['matrix_columns']).T
        fd=[]
        for j,f in enumerate(next(iter(ps.values()))['finite_differences']):
            k=assemble(lambda p:p['finite_differences'][j]['central'])
            diff=k-ad
            fd.append(dict(width=f['width'],relative_to_ad=float(np.linalg.norm(diff)/np.linalg.norm(ad)),
                weak_mode_action_error_over_ad_sigma_min=float(np.linalg.norm(diff@v)/s[-1]),
                branch_crossing_points=sum(any(a!=b for a,b in p['finite_differences'][j]['plastic_substeps_plus_minus']) for p in ps.values()),
                spectral=spectral(k,rhs),
                relative_to_native=None if native is None else float(np.linalg.norm(k-native)/np.linalg.norm(native))))
            matrices[f'{label}:fd:{f["width"]}']=k.tolist()
        matrices[label+':ad']=ad.tolist()
        step=[s for p in ps.values() for s in p['substeps']]
        reports.append(dict(label=label,progress=e['progress'],committed_progress=e['committed']['stage_progress'],
            iteration=e['iteration'],relative_force_residual=e['relative_force_residual'],
            active_material_points=len(ps),active_contact_points=sum(p['active'] for p in e['joint_report']['points']),
            joint_report=e['joint_report'],plastic_points=sum(p['plastic_substeps']>0 for p in ps.values()),
            material_partition_counts={str(n):sum(e['partitions'][i]==n for i in ps) for n in sorted(set(e['partitions'][i] for i in ps))},
            final_normalized_yield_range=[min(s['normalized_yield_after'] for s in step),max(s['normalized_yield_after'] for s in step)],
            pc_range_kpa=[min(s['pc_kpa'] for s in step),max(s['pc_kpa'] for s in step)],
            pressure_range_kpa=[min(s['p_kpa'] for s in step),max(s['p_kpa'] for s in step)],
            exact_current_head_material_replay=True,implicit_spectral=spectral(ad,rhs),
            native_spectral=None if native is None else spectral(native,rhs),finite_differences=fd))
    summary=dict(schema='soil-fem-nc-fixed-history-derivative-diagnostic/1',passed=False,history_converged=False,
        independent_reference=False,base_commit='06bff71dcc54f7060932aff57194f04ffe0593f8',
        material_requests=len(points),fields=len(records),reports=reports,
        limitations=['saved Newton fields, not an input-only independent solution',
            'finite differences below production width floor are measurements only',
            'no new accepted increments; no full removal or failure/stability verdict'])
    write(out/'derivative-measurements.json',summary);write(out/'assembled-matrices.json',matrices)
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--evidence',type=Path,required=True);a=p.parse_args()
    m=measure(a.evidence.resolve());print('fields',m['fields'],'requests',m['material_requests'],'passed=false')
