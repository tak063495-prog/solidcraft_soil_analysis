"""All-direction continuum MCC report from an immutable sampled diagnostic."""

def main():
    """Recheck a bound sampled report with the all-direction continuum criterion."""
    from pathlib import Path
    import argparse,hashlib,json,os
    for key in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[key]='1'
    import numpy as np
    from mcc_acoustic_diagnostic import branch,scan_batch,directions
    import mcc_acoustic_diagnostic
    from mcc_acoustic_exact import criterion
    import mcc_acoustic_exact
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('sampled_report',type=Path);parser.add_argument('output',type=Path);args=parser.parse_args()
    assert not args.output.exists();sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    sampled=json.loads(args.sampled_report.read_text(encoding='utf-8-sig'));sources={str(q.resolve()):sha(q) for q in [args.sampled_report,Path(__file__),Path(mcc_acoustic_diagnostic.__file__),Path(mcc_acoustic_exact.__file__)]}
    assert sampled['diagnostic_only'] and sampled['implementation_sha256']==sha(Path(mcc_acoustic_diagnostic.__file__))
    assert len(sampled['inputs_sha256'])==3
    parsed={}
    for name,digest in sampled['inputs_sha256'].items():
     path=Path(name);assert path.stat().st_size<=1024**3
     raw=path.read_bytes();assert hashlib.sha256(raw).hexdigest()==digest
     value=json.loads(raw)
     if 'tetrahedra' in value and 'nodes' in value:kind='mesh'
     elif 'materials' in value:kind='ground'
     elif 'frames' in value and 'integration_points' in value:kind='result'
     else:raise ValueError('unrecognized bound diagnostic input')
     assert kind not in parsed;parsed[kind]=value
    mesh,ground,result=[parsed[k] for k in ['mesh','ground','result']]
    assert result['complete'] and result['frames'][-1]['name']==sampled['frames'][-1]['stage']
    layout=result['integration_points'];owners=np.asarray([q['element_index'] for q in layout]);weights=np.asarray([q['weight_m3'] for q in layout])
    state=result['frames'][-1]['state'];assert len(state['material_states'])==len(layout)
    assert (owners>=0).all() and (owners<len(mesh['tetrahedra'])).all() and np.isfinite(weights).all() and (weights>0).all()
    active=np.asarray(state['active_elements']);assert active.dtype==bool and len(active)==len(mesh['tetrahedra'])
    materials={m['id']:m['mcc'] for m in ground['materials']};regions={r['id']:r['body'] for r in mesh['regions']}
    ids=[];data=[]
    for i in np.flatnonzero(active[owners]):
     s=state['material_states'][i];d=branch(s['effective_stress_compression_positive_kpa'],s['preconsolidation_pressure_kpa'],s['void_ratio'],materials[mesh['material_ids'][owners[i]]])
     if abs(d['normalized_yield'])<=sampled['normalized_yield_tolerance']:ids.append(int(i));data.append(d)
    assert len(ids)==sampled['frames'][-1]['near_yield_candidates']
    scanned=scan_batch(data,directions());exact=[criterion(d) for d in data];groups={};negative=[];sample_negative=[]
    for i,d,s,r in zip(ids,data,scanned,exact,strict=True):
     body=regions[mesh['region_ids'][owners[i]]]
     row=groups.setdefault(body,dict(near_yield_candidates=0,negative_candidates=0,indeterminate_candidates=0,sampled_negative_candidates=0,negative_volume_m3=0.,maximum_rank_one_ratio=0.))
     row['near_yield_candidates']+=1;row['maximum_rank_one_ratio']=max(row['maximum_rank_one_ratio'],r['maximum_rank_one_ratio'])
     if r['decision']=='negative_continuum_candidate':
      row['negative_candidates']+=1;row['negative_volume_m3']+=float(weights[i]);negative.append(i)
     if r['decision']=='near_boundary_indeterminate':row['indeterminate_candidates']+=1
     if s is not None and s['normalized_by_shear'] < -1e-8:
      assert r['decision']=='negative_continuum_candidate';row['sampled_negative_candidates']+=1;sample_negative.append(i)
    last=sampled['frames'][-1];assert len(sample_negative)==last['sampled_negative_acoustic_candidates']
    assert abs(float(weights[sample_negative].sum())-last['sampled_negative_volume_m3'])<1e-12
    maximum=max(range(len(exact)),key=lambda i:exact[i]['maximum_rank_one_ratio']) if exact else None
    witness=None if maximum is None else dict(global_gp=ids[maximum],position_m=layout[ids[maximum]]['position_m'],**exact[maximum])
    assert all(sha(Path(n))==h for n,h in {**sources,**sampled['inputs_sha256']}.items())
    report=dict(diagnostic_only=True,stage=last['stage'],elements=len(mesh['tetrahedra']),near_yield_candidates=len(ids),
                sampled_negative_candidates=len(sample_negative),all_direction_negative_candidates=len(negative),
                additional_candidates_not_detected_by577_directions=len(set(negative)-set(sample_negative)),
                negative_volume_m3=float(weights[negative].sum()),cad_body_groups=groups,maximum_ratio_witness=witness,
                active_loading_branch_not_determined=True,is_not_discrete_mixed_algorithmic_or_global_tangent=True,
                mathematical_maximization_with_floating_point_indeterminate_margin=True,formal_interval_certificate=False,
                spatial_convergence_accepted=False,p6_complete=False,sha256={**sources,**sampled['inputs_sha256']})
    serialized=json.dumps(report,indent=2,allow_nan=False)
    with args.output.open('x',encoding='utf-8') as stream:stream.write(serialized)
    print(json.dumps({k:report[k] for k in ['elements','near_yield_candidates','sampled_negative_candidates','all_direction_negative_candidates','additional_candidates_not_detected_by577_directions','negative_volume_m3']}))

if __name__ == "__main__":
    main()
