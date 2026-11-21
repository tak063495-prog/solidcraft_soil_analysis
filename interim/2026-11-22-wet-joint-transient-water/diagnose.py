"""Saved-field cause diagnostics, explicitly separate from fresh solves."""
import argparse
import importlib.util
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('c04_diagnose_series',HERE/'saturated_control.py')
ss=importlib.util.module_from_spec(spec);spec.loader.exec_module(ss)
s=ss.s

def diagnose(out):
    common=sorted(t for t,_ in ss.measure.hydraulic_fields(out,out/'saturated-controls/h4','native'))
    shared={}
    for n in (4,8,16):
     case=s.read(out/f'saturated-controls/h{n}/case.json');series=ss.Series(case);run=s.read(out/f'saturated-controls/h{n}/native.json');store={}
     for frame in run['frames'][1:]:
      t=frame['state']['time_s']
      if t not in common:continue
      z=np.array([p['position_m'][2] for p in frame['points']]);head=series.fields(z,t)[0]
      ss.measure.add(store,'hydraulic_gp_head_m',np.array([p['total_head_m'] for p in frame['points']]),head)
     shared[str(n)]=ss.measure.finish(store)
    rows=[]
    for n in (8,16,32):
     a=list(ss.measure.fields(out/f'native-w{n}'));b=list(ss.measure.fields(out/f'independent-w{n}'));prev_a=np.zeros((6,3));prev_b=prev_a.copy()
     for (ra,fa),(rb,fb) in zip(a,b):
      mask=fa['contact_status']!=fb['contact_status']
      for gp in np.flatnonzero(mask):
       rows.append(dict(water_steps=n,stage=ra['stage'],progress=ra['progress'],point=int(gp),native_status=int(fa['contact_status'][gp]),independent_status=int(fb['contact_status'][gp]),
        gap_m=float(fa['gap'][gp]),native_plastic_increment_m=float(np.linalg.norm(fa['plastic_slip'][gp]-prev_a[gp])),independent_plastic_increment_m=float(np.linalg.norm(fb['plastic_slip'][gp]-prev_b[gp])),
        native_incremental_slip_tolerance_applied=False,traction_difference_kpa=float(np.max(abs(fa['traction'][gp]-fb['traction'][gp])))))
      prev_a=fa['plastic_slip'];prev_b=fb['plastic_slip']
    indparts={n:dict(sorted(__import__('collections').Counter(x for r in s.read(out/f'independent-w{n}/result.json')['rows'] for x in r['contact_partitions']).items())) for n in (8,16,32)}
    diag=dict(C04='未達',passed=False,saved_output_diagnostic_not_new_solve=True,contact_status_mismatch_details=rows,independent_contact_partition_histograms=indparts,
     maximum_status_mismatch_native_plastic_increment_m=max(x['native_plastic_increment_m'] for x in rows),maximum_status_mismatch_independent_plastic_increment_m=max(x['independent_plastic_increment_m'] for x in rows),
     continuum_series_error_at_common_28_times=shared,pressure_transfer_error_separated=True,
     cause_scope='Same-old-state/same-partition local replay matches; fresh global solvers choose different contact partitions. No full global partition-policy swap experiment was performed; no unique general cause or continuum acceptance is claimed.')
    return diag

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--report',type=Path,required=True)
    a=p.parse_args();s.write(a.report,diagnose(a.output.resolve()))
