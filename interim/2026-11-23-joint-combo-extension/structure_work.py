"""Physical soil+structure endpoint work from both freshly solved histories.

Includes frame force/rotation-moment work. Subtracts the attached structural
load already present in the soil RHS, so it is not counted twice. Production
nodal loads come from its existing local kernels, independent nodal loads from
frame flexibility/plate quadrature. This derived observation is not another
3D solve and is never substituted for one.
"""
import argparse
import itertools
import json
from pathlib import Path
import numpy as np
from fixtures import ROOT,read,write
from reference import Model
from fields import load_fields
from study import invoke,sha


def measure(folder,out,producer,timeout):
    m=Model(*(read(folder/(key+'.json')) for key in ('mesh','ground','case')))
    out.mkdir(exist_ok=False);raw=[]
    if m.frame:
        process=invoke([producer,folder/'mesh.json',folder/'case.json',folder/'native/accepted.jsonl',out/'native-loads.jsonl'],out/'production-loads.log',timeout)
        write(out/'process.json',process);assert process['exit_code']==0 and not process['timed_out']
        raw=[json.loads(line) for line in (out/'native-loads.jsonl').open()]
        bs=m.frame.initial_state()
    previous_n=np.zeros(m.nd);previous_r=np.zeros(m.nd)
    previous_dn=previous_dr=None;wn=wr=0.;maximum_load=0.;rows=[]
    for i,(left,right) in enumerate(itertools.zip_longest(load_fields(folder/'native-fields.jsonl'),load_fields(folder/'independent-fields.jsonl'))):
        assert left is not None and right is not None
        a,n=left;b,r=right;assert (a['stage'],a['progress'])==(b['stage'],b['progress'])
        un=n['u'].ravel();ur=r['u'].ravel();fn=n['applied_load'].ravel().copy();fr=r['applied_load'].ravel().copy()
        dn=dr=ln=lr=np.array([]);attached_n=np.zeros(m.nd);attached_r=np.zeros(m.nd)
        if m.frame:
            assert raw[i]['factor']==a['factor']==b['factor']
            ln=np.asarray(raw[i]['structure_nodal_force_moment']);attached_n=np.asarray(raw[i]['attached_nodal_load_kn'])
            _,_,lr,eff,*_=m.frame.assemble(bs,b['factor'])
            for node,soil in enumerate(m.frame.attach):
                if soil is not None:attached_r[3*soil:3*soil+3]+=eff[node]
            dn=np.concatenate((n['structure_u'],n['structure_rotation']),axis=1).ravel()
            dr=np.concatenate((r['structure_u'],r['structure_rotation']),axis=1).ravel()
            if previous_dn is None:previous_dn=np.zeros_like(dn);previous_dr=np.zeros_like(dr)
            maximum_load=max(maximum_load,float(np.max(np.abs(ln-lr),initial=0)),float(np.max(np.abs(attached_n-attached_r),initial=0)))
            fn-=attached_n;fr-=attached_r
        soil_n=float(fn@(un-previous_n));soil_r=float(fr@(ur-previous_r))
        structure_n=float(ln@(dn-previous_dn)) if m.frame else 0.
        structure_r=float(lr@(dr-previous_dr)) if m.frame else 0.
        wn+=soil_n+structure_n;wr+=soil_r+structure_r
        rows.append(dict(stage=a['stage'],progress=a['progress'],native=dict(soil_own_load_kn=fn.tolist(),structure_force_moment=ln.tolist(),attached_load_subtracted_kn=attached_n.tolist(),
            soil_endpoint_work_kn_m=soil_n,structure_endpoint_work_kn_m=structure_n,total_endpoint_work_cumulative_kn_m=wn),
            independent=dict(soil_own_load_kn=fr.tolist(),structure_force_moment=np.asarray(lr).tolist(),attached_load_subtracted_kn=attached_r.tolist(),
            soil_endpoint_work_kn_m=soil_r,structure_endpoint_work_kn_m=structure_r,total_endpoint_work_cumulative_kn_m=wr)))
        previous_n=un;previous_r=ur;previous_dn=dn;previous_dr=dr
    expected=read(folder/'independent-result.json')['accepted_states']
    assert len(rows)==expected and (not m.frame or len(raw)==len(rows))
    assert maximum_load<1e-10,'existing all-node load gate'
    difference=np.array([x['native']['total_endpoint_work_cumulative_kn_m']-x['independent']['total_endpoint_work_cumulative_kn_m'] for x in rows])
    reference=np.array([x['independent']['total_endpoint_work_cumulative_kn_m'] for x in rows])
    result=dict(accepted_pairs=len(rows),derived_from_selected_fresh_solutions=True,counted_as_additional_3d_solve=False,
        native_local_kernel_sha256=sha(producer),maximum_structure_and_attached_nodal_load_difference=maximum_load,
        maximum_total_cumulative_work_difference_kn_m=float(np.max(np.abs(difference))),
        total_cumulative_work_relative_l2=float(np.linalg.norm(difference)/max(np.linalg.norm(reference),1e-12)),
        physical_structure_force_and_moment_work_included=True,attached_load_not_double_counted=True,numerical_acceptance=False,rows=rows)
    write(out/'verification.json',result);return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--matrix',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--producer',type=Path,default=ROOT/'target/release/examples/joint_combo_structure_loads');p.add_argument('--timeout',type=int,default=900);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False)
    result={folder.name:measure(folder,a.output/folder.name,a.producer.resolve(),a.timeout) for folder in sorted(a.matrix.iterdir()) if (folder/'comparison.json').is_file()}
    write(a.output/'measurements.json',result)


if __name__=='__main__':main()
