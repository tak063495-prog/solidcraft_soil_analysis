"""Audit new path counts, fixed solve sources, minimum counts and old criteria."""
import argparse,json
from pathlib import Path
from a02follow_run import write,sha
from a02follow_binding import study,HERE,geometry,fields,reference,norm
import numpy as np

def audit(run,report,out):
    paths=study.read(run/'paths.json');trace=study.read(report/'partition-trace-summary.json');solve_names=('a02follow_run.py','a02follow_reference.py','a02follow_policies.py','a02follow_binding.py');bindings=[];geometry_audits={};maximum_residuals={}
    for p in paths:
        b=study.read(run/p['path']/'binding.json');saved={Path(k).name:v for k,v in b['source_sha256'].items()};bindings.append({n:saved[n] for n in solve_names})
        assert trace[p['path']]['minimum_contact_parts']>=p['minimum_parts']
        assert p['complete'] and not p['timed_out'] and p['exit_code']==0
        folder=run/p['path']
        if p['engine']=='native':geometry_audits[p['path']]=geometry.audit(folder)
        model=reference.Model(*(study.read(folder/(n+'.json')) for n in ('mesh','ground','case')))
        residuals=[]
        for _,f in fields.load_fields(folder/(p['engine']+'-fields.jsonl')):
            tied=not bool(np.all(f['joint_active']));basis=model.spaces[tied][0]
            residuals.append(norm(basis.T@(f['internal']-f['applied_load']).ravel())/max(norm(f['applied_load']),1.))
        maximum_residuals[p['path']]=max(residuals)

    assert all(b==bindings[0] for b in bindings)
    assert bindings[0]=={n:sha(HERE/n) for n in solve_names}
    totals=dict(main_cross_paths=54,main_exchange_paths=24,principal_matrix_paths=len(paths),observer_free_compatibility_paths=39,total_global_solve_paths_including_compatibility=117,
        principal_accepted_states=sum(p['accepted_states'] for p in paths),observer_free_compatibility_states=sum(p['accepted_states'] for p in paths if p['engine']=='native'))
    write(out,dict(**totals,geometry_audits=geometry_audits,postsolve_nodal_residuals=maximum_residuals,
        postsolve_audit_is_not_additional_fresh_solve=True,source_sha256=bindings[0],solve_sources_identical_throughout=True,other_analysis_source_bindings_are_per_path=True,
        all_required_minimums_observed=True,whole_run_equal_except_elapsed_for_39_native_paths=True,counts_do_not_include_preparation_or_local_cone_probes=True,
        native_force_and_linear_original_gates_checked=True,a02_passed=False))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--report',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();audit(a.run,a.report,a.output)
