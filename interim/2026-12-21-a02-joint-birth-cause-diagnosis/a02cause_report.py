"""Rebuild derived diagnosis from sealed observations, without a new solve."""
import argparse
import json
from pathlib import Path
import math
import numpy as np
from a02cause_binding import fields,study
from a02cause_driver import CASES,write
from a02cause_metrics import norm


def report(run,active,partitions,out):
    out.mkdir(parents=True,exist_ok=True)
    m=study.read(run/'measurements.json');p=study.read(run/'projection-observations/diagnosis.json')
    a=study.read(active/'active-contact-controls.json');c=study.read(partitions/'partition-global-controls.json')
    write(out/'measurements.json',m);write(out/'projection-diagnosis.json',p)
    write(out/'active-contact-controls.json',a);write(out/'partition-global-controls.json',c)
    result=dict(input_only_reference=True,continuum_reference_established=False,spatial_convergence_established=False,
        all_field_acceptance_defined=False,cloud_diagnostics_remaining=True,extrapolation_used=False,
        quadrature_changes_not_mesh_refinement=True,cases={})
    for name,_,_ in CASES:
        records=m['cases'][name];rows={}
        for source,folder in [('baseline',run/'n4'/name),('minimum2',partitions/name)]:
            error=0.
            for row,f in fields.load_fields(folder/('independent-fields.jsonl')):
                v=f['initial_energy']+f['endpoint_work']-f['elastic_energy']-f['plastic_dissipation']-f['numerical_dissipation']
                error=max(error,float(np.max(abs(v))))
            rows[source]=dict(maximum_contact_work_identity_error_kpa_m=error)
        result['cases'][name]=dict(native_grid=records['native_grid_sensitivity'],
            reference_grid=records['reference_grid_sensitivity'],independent_controls=records['independent_controls'],
            active_increment_control=a[name],minimum_partition_control=c[name],work_identity=rows,
            local_replay_not_global_equilibrium=True)
    write(out/'reference-assessment.json',result)
    accounting={}
    for name,_,_ in CASES:
        rows=list(fields.load_fields(run/'n4'/name/'native-fields.jsonl'))
        detail=[]
        previous=None
        for row,f in rows:
            # Preserve the original generation-based identity. The independent
            # global control above determines whether this accounting explains
            # a discretization effect rather than an equilibrium failure.
            err=f['initial_energy']+f['endpoint_work']-f['elastic_energy']-f['plastic_dissipation']-f['numerical_dissipation']
            detail.append(dict(stage=row['stage'],progress=row['progress'],generation=f['generation'].tolist(),
                work_identity_error_kpa_m=err.tolist(),maximum_absolute=float(np.max(abs(err)))))
        accounting[name]=dict(saved_endpoint_work_definition_unchanged=True,rows=detail,
            local_replay=m['cases'][name]['local_replay'])
    write(out/'work-accounting.json',accounting)
    write(out/'cause-assessment.json',dict(passed=False,a02_complete=False,single_cause_fully_identified=False,
        S1=dict(main='activation-domain/initial-state/equilibrium specification missing',numerical_contributions_identifiable=False),
        S2=dict(main='partial-owner GP/generation specification missing',numerical_contributions_identifiable=False),
        S3=dict(main='initial projected traction violates strict no-tension; not projection nonreproducibility',local_repair_defined=False),
        S4=dict(main='initial projected traction violates strict Coulomb; not a subsequent-increment trial state',local_repair_defined=False),
        S5=dict(main='nodal residual outside contact range; shared-tip rows identically zero',maximum_range_component_kn=p['unreproducible_range_component_maximum_kn']),
        S6=dict(main='tiny absolute roundoff plus existing denominator floor; broad field criterion undefined',new_threshold_defined=False),
        S7=dict(main='endpoint work depends on backward-Euler subdivision; cycle forces/slip depend on global equilibrium path',
                local_return_defect_proved=False,independent_reference_step_sensitive=True,mesh_convergence_established=False),
        S8=dict(main='retained infrastructure/measurement/authentication failures stay distinct',past_failures_removed=False),
        cloud_diagnostics_remaining=True,pc_long_series_optional=True,user_specification_decisions_needed=True))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for flag in ('run','active','partitions','output'):p.add_argument('--'+flag,type=Path,required=True)
    a=p.parse_args();report(a.run,a.active,a.partitions,a.output)


if __name__=='__main__':main()
