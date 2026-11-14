"""Supplemental audit of NEWLY solved A-04 fields; never counted as fresh solve."""
import argparse
import itertools
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE))
from study import read,write,fields,exact
from soil_reference import constraints


def audit(folder):
    native_input_mesh=read(folder/'mesh.json')
    native=read(folder/'native.json');comparison=read(folder/'comparison.json');independent=list(fields(folder/'independent-fields.jsonl'));observed=list(fields(folder/'native-fields.jsonl'))
    frames=native['initialization']['frames'][1:]+native['frames'];assert len(frames)==len(independent)==len(observed)
    maxima={};records=[];previous=None;unsaturated=0;native_mass=independent_mass=0.
    def error(key,a,b):
        value=float(np.max(np.abs(np.asarray(a)-np.asarray(b)),initial=0));maxima[key]=max(maxima.get(key,0.),value);return value
    for frame,(row,actual),(refrow,reference),point in itertools.zip_longest(frames,observed,independent,comparison['all_points']):
        assert row==refrow==point['row']
        for k in ('u','stress'):assert point['fields'][k]['relative_l2']<1e-7,(row,k,point['fields'][k])
        for k in ('effective_internal','fluid_internal','water_body','exterior','physical_load','equivalent_rhs','applied_load','reaction','holding_load'):
            assert point['fields'][k]['maximum_absolute']<1e-10,(row,k,point['fields'][k])
        assert all(v.get('mismatches',0)==0 for v in point['fields'].values())
        _,fixed=constraints(native_input_mesh,actual['active_elements'],'side_rollers')
        actual_support=frame['mechanical_audit']['support_reaction_kn'];ref_support=np.where(fixed,reference['reaction'],0).sum(axis=0)
        support_error=error('production_support_vs_independent_kn',actual_support,ref_support) # resultant is recorded; no new tolerance is inferred from a nodal gate
        actual_reassembled_support=np.where(fixed,actual['reaction'],0).sum(axis=0)
        kernel_error=error('production_support_vs_saved_stress_reassembly_kn',actual_support,actual_reassembled_support);assert kernel_error<1e-10
        physical=reference['physical_load']+(reference['holding_load'] if row['phase']=='excavation_release' else 0)
        balance=physical.sum(axis=0)+ref_support
        balance_error=error('production_force_balance_vs_independent_kn',frame['mechanical_audit']['physical_force_balance_error_kn'],balance) # separate physical resultant diagnostic
        for tag,f in (('native',actual),('independent',reference)):
            assert abs(float(f['mass_defect_m3']))<1e-10
            error(tag+'_water_load_decomposition_kn',f['physical_load']-f['fluid_internal'],f['equivalent_rhs'])
            error(tag+'_holding_decomposition_kn',f['equivalent_rhs']+f['holding_load'],f['applied_load'])
            # Source-owned geometry/active masks, no filling undefined cells.
            assert np.array_equal(f['cell_ids'],np.flatnonzero(f['active_elements']))
            if tag=='native':native_mass=max(native_mass,abs(float(f['mass_defect_m3'])))
            else:independent_mass=max(independent_mass,abs(float(f['mass_defect_m3'])))
            if previous is not None:
                prev=previous[tag];dead=~f['active_elements']
                for k in ('stress','strain','void','plastic_strain','equivalent_plastic_strain','plastic_work'):
                    assert np.array_equal(f[k][dead],prev[k][dead]),(tag,row,k,'inactive material history changed')
                used=set(int(n) for e in f['cell_ids'] for n in native_input_mesh['tetrahedra'][e]);unused=sorted(set(range(len(f['u'])))-used)
                assert np.array_equal(f['u'][unused],prev['u'][unused]),(tag,row,'inactive nodal state changed')
                if row['phase']=='excavation_release' and previous['row']['phase']=='excavation_release':
                    assert row['time_s']==previous['row']['time_s']
                    for k in ('cell_head_m','face_head_m','outflow_m3_s','storage_offset_m3','water_m3','removed_water_m3','net_inflow_m3'):assert np.array_equal(f[k],prev[k]),(tag,row,k,'release advanced water')
        unsaturated+=int(np.count_nonzero(actual['hydraulic_effective_saturation']<1))
        records.append(dict(row=row,actual_support_reaction_kn=actual_support,independent_support_reaction_kn=ref_support.tolist(),support_error_kn=support_error,
            actual_force_balance_error_kn=frame['mechanical_audit']['physical_force_balance_error_kn'],independent_force_balance_error_kn=balance.tolist(),
            native_mass_defect_m3=float(actual['mass_defect_m3']),independent_mass_defect_m3=float(reference['mass_defect_m3'])))
        previous=dict(native=actual,independent=reference,row=row)
    return dict(saved_new_solve_audit_only=True,new_solve_count=0,accepted_states=len(records),existing_pointwise_gates_passed=True,maximum_errors=maxima,
        native_maximum_mass_defect_m3=native_mass,independent_maximum_mass_defect_m3=independent_mass,unsaturated_hydraulic_point_states=unsaturated,
        all_accepted_points=records,physical_seconds_not_release_progress=True,removed_histories_frozen=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args();out=a.run.resolve();results={}
    for folder in sorted(out.glob('*-n*')):
        native_input_mesh=read(folder/'mesh.json');results[folder.name]=audit(folder)
    write(out/'supplemental-audit.json',results)
