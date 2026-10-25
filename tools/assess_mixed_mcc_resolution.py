"""Apply unchanged MCC field gates to explicitly listed resolution cases.

This numerical assessment is not a full feature acceptance. It records all
failed targets and the actual solver/reference inputs instead of dropping an
unfavorable case or replacing a target with a solver's convergence tolerance.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import numpy as np
import fem_reference as fem
from compare_mcc_rate_fem import compare
from mixed_mcc_io import validate
import compare_mcc_rate_fem
import mixed_mcc_io


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def physics(case):
    value=fem.normalize_stage_case(copy.deepcopy(case))
    for key in ('initial_increment','material_maximum_substep','linear_solver'):
        value['settings'].pop(key,None)
    return value


def reference_changes(coarse,fine,n):
    assert len(coarse['frames'])==len(fine['frames'])
    rows=[]
    for a,b in zip(coarse['frames'],fine['frames']):
        assert a['name']==b['name'] and a['active_elements']==b['active_elements']
        u=np.asarray(a['displacement_m'])-b['displacement_m']
        s=np.asarray(a['material_state']);t=np.asarray(b['material_state'])
        rows.append(dict(vertical_subdivisions=n,name=b['name'],
            maximum_displacement_change_m=float(abs(u).max()),
            maximum_pc_change_kpa=float(abs(np.exp(s[...,6])-np.exp(t[...,6])).max()),
            relative_stress_change=float(np.linalg.norm(s[...,:6]-t[...,:6])/max(np.linalg.norm(t[...,:6]),1e-30)),
            maximum_void_ratio_change=float(abs(np.exp(s[...,7])-np.exp(t[...,7])).max()),
            maximum_plastic_strain_change=float(abs(s[...,8:]-t[...,8:]).max())))
    return rows


def assess_references(study,output,resolutions):
    rows=[];gates=[];inputs={str(Path(__file__)):sha(Path(__file__))}
    for n in resolutions:
        paths=[study/f'n{n}/reference-{count}.json' for count in (400,800)]
        coarse,fine=map(read,paths)
        assert coarse['reference_steps']==400 and fine['reference_steps']==800
        for reference in (coarse,fine):
            for name,digest in reference['input_sha256'].items():
                assert sha(Path(name))==digest,'reference source or input changed: '+name
                inputs[name]=digest
            original=[Path(p) for p in list(reference['input_sha256'])[:3]]
            assert read(original[0])==read(study/f'n{n}/mesh.json')
            assert read(original[1])==read(study/'ground.json')
            assert physics(read(original[2]))==physics(read(study/f'n{n}/stages-ultrafine.json'))
            # The original CLI fingerprints its four direct numerical modules.
            # The preparation additionally freezes the imported point_state
            # helper before execution; verify this transitive module as well.
            source_root=next(Path(p).parent for p in reference['input_sha256'] if Path(p).name=='mcc_rate_fem.py')
            preparation=source_root.parent/'preparation.json'
            manifest=read(preparation)['created_inputs_sha256']
            for p in source_root.glob('*.py'):
                assert sha(p)==manifest[p.relative_to(source_root.parent).as_posix()]
                inputs[str(p)]=sha(p)
            inputs[str(preparation)]=sha(preparation)
            for frame in reference['frames']:
                residual=frame['relative_free_force_residual']
                gates.append(dict(name=f'n{n}:{reference["reference_steps"]}:{frame["name"]}:force',observed=residual,limit=1e-9,passed=bool(np.isfinite(residual) and residual<1e-9)))
        changes=reference_changes(coarse,fine,n);rows.extend(changes)
        for row in changes:
            for key,limit in (('maximum_displacement_change_m',1e-8),('maximum_pc_change_kpa',2.5e-5),
                              ('relative_stress_change',2.5e-6),('maximum_void_ratio_change',1.25e-6),('maximum_plastic_strain_change',2.5e-7)):
                value=row[key];gates.append(dict(name=f'n{n}:{row["name"]}:{key}',observed=value,limit=limit,passed=bool(np.isfinite(value) and value<limit)))
        inputs.update({str(p):sha(p) for p in paths})
    report=dict(schema='soil-fem-mcc-reference-integration-refinement/1',passed=all(g['passed'] for g in gates),
        goal_complete=False,p6_complete=False,native_accuracy_accepted=False,
        scope='Independent rate-form FEM 400->800 integration convergence only; no native Rust field acceptance',
        changes=rows,gates=gates,input_and_source_sha256=inputs)
    with output.open('x',encoding='utf-8') as stream:json.dump(report,stream,indent=2,allow_nan=False)
    print(json.dumps(dict(passed=report['passed'],changes=rows,failed_gates=[g for g in gates if not g['passed']])),flush=True)
    if not report['passed']:raise SystemExit(1)


def assessment(study,output,native_name,case_name,resolutions):
    if not resolutions or resolutions!=sorted(set(resolutions)) or any(n<2 or n%2 for n in resolutions):
        raise ValueError('resolutions must be increasing unique positive even counts')
    sources=[Path(__file__),Path(fem.__file__),Path(compare_mcc_rate_fem.__file__),Path(mixed_mcc_io.__file__)]
    gates=[];reports=[];refinements=[];inputs={str(p):sha(p) for p in sources}
    solver=study/'solver/soil-fem-3d.exe';inputs[str(solver)]=sha(solver)
    shared_physics=None
    def gate(name,value,limit,required=True):
        gates.append(dict(name=name,observed=float(value),limit=limit,
                          required_for_finest_resolution_acceptance=required,
                          passed=bool(np.isfinite(value) and value<limit)))
    for n in resolutions:
        where=study/f'n{n}';paths=[where/'mesh.json',study/'ground.json',where/case_name,where/native_name,
                                 where/'reference-400.json',where/'reference-800.json']
        mesh,ground,case,native,coarse,reference=map(read,paths)
        assert coarse['reference_steps']==400 and reference['reference_steps']==800,'incorrect reference integration levels'
        if not native['complete']:raise ValueError(f'n{n}: incomplete native result')
        assert native['provenance']['solver_binary_sha256']==sha(solver),'solver changed'
        assert fem.normalize_stage_case(copy.deepcopy(native['case']))==fem.normalize_stage_case(copy.deepcopy(case)), 'native case differs'
        assert fem.normalize_ground(copy.deepcopy(native['ground']))==fem.normalize_ground(copy.deepcopy(ground)), 'native ground differs'
        for result in (coarse,reference):
            for name,digest in result['input_sha256'].items():
                path=Path(name);assert sha(path)==digest, 'reference source or input changed: '+name
                inputs[name]=digest
            original_mesh,original_ground,original_case=[Path(p) for p in list(result['input_sha256'])[:3]]
            assert read(original_mesh)==mesh and read(original_ground)==ground,'reference geometry or materials differ'
            assert physics(read(original_case))==physics(case),'reference physical case differs'
        if shared_physics is None:shared_physics=physics(case)
        assert physics(case)==shared_physics,'spatial refinement changed physics'
        validate(mesh,native)
        comparison=compare(reference,native,mesh);comparison['vertical_subdivisions']=n;reports.append(comparison)
        changes=reference_changes(coarse,reference,n)
        for i,row in enumerate(comparison['frames']):
            prefix=f'n{n}:{row["name"]}:'
            for key,limit in (('displacement_relative_all_nodes',3e-5),('displacement_maximum_component_error_m',1e-6),
                              ('independent_relative_force_residual',1e-9)):
                gate(prefix+key,row[key],limit,required=n==resolutions[-1] or key=='independent_relative_force_residual')
            for field,key,limit in (('stress','relative_all_points',1e-5),('pc','maximum_absolute_all_points',1e-4),
                                    ('void_ratio','maximum_absolute_all_points',5e-6),('plastic_strain','maximum_absolute_all_points',1e-6)):
                gate(prefix+field,row['fields'][field][key],limit,required=n==resolutions[-1])
            frame=native['frames'][i];tolerance=case['settings']['relative_force_tolerance']*1.01
            gate(prefix+'native_force',frame['relative_force_residual'],tolerance)
            for name in ('relative_volume_equation_residual','relative_pressure_consistency_residual'):
                gate(prefix+name,frame['state']['mixed_mcc'][name],tolerance)
            refinement=changes[i]
            refinements.append(refinement)
            gate(prefix+'reference_400_to_800_displacement',refinement['maximum_displacement_change_m'],1e-8)
            gate(prefix+'reference_400_to_800_pc',refinement['maximum_pc_change_kpa'],2.5e-5)
            gate(prefix+'reference_400_to_800_stress',refinement['relative_stress_change'],2.5e-6)
            gate(prefix+'reference_400_to_800_void',refinement['maximum_void_ratio_change'],1.25e-6)
            gate(prefix+'reference_400_to_800_plastic_strain',refinement['maximum_plastic_strain_change'],2.5e-7)
        inputs.update({str(p):sha(p) for p in paths})
    if len(reports)>1:
        a,b=reports[-2:]
        for i in range(1,len(b['frames'])):
            error=lambda c:c['frames'][i]['fields']['pc']['maximum_absolute_all_points']
            gate(f'spatial_pc_error_ratio:{a["vertical_subdivisions"]}->{b["vertical_subdivisions"]}:frame{i}',error(b)/max(error(a),1e-30),1.)
    required=[g for g in gates if g['required_for_finest_resolution_acceptance']]
    report=dict(schema='soil-fem-mixed-mcc-resolution-assessment/1',passed=all(g['passed'] for g in required),
        all_coarse_and_fine_field_targets_passed=all(g['passed'] for g in gates),candidate_resolution=resolutions[-1],
        goal_complete=False,p6_complete=False,development_only=True,
        scope='Dry NC column loading/unloading/upper-half excavation; all node/GP accuracy and independent integration refinement',
        limits_origin='Existing MCC cycle field gates unchanged, applied to the finest candidate as in prior cycle refinement; all coarse failures retained. Reference 400->800 and equation checks required at every resolution.',
        comparisons=reports,reference_refinements=refinements,gates=gates,input_and_source_sha256=inputs,
        limitations='No acceptance for general terrain, structures, water, native UI or the full P6 scope')
    with output.open('x',encoding='utf-8') as stream:json.dump(report,stream,indent=2,allow_nan=False)
    print(json.dumps(dict(passed=report['passed'],failed_required_gates=[g for g in required if not g['passed']],
                         coarse_diagnostic_failures=[g for g in gates if not g['passed'] and not g['required_for_finest_resolution_acceptance']],
                         spatial_maximum_pc_error_kpa=[dict(n=c['vertical_subdivisions'],error=max(r['fields']['pc']['maximum_absolute_all_points'] for r in c['frames'])) for c in reports])),flush=True)
    if not report['passed']:raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('study',type=Path);parser.add_argument('output',type=Path)
    parser.add_argument('--native-name',default='native-ultrafine.json');parser.add_argument('--case-name',default='stages-ultrafine.json')
    parser.add_argument('--resolutions',type=int,nargs='+',default=[16,32])
    parser.add_argument('--reference-only',action='store_true')
    a=parser.parse_args()
    if a.reference_only:assess_references(a.study.resolve(),a.output.resolve(),a.resolutions)
    else:assessment(a.study.resolve(),a.output.resolve(),a.native_name,a.case_name,a.resolutions)
