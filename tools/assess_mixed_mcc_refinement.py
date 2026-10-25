"""Apply the existing MCC field-accuracy targets to mixed development results.

This assessment deliberately records failed targets. It is neither a solver
nor a scoped acceptance producer; lifecycle, export and source gates are separate.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import numpy as np
from compare_mcc_rate_fem import compare


def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def assess(base,output):
    gates=[];comparisons=[];refinements=[]
    def gate(name,value,limit):gates.append(dict(name=name,observed=float(value),limit=limit,passed=bool(value<limit)))
    for n in (2,4):
        where=base/f'n{n}';mesh=read(where/'mesh.json')
        reference=read(where/('reference-400.json' if n==4 else 'reference-100.json'))
        for name,digest in reference['input_sha256'].items():
            path=Path(name)
            if path.suffix=='.py':path=base/'reference-v5'/path.name
            assert sha(path)==digest,name
        cases=[]
        for tag in ('','-refined','-finest'):
            native=read(where/('native'+tag+'.json'));assert native['complete']
            assert native['provenance']['solver_binary_sha256']==sha(base/'solver/soil-fem-3d.exe')
            case=read(where/('stages'+tag+'.json'));physics=copy.deepcopy(case)
            for key in ('initial_increment','material_maximum_substep'):physics['settings'].pop(key)
            cases.append(physics)
            result=compare(reference,native,mesh);result['resolution']=n;result['numerical_level']=tag or 'coarse'
            comparisons.append(result)
            if tag=='-finest':
                for row in result['frames']:
                    prefix=f'n{n}:{row["name"]}:'
                    for key,limit in (('displacement_relative_all_nodes',3e-5),('displacement_maximum_component_error_m',1e-6),('independent_relative_force_residual',1e-9)):
                        gate(prefix+key,row[key],limit)
                    for field,key,limit in (('stress','relative_all_points',1e-5),('pc','maximum_absolute_all_points',1e-4),
                                             ('void_ratio','maximum_absolute_all_points',5e-6),('plastic_strain','maximum_absolute_all_points',1e-6)):
                        gate(prefix+field,row['fields'][field][key],limit)
                for i,frame in enumerate(native['frames']):
                    mixed=frame['state']['mixed_mcc'];tolerance=case['settings']['relative_force_tolerance']*1.01
                    gate(f'n{n}:frame{i}:volume',mixed['relative_volume_equation_residual'],tolerance)
                    gate(f'n{n}:frame{i}:pressure',mixed['relative_pressure_consistency_residual'],tolerance)
        assert cases[0]==cases[1]==cases[2],'numerical refinement changed physics'
    levels=[read(base/f'n4/reference-{k}.json') for k in (50,100,200,400)]
    for i in range(4):
        u=[np.asarray(d['frames'][i]['displacement_m']) for d in levels]
        changes=[float(np.max(abs(a-b))) for a,b in zip(u,u[1:])]
        refinements.append(dict(frame=i,reference_maximum_displacement_changes_m=changes))
        if i:
            gate(f'reference:frame{i}:latest_displacement_change_m',changes[-1],1e-8)
            gate(f'reference:frame{i}:latest_change_ratio',changes[-1]/max(changes[-2],1e-30),1.)
    report=dict(schema='soil-fem-mixed-mcc-field-refinement-assessment/1',passed=all(g['passed'] for g in gates),
        goal_complete=False,p6_complete=False,development_only=True,
        scope='Common dry NC two-layer column: force loading/unloading/upper-layer removal; numerical field targets only',
        limits_origin='Same field targets used in assess_mcc_cycle_refinement.py; not relaxed for the new formulation',
        comparisons=comparisons,reference_refinements=refinements,gates=gates,
        verifier_sha256=sha(Path(__file__)),
        input_and_evidence_sha256={p.relative_to(base).as_posix():sha(p) for p in base.rglob('*') if p.is_file()})
    with output.open('x',encoding='utf-8') as f:json.dump(report,f,indent=2)
    print(json.dumps(dict(passed=report['passed'],failed_gates=[g for g in gates if not g['passed']],reference_refinements=refinements)),flush=True)
    if not report['passed']:raise SystemExit(1)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('base',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();assess(a.base.resolve(),a.output.resolve())
