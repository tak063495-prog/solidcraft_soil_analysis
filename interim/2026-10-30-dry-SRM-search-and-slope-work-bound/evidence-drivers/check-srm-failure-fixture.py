"""Reproduce the failed search-test origin using unchanged adopted mechanics."""
from pathlib import Path
import copy, datetime, hashlib, json, os, subprocess, time
p=Path(__file__).resolve().parent;root=p.parents[1]
out=p/'srm-adaptive-search-v1/failure-fixture-diagnostic-v1';out.mkdir()
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def dump(q,v):
    with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
assert read(p/'srm-adaptive-search-v1/quality-v1/terminal.json')['failed_step']=='new-ten-search-tests'
assert not (p/'pause-requested').exists()
dump(out/'driver-child.json',dict(pid=os.getpid(),script=str(Path(__file__).resolve()),
    executable=os.sys._base_executable,created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
mesh=read(root/'examples/mc-column-mesh.json');linear_nodes=copy.deepcopy(mesh['nodes']);edges={};high=[]
for tet in mesh['tetrahedra']:
    cell=tet.copy()
    for i,j in [(0,1),(1,2),(2,0),(0,3),(1,3),(2,3)]:
        edge=tuple(sorted((tet[i],tet[j])))
        if edge not in edges:
            edges[edge]=len(mesh['nodes'])
            mesh['nodes'].append([(linear_nodes[edge[0]][a]+linear_nodes[edge[1]][a])/2. for a in range(3)])
        cell.append(edges[edge])
    high.append(cell)
for name, original in list(mesh['node_sets'].items()):
    selected=set(original)
    mesh['node_sets'][name]=sorted(selected|{n for (a,b),n in edges.items() if a in selected and b in selected})
mesh['schema']='soil-fem-mesh/3';mesh['tetrahedra10']=high
dump(out/'mesh.json',mesh)
ground=read(root/'examples/ground-geostatic-mc-column.json')
for material in ground['materials']:
    material['unit_weight_kn_m3']=1e-10
    material['mohr_coulomb'].update(cohesion_kpa=40.,friction_deg=0.,dilation_deg=0.,tension_cutoff_kpa=None)
case=read(root/'examples/stages-geostatic-column.json')
case['stages']=[dict(kind='pressure',name='load',target_pressure_kpa=100.)]
case['verification_reject_once_stages']=[]
case['settings'].update(maximum_newton_iterations=2,initial_increment=1.,minimum_increment=1.)
solver=root/'target/release/soil-fem-3d.exe'
assert hashlib.sha256(solver.read_bytes()).hexdigest()=='a4b5c173be10dca66ef90a06db3a4f1dd7246c3f614bcca41cc0eb58dc63e114'
rows=[]
def trial(label,ground,case):
    g=out/(label+'-ground.json');c=out/(label+'-case.json');result=out/(label+'-result.json')
    dump(g,ground);dump(c,case)
    start=time.monotonic()
    args=[str(solver),'solve-stages',str(out/'mesh.json'),str(g),str(c),str(result)]
    with (out/(label+'.log')).open('x',encoding='utf-8') as log:
        child=subprocess.Popen(args,cwd=root,stdout=log,stderr=subprocess.STDOUT,
            creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        dump(out/(label+'-child.json'),dict(pid=child.pid,parent_pid=os.getpid(),executable=str(solver),
            arguments=args[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        code=child.wait()
    r=read(result)
    row=dict(label=label,exit_code=code,complete=r['complete'],failure_reason=r['failure_reason'],
        geostatic_complete=r['geostatic']['complete'],geostatic_reason=r['geostatic']['failure_reason'],
        geostatic_attempts=[dict(accepted=a['accepted'],reason=a['reason'],restored=a['committed_state_unchanged_on_rejection']) for a in r['geostatic']['attempts']],
        construction_attempts=[dict(accepted=a['accepted'],reason=a['reason'],restored=a['committed_state_unchanged_on_rejection']) for a in r['attempts']],
        elapsed_seconds=time.monotonic()-start)
    rows.append(row);dump(out/(label+'-terminal.json'),row)
    return row
first=trial('original-failing-origin',ground,case)
assert not first['complete'] and not first['geostatic_complete']
assert any('verification' in str(a['reason']).lower() for a in first['geostatic_attempts'])
ground['initialization']['geostatic']['verification_reject_once']=False
chosen=None
for cap in [2,3,4,5]:
    case['settings']['maximum_newton_iterations']=cap
    normal=trial(f'cap{cap}-F1',ground,case)
    reduced=copy.deepcopy(ground)
    for material in reduced['materials']:material['mohr_coulomb']['cohesion_kpa']/=4.
    weak=trial(f'cap{cap}-F4',reduced,case)
    if normal['complete'] and not weak['complete']:
        chosen=cap;break
assert chosen is not None,'no verified strength-dependent nonlinear failure control found'
dump(out/'terminal.json',dict(passed=True,actual_original_failure_reproduced=True,
    original_failure_was_geostatic_verification_rejection_with_no_cutback_reserve=True,
    no_original_physical_upper_was_created=True,corrected_only_test_injection_and_iteration_budget=True,
    all_force_yield_strain_and_linear_tolerances_unchanged=True,
    selected_maximum_newton_iterations=chosen,actual_F1_complete_and_F4_incomplete=True,
    rows=rows,unchanged_adopted_A4_mechanics_used=True,not_an_independent_theory_validation=True,p6_complete=False))
print(json.dumps(dict(passed=True,chosen_newton_cap=chosen,actual_F1_complete_F4_incomplete=True)),flush=True)
