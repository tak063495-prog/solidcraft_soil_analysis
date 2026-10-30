"""Inspect complete actual cut/restart histories; never rerun the FEM."""
from pathlib import Path
import hashlib,json,struct
p=Path(__file__).resolve().parent;source=p/'mcc-local-cut-stop-verification-v1';out=p/'mcc-local-cut-history-bits-v1';out.mkdir()
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
    with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def compare(a,b,path='$',exclude_executable=False):
    assert type(a) is type(b),(path,type(a),type(b))
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for key in a:
            if key=='elapsed_seconds':continue
            if exclude_executable and key=='solver_binary_sha256':continue
            compare(a[key],b[key],path+'.'+key,exclude_executable)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(av,bv) in enumerate(zip(a,b)):compare(av,bv,f'{path}[{i}]',exclude_executable)
    elif isinstance(a,float):assert struct.pack('>d',a)==struct.pack('>d',b),(path,a,b)
    else:assert a==b,(path,a,b)
assert read(source/'terminal.json')['passed']
rows=[];binding={}
for family in ['mixed','projected']:
    new=source/'results-656'/family;old=source/'results-653'/family
    for name in ['mesh.json','ground.json','case.json']:assert sha(new/name)==sha(old/name)
    full=read(new/'uninterrupted.json');before=read(new/'before.json');late=read(new/'late-stop.json');mesh=read(new/'mesh.json')
    compare(read(old/'uninterrupted.json'),full,exclude_executable=True)
    assert read(old/'uninterrupted.json')['provenance']['solver_binary_sha256']!=full['provenance']['solver_binary_sha256']
    assert (new/'before.json').read_bytes()==(new/'stopped.json').read_bytes()
    for name in ['resumed.json','late-stop-resumed.json']:compare(full,read(new/name))
    cut=mesh['excavation_sets']['LocalCut'];assert len(cut)==1
    dormant=[i for i,q in enumerate(before['integration_points']) if q['element_index'] in cut]
    assert len(dormant)==27
    for value in [before,late,full]:
        s=value['last_committed'];assert s['active_elements'].count(False)==1
        assert all(not s['active_elements'][ei] for ei in cut)
        for i in dormant:compare(before['last_committed']['material_states'][i],s['material_states'][i])
    for run in [before,late]:
        state=run['last_committed'];holding=run['excavation_audits'][0]['initial_holding_load_kn']
        assert state['surface_pressure_kpa']==0. and len(holding)==len(state['excavation_release_load_kn'])==3*len(mesh['nodes'])
        assert max(abs(x) for x in holding)>1e-8
        for actual,initial in zip(state['excavation_release_load_kn'],holding):
            assert struct.pack('>d',actual)==struct.pack('>d',(1-state['stage_progress'])*initial)
    vertices=[mesh['nodes'][i] for i in mesh['tetrahedra'][cut[0]]]
    columns=[[vertices[j+1][i]-vertices[0][i] for i in range(3)] for j in range(3)]
    a,b,c=columns
    determinant=a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
    volume=abs(determinant)/6
    assert abs(before['excavation_audits'][0]['removed_volume_m3']-volume)<1e-14
    for q in [*new.glob('*.json'),*old.glob('*.json')]:binding[str(q)]=sha(q)
    rows.append(dict(formulation=family,elements=12,integration_points=324,dormant_points_bitwise_frozen=27,
        removed_volume_m3=volume,stopped_Run_entire_bytes_identical=True,
        early_and_late_resumed_entire_histories_types_and_binary64_bits_equal_except_elapsed_seconds=True,
        release_force_every_component_exact_to_remaining_accepted_progress=True,
        actual_full_and_partial_local_cut_active_masks_and_zero_floor_surcharge_checked=True,
        old_and_new_uninterrupted_physics_bitwise_equal_except_elapsed_and_executable=True))
bad=read(source/'results-656/mixed/late-stop.json');bad['last_committed']['excavation_release_load_kn'][0]+=1e-12
try:compare(read(source/'results-656/mixed/late-stop.json'),bad)
except AssertionError as error:negative=str(error)
else:raise AssertionError('changed release history was not detected')
binding[str(source/'source-binding.json')]=sha(source/'source-binding.json');binding[str(Path(__file__))]=sha(Path(__file__))
for name,value in [('source-binding.json',binding),('terminal.json',dict(passed=True,rows=rows,
    release_component1e_minus12_corruption_rejected=negative,no_FEM_repeated=True,
    history_and_bookkeeping_acceptance_not_independent_nonlinear_solution_or_spatial_convergence=True,p6_complete=False))]:
    with (out/name).open('x',encoding='utf-8') as f:json.dump(value,f,indent=2,allow_nan=False)
print('PASS actual local MCC: cut history, early/late stop, dormant GP and release force bits',flush=True)
