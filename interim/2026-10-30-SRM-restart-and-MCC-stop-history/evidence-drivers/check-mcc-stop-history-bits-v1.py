"""Check the already computed real MCC histories without repeating a solve."""
from pathlib import Path
import hashlib,json,struct
p=Path(__file__).resolve().parent;source=p/'mcc-stop-boundary-verification-v1'
out=p/'mcc-stop-history-bits-v1';out.mkdir()
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
binding={};rows=[]
for formulation in ['displacement','mixed']:
    old=source/'results-653'/formulation;new=source/'results-656'/formulation
    for name in ['mesh.json','ground.json','case.json']:
        assert sha(old/name)==sha(new/name)
        binding[str(old/name)]=sha(old/name);binding[str(new/name)]=sha(new/name)
    full=read(new/'uninterrupted.json')
    compare(read(old/'uninterrupted.json'),full,exclude_executable=True)
    assert read(old/'uninterrupted.json')['provenance']['solver_binary_sha256']!=full['provenance']['solver_binary_sha256']
    for boundary in ['completed-geostatic','plastic-construction']:
        before=new/(boundary+'-before.json');stopped=new/(boundary+'-stopped.json');resumed=new/(boundary+'-resumed.json')
        assert before.read_bytes()==stopped.read_bytes(),(formulation,boundary,'stopped bytes')
        compare(full,read(resumed))
        for q in [before,stopped,resumed,new/'uninterrupted.json',old/'uninterrupted.json']:binding[str(q)]=sha(q)
        rows.append(dict(formulation=formulation,boundary=boundary,stopped_entire_JSON_bytes_identical=True,
            resumed_all_fields_types_and_binary64_bits_match_uninterrupted_except_elapsed_seconds=True))
negative=read(source/'results-656/mixed/uninterrupted.json')
negative['last_committed']['material_states'][0]['void_ratio']+=1e-12
try:compare(read(source/'results-656/mixed/uninterrupted.json'),negative)
except AssertionError as error:negative_diagnostic=str(error)
else:raise AssertionError('changed hardening history was not detected')
binding[str(source/'source-binding.json')]=sha(source/'source-binding.json')
binding[str(Path(__file__))]=sha(Path(__file__))
for q,v in [('source-binding.json',binding),('terminal.json',dict(passed=True,rows=rows,
    old_and_new_uninterrupted_fields_also_bitwise_equal_except_elapsed_and_exact_executable_provenance=True,
    changed_void_ratio_1e_minus12_rejected=negative_diagnostic,no_FEM_repeated=True,p6_complete=False))]:
    with (out/q).open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
print('PASS actual MCC: stopped bytes, complete resumed history bits, old/new unchanged physics',flush=True)
