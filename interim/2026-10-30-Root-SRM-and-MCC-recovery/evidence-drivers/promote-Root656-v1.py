"""Adopt the fully qualified source after preserving635 and reviewing readers."""
from pathlib import Path,PurePosixPath
import datetime,hashlib,json,os,shutil,struct,subprocess,sys,time,zipfile
p=Path(__file__).resolve().parent;root=p.parents[1];out=p/'Root-source656-promotion-v1';out.mkdir()
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
    with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
    q.parent.mkdir(parents=True,exist_ok=True)
    with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
def active():
    if (p/'pause-requested').exists():raise SystemExit(2)
def compare(a,b,old_binary,new_binary,path='$'):
    assert type(a)is type(b),(path,type(a),type(b))
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a:
            if k=='elapsed_seconds':continue
            if k=='solver_binary_sha256':assert a[k]==old_binary and b[k]==new_binary,(path,a[k],b[k]);continue
            compare(a[k],b[k],old_binary,new_binary,path+'.'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,old_binary,new_binary,f'{path}[{i}]')
    elif isinstance(a,float):assert struct.pack('>d',a)==struct.pack('>d',b),(path,a,b)
    else:assert a==b,(path,a,b)
active();dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
candidate=p/'srm-durable-job-v1/runtime';quality=p/'srm-durable-job-v1/quality-v3'
assert read(quality/'terminal.json')['passed'] and read(quality/'terminal.json')['Rust_passed']==347
base=p/'Root-source635-promotion-v1/quality-v1';old=read(base/'source-snapshot.json');new=read(quality/'source-snapshot.json')
sys.path.insert(0,str(root/'tools'));from freeze_sources import manifest
assert manifest()==old and len(old)==635 and len(new)==656
assert set(old)<=set(new) and all(sha(candidate/n)==h for n,h in new.items())
assert {n for n in old if old[n]!=new[n]}=={'src/lib.rs','src/main.rs','src/stages.rs'}
reader=read(p/'Root656-reader-audit-v1/inventory.json');assert not reader['direct_production_tool_processes']
# Final current-state check: every listed source consumer lives in a frozen run.
ps="Get-CimInstance Win32_Process | Select-Object ProcessId,ParentProcessId,ExecutablePath,CommandLine,@{n='Created';e={$_.CreationDate.ToUniversalTime().ToString('o')}} | ConvertTo-Json -Depth 2 -Compress"
current=json.loads(subprocess.check_output(['pwsh','-NoProfile','-Command',ps],encoding='utf-8-sig',creationflags=subprocess.CREATE_NO_WINDOW))
production_tools=str(root/'tools').lower().replace('\\','/');production_exe=root/'target/release/soil-fem-3d.exe'
assert not any(str(r.get('ExecutablePath') or '').lower().replace('\\','/')==str(production_exe).lower().replace('\\','/') for r in current)
assert not any(production_tools+'/' in (r.get('CommandLine') or '').lower().replace('\\','/') for r in current)
dump(out/'reader-review.json',dict(passed=True,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    previous_inventory=str(p/'Root656-reader-audit-v1/inventory.json'),previous_inventory_sha256=sha(p/'Root656-reader-audit-v1/inventory.json'),
    reviewed34350_future_runtime_extracted_from_frozen635_ZIP=True,reviewed_slope_auditor_uses_frozen_offline635=True,
    reviewed_global_comparison_uses_frozen_python_blocked_and_left_runtime=True,
    reviewed_owned4_native_host_uses_its_frozen_runtime=True,no_live_direct_Root_tools_or_Root_solver=True,
    no_existing_live_source_or_frozen_binary_modified=True,p6_complete=False))
before=out/'before';before.mkdir()
for name in ['source-snapshot.json','source-snapshot.zip']:shutil.copyfile(base/name,before/name);assert sha(base/name)==sha(before/name)
with zipfile.ZipFile(before/'source-snapshot.zip') as z:
    assert z.testzip() is None and set(z.namelist())==set(old) and len(z.namelist())==len(old)
    assert all(hashlib.sha256(z.read(n)).hexdigest()==h==sha(root/n) for n,h in old.items())
assert sha(production_exe)=='a4b5c173be10dca66ef90a06db3a4f1dd7246c3f614bcca41cc0eb58dc63e114'
old_exe=before/'Root635-soil-fem-3d.exe';shutil.copyfile(production_exe,old_exe);assert sha(old_exe)==sha(production_exe)
records=[]
for family in ['mixed','projected']:
    active();source=p/'mcc-local-cut-stop-verification-v1/results-653'/family;folder=before/family;folder.mkdir()
    args=[old_exe,'solve-stages',source/'mesh.json',source/'ground.json',source/'case.json',folder/'result.json']
    with (folder/'solve.log').open('x',encoding='utf-8') as log:
        child=subprocess.Popen(list(map(str,args)),cwd=root,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        dump(folder/'child.json',dict(pid=child.pid,parent_pid=os.getpid(),executable=str(old_exe),arguments=list(map(str,args[1:])),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        code=child.wait()
    assert code==0 and read(folder/'result.json')['complete']
    original=read(source/'uninterrupted.json');actual=read(folder/'result.json')
    compare(original,actual,original['provenance']['solver_binary_sha256'],sha(old_exe))
    records.append(dict(formulation=family,passed=True,old_Root635_CLI_all_history_bits_match_old653_probe_except_elapsed_and_verified_producer=True,report_sha256=sha(folder/'result.json')))
dump(out/'old-Root635-baselines.json',dict(passed=True,records=records,old_binary_sha256=sha(old_exe),p6_complete=False))
active();assert manifest()==old and all(sha(candidate/n)==h for n,h in new.items())
changed=sorted(n for n in new if old.get(n)!=new[n]);dump(out/'source-adoption-plan.json',dict(old_source_files=635,new_source_files=656,added=sorted(set(new)-set(old)),changed_existing=sorted(n for n in old if old[n]!=new[n]),all_existing_user_files_checked=True,old_source_and_binary_verified_backup=True))
for name in changed:
    path=PurePosixPath(name);assert not path.is_absolute() and '..' not in path.parts
    target=root/name;assert target.resolve().is_relative_to(root) and not target.is_symlink()
    if name in old:assert sha(target)==old[name]
    else:assert not target.exists()
    raw=(candidate/name).read_bytes();assert hashlib.sha256(raw).hexdigest()==new[name];target.parent.mkdir(parents=True,exist_ok=True)
    temporary=target.with_name(target.name+'.Root656-promotion.tmp')
    with temporary.open('xb') as f:f.write(raw);f.flush();os.fsync(f.fileno())
    os.replace(temporary,target)
assert manifest()==new
# Documentation reflects the adopted limited feature, preserving all35 items.
text='**2026-10-30追補：SRM探索・仮想仕事による離散釣合い排除・受理増分の保存/復旧を本体ソース656版へ取り込みました。同一の数値ソースはRust347件・Python219件に合格しました。本体配置のビルド・MCC履歴・SRM復旧・固定斜面の全GP検証を追加確認します。斜面の独立上端はF=1.45、元の2D変位誤差5％基準・一般収束・全35項目のP6完成は未達です。** [SRM手順](docs/SRM_SEARCH_JOBS.ja.md)\n\n'
for name in ['README.md','docs/P6_STATUS.ja.md']:
    prior=root/'runs/interim-upload-20261030-v53/publication'/name
    head,rest=prior.read_text('utf-8-sig').split('\n',1);intro=text if name=='README.md' else text.replace('(docs/SRM_SEARCH_JOBS.ja.md)','(SRM_SEARCH_JOBS.ja.md)')
    (root/name).write_text(head+'\n\n'+intro+rest.lstrip('\r\n'),'utf-8',newline='\n')
replacements={
 'docs/SRM_CONTROLLED_SEARCH.ja.md':('別の開発候補です。実行した検証ログが合格するまで採用版とは扱いません。現時点ではRustライブラリAPIであり、所有ジョブ・原子的ファイル保存・復旧用ディレクトリー・SolidCraftの操作画面は次工程です。',
    '本体656版へ接続したRust APIです。同じ数値ソースの全Rust347件・Python219件と、MCCの塑性/局所掘削中の停止・再開が合格しました。所有ジョブ・原子的ファイル保存・復旧用ディレクトリーは[保存型CLI](SRM_SEARCH_JOBS.ja.md)から使用できます。SolidCraftの途中工程ゲート操作画面は未接続です。'),
 'docs/SRM_SEARCH_JOBS.ja.md':('別の開発候補です。検証ログが合格するまで採用版とは扱いません。','本体656版へ接続した保存型CLIです。同一数値ソースの全Rust347件・Python219件と、実CLIの協調停止・改変拒否・保存障害・新規復旧が合格しました。'),
 'docs/SRM_SEARCH.ja.md':('開発中の別候補です。一般斜面・水圧・施工安全認定・P6全体の完了を意味しません。実行結果は版を固定した検証ログで確認します。',
    '本体656版へ接続した乾燥MCの限定機能です。同一数値ソースのRust347件・Python219件が合格しました。一般斜面・水圧・施工安全認定・P6全体の完了を意味しません。実行結果は版を固定した検証ログで確認します。')}
for name,(a,b) in replacements.items():
    q=root/name;v=q.read_text('utf-8-sig');assert v.count(a)==1;q.write_text(v.replace(a,b),'utf-8',newline='\n')
q=root/'docs/SRM_SEARCH.ja.md';v=q.read_text('utf-8');a='保存型の協調中断/再開、施工後の材料epochから始める別方式、固定水圧SRM、一般支持条件の独立上端、塑性斜面の2D精度・一般3D収束・SolidCraft入力UI/途中ゲートは別の受入項目です。'
assert v.count(a)==1;v=v.replace(a,'保存型の協調中断/再開は[保存型CLI](SRM_SEARCH_JOBS.ja.md)へ接続しました。施工後の材料epochから始める探索との統合、固定水圧SRM、斜め拘束を含む一般支持条件、塑性斜面の2D精度・一般3D収束・SolidCraft入力UI/途中ゲートは別の受入項目です。');q.write_text(v,'utf-8',newline='\n')
for name,a,b in [('docs/MC_VIRTUAL_WORK_BOUND.ja.md','（開発候補）',''),('docs/SRM_MEAN_STRESS_BOUND.ja.md','（開発中）',''),
    ('docs/SRM_CONTROLLED_SEARCH.ja.md','受理増分ごとの保存候補','受理増分ごとの保存'),('docs/SRM_SEARCH_JOBS.ja.md','保存・復旧できる乾燥MCのSRM探索候補','保存・復旧できる乾燥MCのSRM探索'),
    ('docs/SRM_SEARCH.ja.md','乾燥MCのSRM探索候補','乾燥MCのSRM探索')]:
    q=root/name;v=q.read_text('utf-8');assert v.count(a)>=1;q.write_text(v.replace(a,b,1),'utf-8',newline='\n')
q=root/'docs/SRM_MEAN_STRESS_BOUND.ja.md';v=q.read_text('utf-8');v=v.replace('後続の隔離した645ファイル候補では、物理荷重を結び付けた乾燥MCの限定的な探索へ接続して検証します。','本体656版では、物理荷重を結び付けた乾燥MCの限定探索へ接続し、同じ数値ソースの全347 Rust・219 Python試験が合格しました。');q.write_text(v,'utf-8',newline='\n')
q=root/'docs/P6_STATUS.ja.md';v=q.read_text('utf-8');lines=v.splitlines()
for i,line in enumerate(lines):
    if line.startswith('|P3c-01|'):lines[i]='|P3c-01|MC系SRMの独立履歴再生、Fs探索、途中工程ゲート|部分（探索・保存CLIを本体接続）|656版は各率の新規初期化・施工接頭部再生、要求率の実評価、受理増分の原子的保存・協調停止・新規復旧を実装。同一数値ソースの347 Rust/219 Pythonと実CLIが合格。完了状態と要求率の合否を区別し、不収束を上端にしない。無低減確定epochを起点にする探索との統合・途中工程の自動継続ゲート/Native UIは未完。|'
    if line.startswith('|P3c-02|'):lines[i]='|P3c-02|物理的上端候補と数値エラーの区別、2D/3D照合|部分（乾燥MCの独立排除）|独立6反力の平均応力と、分布拘束の全GP仮想仕事を接続。全点の有限支持・正重み・3Dテンソル・拘束と誤差余裕を検査。固定1200 Tet10の二次仮想場でF=1.45の離散釣合い排除を独立確認した。正確なFs・連続体安定性・元の2D変位5％基準は未達。一般3D/水/構造との統合・収束は残る。|'
q.write_text('\n'.join(lines)+'\n','utf-8',newline='\n')
current_sources=manifest();assert len(current_sources)==656 and current_sources.keys()==new.keys()
differences={n for n in new if new[n]!=current_sources[n]};assert differences and all(Path(n).suffix=='.md' for n in differences)
assert all(sha(candidate/n)==h for n,h in new.items())
dump(out/'adopted-source-manifest.json',current_sources)
dump(out/'adoption-source-equivalence.json',dict(passed=True,Root_files=656,qualified_files=656,
    every_nonMarkdown_source_byte_exact_to_full347_Rust219_Python_candidate=True,
    only_documentation_differences=sorted(differences),all_old635_sources_and_binary_deeply_backed_up=True,
    live_source_runtimes_and_solver_images_and_jobs_unchanged=True,
    actual_Root_CLI_build_and_qualification_pending=True,p6_complete=False))
print('Adopted qualified656 sources; actual Root build and CLI/Python validation required',flush=True)
