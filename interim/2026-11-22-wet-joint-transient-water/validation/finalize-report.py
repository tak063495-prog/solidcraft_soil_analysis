import hashlib,json,shutil,subprocess
from pathlib import Path
root=Path.cwd();h=root/'interim/2026-11-22-wet-joint-transient-water';vd=h/'validation'
read=lambda p:json.loads(Path(p).read_text())
def write(p,v):Path(p).write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
raw_monitor=read(h/'validation.json');assert raw_monitor['source_hashes_unchanged'] and raw_monitor['python_id_union_exact']
write(vd/'initial-driver-summary.json',raw_monitor)
retry=read(vd/'python-core-0-retry-process.json');assert retry['exit_code']==0 and retry['source_hashes_unchanged']
assert [x['name'] for x in raw_monitor['commands'] if x['exit_code']!=0]==['python-core-0']
v=dict(raw_monitor);v['commands']=[retry if x['name']=='python-core-0' else x for x in raw_monitor['commands']]
v['complete']=all(x['exit_code']==0 for x in v['commands'])
selected=[read(vd/g/'python-test-inventory.json') for g in ['python-core-0-retry','python-core-1','python-mean-evidence','python-spatial-evidence']]
assert all(x['exit_code']==0 for x in selected)
v['python_total']=sum(x['tests_run'] for x in selected)
v['python_group_totals']={x['group']:x['tests_run'] for x in selected}
for key,field in [('python_failures','failures'),('python_errors','errors'),('python_skipped','skipped'),('python_expected_failures','expected_failures'),('python_unexpected_successes','unexpected_successes')]:v[key]=sum(x[field] for x in selected)
v['initial_incomplete_attempts_retained']=[x for x in raw_monitor['commands'] if x['exit_code']!=0]
v['initial_driver_summary']='validation/initial-driver-summary.json'
v['recovery']='Two old archived ZIPs absent from the sparse materialization caused two errors in the initial core0. Unchanged original Git blobs were downloaded/SHA-verified; the complete 274-test group was rerun successfully in a new process. No test/source/tolerance change.'
assert v['complete'] and v['rust_total_passed']==451 and v['python_total']==554
assert all(v[k]==0 for k in ['rust_failed','rust_ignored','rust_filtered_out','python_failures','python_errors','python_skipped','python_expected_failures','python_unexpected_successes'])
write(vd/'driver-summary.json',v)
core=read(vd/'source-before.json');late={}
for p in vd.glob('supplemental*before.json'):late.update(read(p)['hashes'])
assert set(core).isdisjoint(late)
combined={**core,**late};assert all(sha(root/p)==digest for p,digest in combined.items())
write(vd/'supplemental-source-after.json',late)
base_ids=read(root/'interim/2026-11-20-wet-joint-plastic-bulk-integration/validation/python-core-0/python-test-inventory.json')['discovered_test_ids']
groups=['python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence'];invs=selected
ids=invs[0]['discovered_test_ids'];assert all(x['discovered_test_ids']==ids for x in invs)
assert set(base_ids)<=set(ids);new=sorted(set(ids)-set(base_ids));assert len(new)==11
write(h/'python-coverage.json',dict(baseline_tests=len(base_ids),baseline_ids_preserved=True,new_tests=len(new),new_test_ids=new,discovered_total=len(ids),executed_total=sum(x['tests_run'] for x in invs),all_ids_partitioned_once=True,groups={x['group']:x['tests_run'] for x in invs},tools_tests_included=True))
manifest=read(h/'archive-sha256.json');binding={}
for label,live in [('probe',root/'target/release/examples/wet_joint_transient_probe'),('solver',root/'target/release/soil-fem-3d')]:
 sealed=read(root/'work/c04/formal/identity.json')['binaries'][label];assert sha(live)==sealed==sha(root/'work/c04/formal'/label)
 binding[label]=dict(sha256=sealed,validation_build_matches_fresh_solve=True)
source_archive={p:d for p,d in combined.items() if 'source/'+p in manifest['entries']}
assert all(manifest['entries']['source/'+p]['sha256']==digest for p,digest in source_archive.items())
write(h/'source-binding.json',dict(core_sources_frozen_before_fmt_clippy_Rust=len(core),supplemental_sources_frozen_during_Rust_before_Python=len(late),total_sources_unchanged=len(combined),all_hashes_unchanged=True,hashes=combined,binaries=binding,archive_source_bindings_checked=len(source_archive),capture_scope='No claim the four later Python observation/series files were frozen before Rust began; all were frozen before complete Python discovery.',fresh_input_sha256={str(p.relative_to(root/'work/c04/formal')):sha(p) for p in (root/'work/c04/formal').rglob('case.json')}))
v.update(core_source_files=len(core),supplemental_source_files=len(late),final_numeric_test_observation_sources_unchanged=len(combined),baseline_python_ids_preserved=543,new_python_tests=11,source_binding='source-binding.json')
write(h/'validation.json',v)
shutil.copy2(root/'work/c04/required-evidence-download.json',vd/'required-evidence-download.json')
text=f'''\n\n## Linux全件検証・封印・Git保持の確定結果\n\nfmt、全ターゲットclippy（`--locked --all-targets -- -D warnings`）、locked release build/testはすべて終了コード0。**Rust451件（69集計群）、Python554件（既存543＋新規11、toolsの2件を内包）**が成功した。Python4群は{v['python_group_totals']}。発見/実行554 IDに欠落・重複0、failure/error/ignored/skipped/filtered-out/expected-failure/unexpected-successは0。新規11件の単独スモークを全件数へ重複加算していない。初回core0は旧ZIP2件の配置不足で2エラーだったため、元Git blob/SHAで取得後に274件全体を別PID/別フォルダーで再実行した。初回の失敗ログ・inventory・driver summaryは保持し、最終成功件数へ混ぜていない。試験・閾値・旧参照値の変更はない。\n\nfmt/Rust前の1,147ソースと、Rust検証中・全Python発見前に追加固定した4ソースの計**1,151ファイル**は最後までSHA一致した。4件がRust開始前から固定済みだったとは主張しない。封印内ソース・実行版と新規求解/最終検証の結び付けは `source-binding.json`。Linux環境、固定依存版、スレッド数と各実PID/コマンド/終了コード/実測時間は `validation/`。\n\n生証拠は **82,725,393 bytes・5 part・1,462エントリー**。全part・ZIP・各エントリーのサイズ/CRC/SHAを新規展開で照合し、封印ソースの隔離再計測で測定JSONと原因診断のオブジェクト完全一致を確認した。`sealed-replay-verification.json`。保存再監査を18新規経路・実checkpoint再開へ加算しない。\n\n起点5,707追跡パス中**5,704 Git blobを保持**し、変更した既存3件は指定残件文書への追記だけで元バイト列を接頭部として保持。旧interim・src・既存tools/試験/入力/CI・NC実装/設定・旧失敗/参照/hash/snapshotは不変。`base-preservation.json`。**C-04/A-04・C-01/C-02/C-03・NC-01・P6は未達を維持し、新先端Windows CIは未確定**とする。\n'''
with (h/'README.ja.md').open('ab') as f:f.write(text.encode())
extra='\n\nC-04追補のLinux全件検証はfmt/clippy/build/testが成功、Rust451件・Python554件（既存543＋新規11、tools2件を内包）。失敗・無視・skip・filter除外・expected failure・unexpected successは0。新先端Windows CIの成否は未確定で、C-04/A-04等の未達判定は維持する。\n'
for p in [root/'docs/P6_STATUS.ja.md',root/'docs/STRUCTURE_WATER_JOINT_REMAINDERS.ja.md',root/'interim/2026-11-16-issue-triage/todo.ja.md']:
 with p.open('ab') as f:f.write(extra.encode())
base='7ae0a2a669e0e1df93f3dd9791bbeda056334c3e';basepaths=subprocess.check_output(['git','ls-tree','-r','--name-only',base],text=True).splitlines();allowed=['docs/P6_STATUS.ja.md','docs/STRUCTURE_WATER_JOINT_REMAINDERS.ja.md','interim/2026-11-16-issue-triage/todo.ja.md'];changed=set(subprocess.check_output(['git','diff','--name-only',base],text=True).splitlines())&set(basepaths);assert changed==set(allowed)
prefix={}
for name in allowed:
 old=subprocess.check_output(['git','show',base+':'+name]);now=Path(name).read_bytes();assert now.startswith(old);prefix[name]=dict(old_sha256=hashlib.sha256(old).hexdigest(),old_bytes=len(old),original_prefix_exact=True,appended_bytes=len(now)-len(old))
write(h/'base-preservation.json',dict(starting_commit=base,base_tracked_paths=len(basepaths),unchanged_existing_Git_blobs=len(basepaths)-len(changed),only_existing_changes=sorted(changed),original_document_prefixes=prefix,src_tools_existing_tests_examples_CI_unchanged=True,old_interim_unchanged_except_todo_append=True,NC_implementation_configuration_unchanged=True,old_evidence_hash_snapshots_not_rewritten=True))
actual=sum(x['elapsed_seconds'] for x in raw_monitor['commands'])/60
with (h/'REPRODUCE.ja.md').open('ab') as f:f.write(f'\n\n今回Linuxの全件検証は初回Rustテストビルド・旧証拠の配置不足による初回群と再検証を含む監視プロセス全体で約{actual:.1f}分。Windowsで4群を直列に実行する場合は60〜120分を初期予算とし、CPU/ディスク/初回ビルド差を見込む。全件実測値はvalidation.jsonを参照。タイムアウトを設けて成功扱いにはしない。\n'.encode())
print(json.dumps(dict(Rust=v['rust_total_passed'],Python=v['python_total'],groups=v['python_group_totals'],source_files=len(combined),unchanged_base=len(basepaths)-3,verification_minutes=actual),ensure_ascii=False,indent=2))
