"""Give the private source its real upstream Git base; never modify Root Git."""
from pathlib import Path
import datetime,hashlib,json,os,shutil,subprocess,sys,zipfile
p=Path(__file__).resolve().parent;root=p.parents[1];native=p/'native';out=p/'tracked-native-preparation-v1';out.mkdir()
read=lambda q:json.loads(q.read_text('utf-8'))
sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
def dump(q,v):
 with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2)
dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
assert read(p/'ordinary-typed-large-v3/terminal.json')['passed'] and read(p/'quality-v5/terminal.json')['passed']
manifest=read(p/'quality-v5/source72.json');assert all(sha(native/name)==h for name,h in manifest['files_sha256'].items())
assert native.resolve().parent==p.resolve() and not (native/'.git').exists()
parent_git={str(q):sha(q) for q in [root/'.git/HEAD',root/'.git/config'] if q.is_file()}
preserved=out/'preserved-reader';preserved.mkdir()
bindings=read(p/'ordinary-typed-large-v3/library-bindings.json')
for original,h in bindings.items():
 q=Path(original);assert sha(q)==h;shutil.copyfile(q,preserved/q.name);assert sha(preserved/q.name)==h
for name in ['ordinary-reader.exe','reader-receipt.json','terminal.json']:
 q=p/'ordinary-typed-large-v3'/name;shutil.copyfile(q,preserved/name);assert sha(q)==sha(preserved/name)
dump(preserved/'identities.json',{q.name:sha(q) for q in preserved.iterdir() if q.is_file()})
ui=native/'crates/ui-egui/src/soil_results.rs';s=ui.read_text('utf-8')
old='ui.collapsing("Material parameters",|ui| {';assert s.count(old)==1
ui.write_text(s.replace(old,'ui_kit::section(ui,"soil_epoch_material_parameters","Material parameters",false,|ui| {'),'utf-8',newline='\n')
def call(label,args):
 if (p/'pause-requested').exists():raise SystemExit(2)
 command=['git',*map(str,args)]
 with (out/(label+'.log')).open('x',encoding='utf-8') as log:
  child=subprocess.Popen(command,cwd=p,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
  dump(out/(label+'-child.json'),dict(pid=child.pid,parent_pid=os.getpid(),executable='git',arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
  code=child.wait()
 dump(out/(label+'-terminal.json'),dict(exit_code=code));assert code==0
call('init',['init',native])
assert (native/'.git').is_dir()
for key,value in [('core.autocrlf','false'),('core.filemode','false'),('core.symlinks','false')]:
 call('config-'+key.split('.')[-1],['-C',native,'config',key,value])
base=manifest['base_commit'];assert base=='edf57cec293dcea301201f1d5173b33e16351466'
call('fetch-base',['-C',native,'fetch','--depth=1','https://github.com/tangency-labs/solidcraft.git',base])
tree=subprocess.check_output(['git','-C',str(native),'ls-tree','-rz','--full-tree',base])
entries={}
for row in tree.split(b'\0'):
 if not row:continue
 info,name=row.split(b'\t',1);mode,kind,oid=info.decode().split();assert kind=='blob'
 entries[name.decode('utf-8')]=oid
exact=0;eol=[]
archive=p.parent/'native-mc-epoch-view-development-20261029-v1/native-base.zip'
with zipfile.ZipFile(archive) as z:
 archive_files={name for name in z.namelist() if not name.endswith('/')}
 assert set(entries)==archive_files,(set(entries)-archive_files,archive_files-set(entries))
 for name,oid in entries.items():
  data=z.read(name)
  def blob(raw):return hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
  if blob(data)==oid:exact+=1
  else:assert blob(data.replace(b'\r\n',b'\n'))==oid,name;eol.append(name)
dump(out/'upstream-base-verification.json',dict(passed=True,commit=base,files=len(entries),byte_exact_files=exact,line_ending_only_files=eol,full_base_tree_compared=True))
call('index-base',['-C',native,'read-tree',base])
call('detached-base',['-C',native,'update-ref','--no-deref','HEAD',base])
head=subprocess.check_output(['git','-C',str(native),'rev-parse','HEAD'],text=True).strip();assert head==base
top=subprocess.check_output(['git','-C',str(native),'rev-parse','--show-toplevel'],text=True).strip();assert Path(top).resolve()==native.resolve()
# Cached tarball build scripts watched only the old parent Git path. Refresh the
# scripts' timestamps so Cargo restamps the next build from this actual checkout.
for name in ['crates/io/build.rs','crates/ui-egui/build.rs']:os.utime(native/name,None)
assert parent_git=={q:sha(Path(q)) for q in parent_git}
assert all(sha(native/name)==h for name,h in manifest['files_sha256'].items() if name!='crates/ui-egui/src/soil_results.rs')
dump(out/'terminal.json',dict(passed=True,real_upstream_base_fetched_and_verified=True,head=base,
 private_native_git_only=True,Root_git_untouched=True,no_new_commit_created=True,
 no_working_source_file_checkout_or_reset_performed=True,working_soil_patch_remains_explicitly_uncommitted=True,
 old_reader_and_libraries_preserved=True,material_parameters_use_standard_ui_kit_section=True,
 unknown_source_cache_trust_guard_unchanged=True,next_build_must_restamp_io_and_UI=True,
 focused_cache_and_UI_gates_and_full_CI_still_required=True,p6_complete=False))
print('Prepared tracked private Native source; Root Git untouched, original cache guard unchanged',flush=True)
