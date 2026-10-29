from pathlib import Path
import datetime,hashlib,json,shutil
p=Path(__file__).resolve().parent
read=lambda q:json.loads(q.read_text('utf-8'))
sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
out=p/'gui-native67/revision4'
obs=read(out/'observation.json')
auto=(out/'uniform-stress-auto-expanded-tree.txt').read_text('utf-8')
manual=(out/'uniform-stress-manual-narrow-range-tree.txt').read_text('utf-8')
assert 'Auto legend expanded for nearly uniform data.' in auto
assert 'Auto legend expanded for nearly uniform data.' not in manual
for text in [auto,manual]:assert 'Range: 6.97234e1 to 6.97234e1 kPa' in text
assert all(str(v) in manual for v in obs['raw_and_manual_narrow_range_kpa'])
assert read(p/'ordinary-large-v4/terminal.json')['passed']
identity=read(p/'desktop-and-reader-v4/binary-identities.json')
preserved=p/'preserved-native67-ui-v4';preserved.mkdir();files={}
for key in ['desktop','engine','serde_json']:
    source=Path(identity[key]);assert sha(source)==identity[key+'_sha256']
    dest=preserved/source.name;shutil.copyfile(source,dest);assert sha(dest)==identity[key+'_sha256']
    files[key]=dict(path=str(dest),sha256=sha(dest))
(preserved/'identities.json').write_text(json.dumps(files,indent=2),'utf-8')
source=read(p/'quality-v4/source67.json')['files_sha256']
assert all(sha(p/'native'/n)==h for n,h in source.items())
capture=read(p/'actual-gui-native-v3/terminal.json')
assert all(sha(Path(q))==h for q,h in capture['design_and_capture_sha256'].items())
with (out/'verification.json').open('x',encoding='utf-8') as f:
    json.dump(dict(passed=True,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        binary=files['desktop'],source67_unchanged=True,
        actual_GUI_auto_uniform_and_manual_narrow_override=True,
        original_GP_range_retained=True,original_GUI_design_and_capture_unchanged=True,
        screenshots_sha256={q.name:sha(q) for q in out.glob('*.png')},
        other_GPU_and_full_CI_and_release_separate=True,p6_complete=False),f,indent=2)
print('PASS Native67-v4 actual GUI auto legend/manual override; immutable binary and libraries preserved')
