from pathlib import Path
import hashlib,json,shutil,sys
sys.path.insert(0,'tools')
from freeze_sources import manifest
root=Path('work/c04-evidence');root.mkdir()
def copy_tree(source,target):
 for p in sorted(source.rglob('*')):
  if not p.is_file() or '__pycache__' in p.parts:continue
  to=target/p.relative_to(source);to.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,to)
copy_tree(Path('work/c04-global-final'),root/'data/global-final')
copy_tree(Path('work/c04-global-controls'),root/'retained/preliminary-controls')
copy_tree(Path('work/c04-validation-preparation'),root/'retained/preparation')
base=Path('work/c04-old-sealed-replay/unpacked/data/final');target=root/'data/old-baseline'
for name in ('mesh.json','ground.json'):target.mkdir(parents=True,exist_ok=True);shutil.copy2(base/name,target/name)
for w in (8,16,32):
 for kind in ('native','independent'):
  src=base/f'{kind}-w{w}';dst=target/src.name;dst.mkdir()
  names=['fields.jsonl','case.json']+(['raw/constitutive.jsonl','raw/run.json'] if kind=='native' else ['result.json'])
  for name in names:
   p=dst/name;p.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src/name,p)
(root/'reports').mkdir()
shutil.copy2('work/c04-old-sealed-replay/verification.json',root/'reports/old-sealed-verification.json')
shutil.copy2('work/c04-reference-assessment.json',root/'reports/reference-assessment.json')
shutil.copy2('work/c04-global-contributions.json',root/'reports/global-contributions.json')
# Exact source bytes used by the production binary and the independent helpers.
for name,digest in manifest().items():
 p=Path(name);to=root/'source'/name;to.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,to)
 if hashlib.sha256(to.read_bytes()).hexdigest()!=digest:raise ValueError(name)
for folder in ('2026-11-20-wet-joint-plastic-bulk-integration','2026-11-22-wet-joint-transient-water','2026-12-19-c04-transient-water-cause-diagnosis'):
 for p in sorted((Path('interim')/folder).glob('*.py')):
  to=root/'source'/p;to.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,to)
print(root)
