"""Promote only the three fully tested worker files; preserve exact previous sources."""
from pathlib import Path
import hashlib,json,shutil,sys
p=Path(__file__).resolve().parent;root=p.parents[1];runtime=p/'runtime'
read=lambda q:json.loads(q.read_text('utf-8-sig'))
sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
assert read(p/'acceptance.json')['passed'] and read(p/'quality/terminal.json')['python_tests']==171
tested=read(p/'quality/source-snapshot.json');before_candidate=read(p/'source-snapshot.json')
assert all(sha(runtime/n)==h for n,h in tested.items())
names=['tools/native_job.py','tools/native_stage_profile.py','tests/test_native_job.py']
sys.path.insert(0,str(root/'tools'));from freeze_sources import manifest,freeze
before=manifest();assert len(before)==611 and set(before)==set(tested)
assert all(before[n]==h for n,h in tested.items() if n not in names and not n.endswith('.md'))
assert all(before[n]==before_candidate[n] for n in names)
preserved=p/'root-before';preserved.mkdir();freeze(preserved)
for n in names:
    target=preserved/n;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(root/n,target)
    shutil.copy2(runtime/n,root/n);assert sha(root/n)==tested[n]
after=manifest();assert all(after[n]==h for n,h in tested.items() if not n.endswith('.md'))
frozen=p/'root-after';frozen.mkdir();freeze(frozen)
with (p/'promotion.json').open('x') as f:json.dump(dict(promoted=True,only_non_markdown_changes=names,production_source_files=611,
    python_tests=171,all_production_non_markdown_exact_to_accepted_candidate=True,rust_sources_unchanged_from272_passed=True,
    actual_native_abrupt_loss_recovered_all8_exact=True,solver_binary_sha256=read(p/'acceptance.json')['solver_sha256'],
    original_sources_preserved=True,actual_root_worker_remaining_excavation_acceptance_pending=True,native_gui_recovery_acceptance_separate=True,p6_complete=False),f,indent=2)
print('PASS promoted exact3 worker files, full171 QA and actual abrupt-loss control; original sources preserved')
