"""Reproduce the complete final measurements from verified sealed observations."""
import hashlib,json,pathlib,sys
root=pathlib.Path.cwd();here=root/'interim/2026-12-25-b01-curved-tri6-cause-diagnosis';sys.path.insert(0,str(here))
from b01cause_evidence import verify
from b01cause_report import summarize
unsealed=root/'work/b01cause-unsealed';manifest=verify(here)
obj=summarize(unsealed/'b01cause-formal')
c=json.loads((unsealed/'b01cause-cli-controls/cli-controls-observations.json').read_text());audit=json.loads((unsealed/'b01cause-cli-controls/cli-exact-history-audit.json').read_text())
obj['production_cli_controls']={k:c[k] for k in ['solver_sha256','refusals','accepted_controls','all_original_contracts_confirmed','prior_all_members_sha256_verified']}
obj['production_cli_controls']['accepted_legacy_exact_comparisons']=[dict(name=r['name'],exact=r['original_run_equal_except_authenticated_binary_and_time']) for r in c['records'] if r['accepted']]
obj['same_input_repeat']=json.loads((unsealed/'b01cause-formal/same-input-repeat.json').read_text())
obj['production_cli_controls']['accepted_legacy_exact_comparisons_after_recursive_time_audit']=audit['records'];obj['production_cli_controls']['initial_false_comparison_retained']=True
raw=(json.dumps(obj,ensure_ascii=False,indent=2)+'\n').encode();target=here/'measurements.json';assert raw==target.read_bytes(),'full final measurements differ'
report=dict(all_parts_and_entries_sha256_verified=manifest['all_sha256_verified'],reproduced_final_measurements_byte_exact=True,measurements_sha256=hashlib.sha256(raw).hexdigest(),numerical_reruns=0)
(here/'validation/full-measurements-reproduction.json').write_text(json.dumps(report,indent=2)+'\n')
(here/'validation/executed-b01cause_reproduce_saved.py').write_bytes(pathlib.Path(__file__).read_bytes())
print(report)
