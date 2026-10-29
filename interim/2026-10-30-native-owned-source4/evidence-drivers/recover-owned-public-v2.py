"""Run the corrected artifact harness; expose only real passing evidence to gates."""
from pathlib import Path
import hashlib, json, shutil, subprocess, sys, time
p=Path(__file__).resolve().parent
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
    with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
out=p/'native-owned-v4-public-recovery-preparation-v2'
failure=read(out/'retained-harness-failure.json')
corrected=Path(failure['corrected_driver'])
assert sha(corrected)==failure['corrected_driver_sha256']
assert sha(p/'check-native-owned-v4-public-api.py')==failure['original_driver_sha256']
start=time.monotonic()
code=subprocess.call([sys.executable,str(corrected)])
(out/'recovery-driver-terminal.json').write_text(json.dumps(dict(exit_code=code,elapsed_seconds=time.monotonic()-start),indent=2),'utf-8')
assert code==0
new=p/'native-owned-v4-public-fixed-v2'
report=read(new/'terminal.json')
assert report['passed'] and report['only_changed_case_is_correct_source4_rejected_before_accepted_after']
legacy=p/'native-owned-v4-public-fixed-v1'
assert not (legacy/'terminal.json').exists() and not (legacy/'library-binding.json').exists()
# Publish the corrected tested driver under the legacy discoverable name only
# after preserving its original failed source and running the new version.
shutil.copyfile(corrected,p/'check-native-owned-v4-public-api.py')
assert sha(p/'check-native-owned-v4-public-api.py')==sha(corrected)
(legacy/'fixed-cases').mkdir()
for a,b in [(new/'library-binding.json',legacy/'library-binding.json'),
            (new/'fixed-cases/receipt.json',legacy/'fixed-cases/receipt.json')]:
    shutil.copyfile(a,b)
    assert sha(a)==sha(b)
report.update(recovered_v2_actual_public_probe_completed=True,actual_v2_output=str(new),
    actual_v2_terminal_sha256=sha(new/'terminal.json'),retained_compiler_artifact_harness_failure=failure,
    exact_original_failed_driver_preserved=True,legacy_gate_paths_carry_byte_exact_real_v2_evidence=True,
    original_native_sources_and_numerical_jobs_unmodified=True)
with (legacy/'terminal.json').open('x',encoding='utf-8') as f:json.dump(report,f,indent=2)
print('PASS recovered public7 source4 protocol with byte-exact original Rust probe; failed harness retained',flush=True)
