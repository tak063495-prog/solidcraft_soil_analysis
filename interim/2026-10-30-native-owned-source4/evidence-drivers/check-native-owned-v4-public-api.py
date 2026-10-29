"""Repeat the exact original public process probe against the fixed engine."""
from pathlib import Path
import datetime, hashlib, json, os, shutil, subprocess, sys, time

p = Path(__file__).resolve().parent
root = p.parents[1]
native = p / 'native-owned-v4'
out = p / 'native-owned-v4-public-fixed-v2'
out.mkdir()
read = lambda q: json.loads(q.read_text('utf-8-sig'))

def sha(q):
    with q.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def dump(q, value):
    with q.open('x', encoding='utf-8') as f:
        json.dump(value, f, indent=2, allow_nan=False)

def pause_check():
    if (p / 'pause-requested').exists():
        raise SystemExit(2)

pause_check()
dump(out / 'driver-child.json', dict(pid=os.getpid(), executable=sys._base_executable,
    script=str(Path(__file__).resolve()), created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
baseline_folder = p / 'native-owned-v4-public-repro-v1'
baseline = read(baseline_folder / 'terminal.json')
assert baseline['passed'] and baseline['baseline_gap_reproduced']
assert sha(p / 'owned-v4-completion-probe.rs') == baseline['probe_sha256']
quality = p / 'native-owned-v4-quality-v1'
gate = quality / 'terminal.json'
while not gate.exists():
    pause_check()
    time.sleep(5)
assert read(gate)['passed']
sources = read(quality / 'source-manifest.json')['files_sha256']
assert len(sources) == 73
binding = {str(native / name): h for name, h in sources.items()}
binding.update({str(q): sha(q) for q in [p / 'owned-v4-completion-probe.rs',
    p / 'owned-v4-protocol-test-worker.py', gate, quality / 'production-library.log', Path(__file__)]})
artifacts = []
for line in (quality / 'production-library.log').read_text('utf-8').splitlines():
    try:
        item = json.loads(line)
    except json.JSONDecodeError:
        continue
    if item.get('reason') == 'compiler-artifact':
        artifacts.append(item)

def libs(name):
    return {f: item for item in artifacts if item['target']['name'] == name
            for f in item['filenames'] if f.endswith('.rlib')}

engines = libs('solidcraft_engine')
assert len(engines) == 1
engine_artifact = Path(next(iter(engines)))
# A direct cargo build --lib reports the public, un-hashed root rlib. Bind its
# actual dependency rlib by full bytes, then use that exact Cargo fingerprint.
matches = [q for q in (native / 'target/ci/deps').glob('libsolidcraft_engine-*.rlib')
           if q.stat().st_size == engine_artifact.stat().st_size and sha(q) == sha(engine_artifact)]
assert len(matches) == 1, 'unique byte-identical engine dependency artifact required'
engine = matches[0]
binding[str(engine_artifact)] = sha(engine_artifact)
fingerprint = native / 'target/ci/.fingerprint' / ('solidcraft-engine-'+engine.stem.rsplit('-',1)[1]) / 'lib-solidcraft_engine.json'
deps = [row for row in read(fingerprint)['deps'] if row[1] == 'serde_json']
assert len(deps) == 1
expected = deps[0][-1]
matched = []
choices = []
for name, item in libs('serde_json').items():
    library = Path(name)
    folder = native / 'target/ci/.fingerprint' / ('serde_json-'+library.stem.rsplit('-',1)[1])
    code = folder / 'lib-serde_json'
    value = int.from_bytes(bytes.fromhex(code.read_text('ascii').strip()), 'little')
    choices.append(dict(path=name, features=item['features'], fingerprint_u64=str(value), matches_engine=value == expected))
    if value == expected:
        matched.append(library)
    binding.update({str(q): sha(q) for q in [library, code, folder / 'lib-serde_json.json']})
assert len(matched) == 1
serde = matched[0]
binding.update({str(q): sha(q) for q in [engine, fingerprint]})
document = p / 'actual14909-blocked-v2/workflow/conformal/source.solidcraft.json'
binding[str(document)] = sha(document)
assert read(baseline_folder / 'original-binding.json')[str(document)] == sha(document)
host = out / 'explicit-protocol-only-host'
(host / 'tools').mkdir(parents=True)
(host / '.venv/Scripts').mkdir(parents=True)
for a, b in [(baseline_folder / 'explicit-protocol-only-host/.venv/Scripts/python.exe', host / '.venv/Scripts/python.exe'),
             (baseline_folder / 'explicit-protocol-only-host/.venv/pyvenv.cfg', host / '.venv/pyvenv.cfg'),
             (baseline_folder / 'explicit-protocol-only-host/tools/native_job.py', host / 'tools/native_job.py')]:
    shutil.copyfile(a, b)
    assert sha(a) == sha(b)
    binding[str(a)] = sha(a)
    binding[str(b)] = sha(b)
dump(host / 'TEST-ONLY-NO-FEM.json', dict(protocol_test_only=True,
    no_numerical_solve_or_real_result_or_independent_numerical_audit=True,
    same_worker_bytes_as_original_baseline=True))
dump(out / 'library-binding.json', dict(passed=True, engine=str(engine), engine_sha256=sha(engine),
    serde=str(serde), serde_sha256=sha(serde), engine_serde_fingerprint_u64=str(expected), choices=choices,
    selected_by_actual_engine_dependency=True, source_manifest_sha256=sha(quality/'source-manifest.json')))
dump(out / 'source-binding.json', binding)
binary = out / 'owned-v4-completion-probe.exe'
env = dict(os.environ, PYTHONUTF8='1', PYTHONDONTWRITEBYTECODE='1', OPENBLAS_NUM_THREADS='1',
    MKL_NUM_THREADS='1', OMP_NUM_THREADS='1')
env.pop('PYTHONPATH', None)
env.pop('PYTHONHOME', None)
records = []

def call(label, command):
    pause_check()
    assert all(sha(Path(q)) == h for q, h in binding.items())
    start = time.monotonic()
    with (out / (label+'.log')).open('x', encoding='utf-8') as log:
        child = subprocess.Popen(list(map(str, command)), cwd=native, env=env, stdout=log,
            stderr=subprocess.STDOUT, creationflags=subprocess.CREATE_NO_WINDOW | subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        dump(out / (label+'-child.json'), dict(pid=child.pid, parent_pid=os.getpid(), executable=str(command[0]),
            arguments=list(map(str, command[1:])), created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        code = child.wait()
    assert all(sha(Path(q)) == h for q, h in binding.items())
    row = dict(step=label, exit_code=code, elapsed_seconds=time.monotonic()-start, sources_unchanged=True)
    records.append(row)
    dump(out / (label+'-terminal.json'), row)
    assert code == 0, label
    print(json.dumps(row), flush=True)

call('compile', ['C:/Users/link_/.cargo/bin/rustc.exe', '--edition', '2024',
    '-L', 'dependency='+str(serde.parent), '--extern', 'solidcraft_engine='+str(engine),
    '--extern', 'serde_json='+str(serde), p / 'owned-v4-completion-probe.rs', '-o', binary])
call('fixed', [binary, document, host, out / 'fixed-cases', 'fixed'])
report = read(out / 'fixed-cases/receipt.json')
assert report['passed'] and report['fixed_behavior_expected']
old_report = read(baseline_folder / 'baseline-cases/receipt.json')
assert report['original_projected_CAD_snapshot'] == old_report['original_projected_CAD_snapshot']
before = {r['mode']: r['observed_status']['result_usable'] for r in old_report['cases']}
after = {r['mode']: r['observed_status']['result_usable'] for r in report['cases']}
assert len(after) == 7 and set(after) == set(before)
assert [name for name in before if before[name] != after[name]] == ['source4']
assert after['legacy'] and after['source4'] and all(not after[name] for name in after if name not in ['legacy','source4'])
dump(out / 'terminal.json', dict(passed=True, records=records, source_files=73,
    same_exact_public_probe_source_and_same_original_projected_CAD_document_and_same_worker_bytes=True,
    real_owned_child_public_start_status_and_existing_receipt_hash_checks_exercised=True,
    only_changed_case_is_correct_source4_rejected_before_accepted_after=True,
    original_legacy3_acceptance_and_all5_invalid_or_future_refusals_preserved=True,
    source_manifest_sha256=sha(quality/'source-manifest.json'),
    numerical_solve_and_actual_large_reader_GUI_and_full_CI_release_and_P6_not_tested_here=True,
    explicit_protocol_stub_only=True, p6_complete=False))
print('PASS same public owned-process probe: correct source4 fixed, legacy and5 refusals retained; no FEM computed', flush=True)
