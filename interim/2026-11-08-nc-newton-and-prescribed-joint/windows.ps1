param([Parameter(Mandatory=$true)][string]$OutputDir)
$ErrorActionPreference = 'Stop'
if (Test-Path $OutputDir) { throw 'Use a new output directory.' }
if (!(Test-Path Cargo.toml)) { throw 'Run from the repository root.' }
$Python = Join-Path (Get-Location) '.venv/Scripts/python.exe'
if (!(Test-Path $Python)) { throw 'Create .venv and install requirements-reference.lock.txt first.' }
& git merge-base --is-ancestor 813b7517d1b810c795330491a67e0c71976cdc1c HEAD
if ($LASTEXITCODE -ne 0) { throw 'HEAD does not contain the specified PR #13 starting commit.' }
New-Item -ItemType Directory $OutputDir | Out-Null
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$env:CARGO_BUILD_JOBS = '2'
$env:RUST_TEST_THREADS = '2'
$env:PYTHONUTF8 = '1'
Start-Transcript -Path (Join-Path $OutputDir 'windows-transcript.log')
try {
    & git rev-parse HEAD
    & rustc --version --verbose
    & $Python -m pip freeze
    & cargo fmt --all --check
    if ($LASTEXITCODE -ne 0) { throw 'fmt failed' }
    & cargo clippy --locked --all-targets -- -D warnings
    if ($LASTEXITCODE -ne 0) { throw 'clippy failed' }
    & cargo build --locked --release
    if ($LASTEXITCODE -ne 0) { throw 'release build failed' }
    & cargo build --locked --release --example nc_newton_trace --example nc_material_probe
    if ($LASTEXITCODE -ne 0) { throw 'diagnostic build failed' }
    & cargo test --locked --release --all-targets
    if ($LASTEXITCODE -ne 0) { throw 'Rust tests failed' }
    & $Python -m unittest discover -s tests -p 'test_*.py' -v
    if ($LASTEXITCODE -ne 0) { throw 'Python tests failed' }
    & $Python -m unittest discover -s tools -p 'test_*.py' -v
    if ($LASTEXITCODE -ne 0) { throw 'Python tools tests failed' }
    $env:NC_NEW_INPUT = Join-Path $OutputDir 'nc-input'
    @'
import os, sys
from pathlib import Path
sys.path.insert(0, 'tools')
import probe_low_preload_nc as low
out = Path(os.environ['NC_NEW_INPUT']); out.mkdir(exist_ok=False)
inputs = low.expected_inputs(1/128, minimum_increment=1/4194304,
                            dense_lu=True, tangent_perturbation=1e-9)
for name, value in zip(low.NAMES, inputs):
    low.write(out/(name+'.json'), value)
'@ | & $Python -
    if ($LASTEXITCODE -ne 0) { throw '5 kPa input generation failed' }
    @{ passed = $false; state = 'started'; complete = $null } | ConvertTo-Json | Set-Content -Encoding utf8 (Join-Path $OutputDir 'nc-status.json')
    & target/release/examples/nc_newton_trace.exe $env:NC_NEW_INPUT (Join-Path $OutputDir 'nc-trace')
    if ($LASTEXITCODE -ne 0) { throw 'diagnostic process failed; retain partial trace without acceptance' }
    $NcRun = Get-Content -Raw (Join-Path $OutputDir 'nc-trace/run.json') | ConvertFrom-Json
    @{ passed = $false; state = 'diagnostic_only'; complete = $NcRun.complete; progress = $NcRun.last_committed.stage_progress; failure_reason = $NcRun.failure_reason } | ConvertTo-Json | Set-Content -Encoding utf8 (Join-Path $OutputDir 'nc-status.json')
    & $Python tools/measure_nc_newton_trace.py --trace (Join-Path $OutputDir 'nc-trace') --inputs $env:NC_NEW_INPUT --material-probe target/release/examples/nc_material_probe.exe --output (Join-Path $OutputDir 'nc-measured')
    if ($LASTEXITCODE -ne 0) { throw 'trace evidence incomplete; no NC acceptance' }
    & $Python tools/verify_wet_structures.py --solver target/release/soil-fem-3d.exe --output (Join-Path $OutputDir 'prescribed') --water-mode prescribed
    if ($LASTEXITCODE -ne 0) { throw 'limited PRESCRIBED reference failed' }
    Write-Host 'Read measurements.json: NC passed=false remains. PRESCRIBED scope is limited; SOIL remains unverified.'
} finally {
    Stop-Transcript
}
