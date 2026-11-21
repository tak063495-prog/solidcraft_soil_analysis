param([Parameter(Mandatory=$true)][string]$OutputDir)
$ErrorActionPreference = 'Stop'
if (Test-Path $OutputDir) { throw 'Use a new output directory.' }
if (!(Test-Path Cargo.toml)) { throw 'Run from the repository root.' }
$Python = Join-Path (Get-Location) '.venv/Scripts/python.exe'
if (!(Test-Path $Python)) { throw 'Create .venv and install requirements-reference.lock.txt first.' }
& git merge-base --is-ancestor 9a15c83cd7de73e42e9344a8b4bfab31918e5d49 HEAD
if ($LASTEXITCODE -ne 0) { throw 'The specified PR #13 starting commit is missing.' }
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
    & cargo test --locked --release --all-targets
    if ($LASTEXITCODE -ne 0) { throw 'Rust tests failed' }
    & $Python -m unittest discover -s tests -p 'test_*.py' -v
    if ($LASTEXITCODE -ne 0) { throw 'Python tests failed' }
    & $Python -m unittest discover -s tools -p 'test_*.py' -v
    if ($LASTEXITCODE -ne 0) { throw 'Python tools tests failed' }
    & $Python tools/verify_wet_structures.py --solver target/release/soil-fem-3d.exe --output (Join-Path $OutputDir 'fixed-prescribed') --water-mode prescribed-fixed
    if ($LASTEXITCODE -ne 0) { throw 'Limited fixed PRESCRIBED integration failed; retain evidence without acceptance.' }
    Write-Host 'Fixed PRESCRIBED scope: elastic, closed sticking, inherited 1 m reservoir; incompatible rise omitted. NC remains unaccepted.'
    $Receipt = @{
        schema = 'soil-fem-fixed-prescribed-windows-validation/1'
        head_commit = (& git rev-parse HEAD)
        operating_system = [Environment]::OSVersion.VersionString
        software_checks_passed = $true
        fixed_prescribed_report = 'fixed-prescribed/verification.json'
        nc_passed = $false
        p3b_p4_rows_complete = $false
        p6_complete = $false
    }
    $Receipt | ConvertTo-Json -Depth 8 | Set-Content -Encoding UTF8 (Join-Path $OutputDir 'windows-validation.json')
} finally {
    Stop-Transcript
}
