param([Parameter(Mandatory=$true)][string]$OutputDir)
$ErrorActionPreference = 'Stop'
if (Test-Path $OutputDir) { throw 'Use a new output directory.' }
if (!(Test-Path Cargo.toml)) { throw 'Run from the repository root.' }
$Python = Join-Path (Get-Location) '.venv/Scripts/python.exe'
if (!(Test-Path $Python)) { throw 'Create .venv and install requirements-reference.lock.txt first.' }
& git merge-base --is-ancestor 9cd1b719e8dbb65f6d4dafc2f53df2df90679075 HEAD
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
    & $Python tools/verify_wet_structures.py --solver target/release/soil-fem-3d.exe --output (Join-Path $OutputDir 'soil') --water-mode soil
    if ($LASTEXITCODE -ne 0) { throw 'Limited SOIL integration failed; retain evidence without acceptance.' }
    Write-Host 'SOIL scope: elastic, closed sticking, inherited bottom head and continuum role only. NC remains unaccepted.'
} finally {
    Stop-Transcript
}
