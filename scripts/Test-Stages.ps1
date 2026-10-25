$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
$env:PATH = "$env:USERPROFILE\.cargo\bin;$env:PATH"
$python = Join-Path $projectRoot '.venv\Scripts\python.exe'
$tag = (Get-Date -Format 'yyyyMMdd-HHmmss') + '-' + [guid]::NewGuid().ToString('N').Substring(0,8)
$resultDir = Join-Path $projectRoot "runs\stages-test-$tag"
New-Item -ItemType Directory -Path $resultDir | Out-Null
Push-Location -LiteralPath $projectRoot
Start-Transcript -Path "$resultDir\verification.log" | Out-Null
try {
    cargo fmt --check
    if ($LASTEXITCODE -ne 0) { throw 'Rust formatting check failed' }
    cargo test --release --locked -- --nocapture
    if ($LASTEXITCODE -ne 0) { throw 'Rust verification failed' }
    cargo clippy --locked --all-targets -- -D warnings
    if ($LASTEXITCODE -ne 0) { throw 'Rust lint failed' }
    cargo build --release --locked
    if ($LASTEXITCODE -ne 0) { throw 'CLI build failed' }
    & $python tools\verify_stages.py --solver target\release\soil-fem-3d.exe --output-dir "$resultDir\independent-small" --small-strain
    if ($LASTEXITCODE -ne 0) { throw 'Independent small-increment comparison failed' }
    & $python tools\verify_stages.py --solver target\release\soil-fem-3d.exe --output-dir "$resultDir\independent-stress-test"
    if ($LASTEXITCODE -ne 0) { throw 'Independent numerical stress-test comparison failed' }
    & "$PSScriptRoot\Invoke-Stages.ps1" -OutputDirectory "$resultDir\cad"
    Write-Host "Staged FEM verification passed: $resultDir"
} finally { Stop-Transcript | Out-Null; Pop-Location }
