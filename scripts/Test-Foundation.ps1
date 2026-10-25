$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
$env:PATH = "$env:USERPROFILE\.cargo\bin;$env:PATH"
$python = Join-Path $projectRoot '.venv\Scripts\python.exe'
$cli = Join-Path $projectRoot 'tools\solidcraft-0.20.0\solidcraft-0.20.0-windows-x64-portable\solidcraft-cli.exe'
$tag = (Get-Date -Format 'yyyyMMdd-HHmmss') + '-' + [guid]::NewGuid().ToString('N').Substring(0,8)
$resultDir = Join-Path $projectRoot "runs\foundation-test-$tag"
New-Item -ItemType Directory -Path $resultDir | Out-Null
Push-Location -LiteralPath $projectRoot
Start-Transcript -Path "$resultDir\verification.log" | Out-Null
try {
    cargo fmt --check
    if ($LASTEXITCODE -ne 0) { throw 'Rust formatting check failed' }
    cargo test --locked
    if ($LASTEXITCODE -ne 0) { throw 'Rust tests failed' }
    cargo clippy --locked --all-targets -- -D warnings
    if ($LASTEXITCODE -ne 0) { throw 'Rust lint failed' }
    & "$PSScriptRoot\Invoke-Ground.ps1" -OutputDirectory "$resultDir\ground"
    & $python tools\verify_mcc.py --solver target\release\soil-fem-3d.exe --output-dir "$resultDir\mcc"
    if ($LASTEXITCODE -ne 0) { throw 'Independent MCC reference check failed' }
    & $python tools\verify_cad_mapping.py --cli $cli --base-run "$resultDir\ground" --output-dir "$resultDir\cad-mapping"
    if ($LASTEXITCODE -ne 0) { throw 'CAD mapping verification failed' }
    Write-Host "Foundation verification passed: $resultDir"
} finally { Stop-Transcript | Out-Null; Pop-Location }
