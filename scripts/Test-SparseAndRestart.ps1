$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
$env:PATH = "$env:USERPROFILE\.cargo\bin;$env:PATH"
$python = Join-Path $projectRoot '.venv\Scripts\python.exe'
$shellExecutable = (Get-Process -Id $PID).Path
$tag = (Get-Date -Format 'yyyyMMdd-HHmmss') + '-' + [guid]::NewGuid().ToString('N').Substring(0,8)
$resultDir = Join-Path $projectRoot "runs\sparse-restart-test-$tag"
New-Item -ItemType Directory -Path $resultDir | Out-Null
Push-Location -LiteralPath $projectRoot
Start-Transcript -Path "$resultDir\verification.log" | Out-Null
try {
    & $python tools\run_logged.py "$resultDir\format.log" cargo fmt --check
    if ($LASTEXITCODE -ne 0) { throw 'Formatting check failed' }
    & $python tools\run_logged.py "$resultDir\rust-tests.log" cargo test --release --locked -- --nocapture
    if ($LASTEXITCODE -ne 0) { throw 'Rust tests failed' }
    & $python tools\run_logged.py "$resultDir\lint.log" cargo clippy --locked --all-targets -- -D warnings
    if ($LASTEXITCODE -ne 0) { throw 'Lint failed' }
    & $python tools\run_logged.py "$resultDir\build.log" cargo build --release --locked
    if ($LASTEXITCODE -ne 0) { throw 'Release build failed' }
    & $python tools\run_logged.py "$resultDir\cad-tests.log" $python -m unittest discover -s tests -p test_cad_excavation.py -v
    if ($LASTEXITCODE -ne 0) { throw 'CAD Boolean tests failed' }
    & $python tools\run_logged.py "$resultDir\linear.log" $python tools\verify_sparse.py --output-dir "$resultDir\linear"
    if ($LASTEXITCODE -ne 0) { throw 'Sparse linear reference/scaling verification failed' }
    & $python tools\run_logged.py "$resultDir\cad.log" $shellExecutable -NoProfile -File "$PSScriptRoot\Invoke-LocalExcavation.ps1" -ElementOrder 1 -OutputDirectory "$resultDir\cad"
    if ($LASTEXITCODE -ne 0) { throw 'CSR local excavation failed' }
    & $python tools\run_logged.py "$resultDir\independent-3d.log" $python tools\verify_local_excavation.py "$resultDir\cad"
    if ($LASTEXITCODE -ne 0) { throw 'Independent 3D CSR comparison failed' }
    & $python tools\run_logged.py "$resultDir\cad-075.log" $shellExecutable -NoProfile -File "$PSScriptRoot\Invoke-LocalExcavation.ps1" -ElementOrder 1 -MeshSizeM 0.75 -OutputDirectory "$resultDir\cad-075"
    if ($LASTEXITCODE -ne 0) { throw 'Refined CSR plastic excavation failed' }
    & $python tools\run_logged.py "$resultDir\checkpoint.log" $python tools\verify_checkpoint.py --ground-run "$resultDir\cad\ground" --output-dir "$resultDir\checkpoint"
    if ($LASTEXITCODE -ne 0) { throw 'Separate-process checkpoint verification failed' }
    & $python tools\run_logged.py "$resultDir\independent-plastic.log" $python tools\verify_stages.py --solver target\release\soil-fem-3d.exe --output-dir "$resultDir\independent-plastic" --small-strain
    if ($LASTEXITCODE -ne 0) { throw 'Independent MCC regression failed' }
    $record = @{
        schema = 'soil-fem-sparse-restart-suite/1'
        scope = 'CSR/MCC/accepted-increment checkpoint integration; P6 overall remains incomplete'
        solver_sha256 = (Get-FileHash -LiteralPath 'target\release\soil-fem-3d.exe' -Algorithm SHA256).Hash.ToLower()
        cargo_lock_sha256 = (Get-FileHash -LiteralPath 'Cargo.lock' -Algorithm SHA256).Hash.ToLower()
        result_directory = $resultDir
        passed = $true
    }
    $record | ConvertTo-Json | Set-Content -LiteralPath "$resultDir\verification.json" -Encoding utf8
    Write-Host "Sparse FEM and checkpoint verification passed: $resultDir"
} finally { Stop-Transcript | Out-Null; Pop-Location }
