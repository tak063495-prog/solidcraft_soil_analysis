$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
$env:PATH = "$env:USERPROFILE\.cargo\bin;$env:PATH"
$python = Join-Path $projectRoot '.venv\Scripts\python.exe'
$shellExecutable = (Get-Process -Id $PID).Path
$tag = (Get-Date -Format 'yyyyMMdd-HHmmss') + '-' + [guid]::NewGuid().ToString('N').Substring(0,8)
$resultDir = Join-Path $projectRoot "runs\local-test-$tag"
New-Item -ItemType Directory -Path $resultDir | Out-Null
Push-Location -LiteralPath $projectRoot
Start-Transcript -Path "$resultDir\verification.log" | Out-Null
try {
    & $python tools\run_logged.py "$resultDir\format.log" cargo fmt --check
    if ($LASTEXITCODE -ne 0) { throw 'Rust formatting check failed' }
    & $python tools\run_logged.py "$resultDir\rust-tests.log" cargo test --release --locked -- --nocapture
    if ($LASTEXITCODE -ne 0) { throw 'Rust verification failed' }
    & $python tools\run_logged.py "$resultDir\lint.log" cargo clippy --locked --all-targets -- -D warnings
    if ($LASTEXITCODE -ne 0) { throw 'Rust lint failed' }
    & $python tools\run_logged.py "$resultDir\build.log" cargo build --release --locked
    if ($LASTEXITCODE -ne 0) { throw 'CLI build failed' }
    & $python tools\run_logged.py "$resultDir\cad-tests.log" $python -m unittest discover -s tests -p test_cad_excavation.py -v
    if ($LASTEXITCODE -ne 0) { throw 'CAD Boolean membership verification failed' }
    & $python tools\run_logged.py "$resultDir\cad.log" $shellExecutable -NoProfile -File "$PSScriptRoot\Invoke-LocalExcavation.ps1" -ElementOrder 1 -OutputDirectory "$resultDir\cad"
    if ($LASTEXITCODE -ne 0) { throw 'SolidCraft local-excavation run failed' }
    & $python tools\run_logged.py "$resultDir\independent-cad.log" $python tools\verify_local_excavation.py "$resultDir\cad"
    if ($LASTEXITCODE -ne 0) { throw 'Independent 3D local-excavation comparison failed' }
    & $python tools\run_logged.py "$resultDir\cad-fine.log" $shellExecutable -NoProfile -File "$PSScriptRoot\Invoke-LocalExcavation.ps1" -ElementOrder 1 -MeshSizeM 1.5 -OutputDirectory "$resultDir\cad-fine"
    if ($LASTEXITCODE -ne 0) { throw 'Refined SolidCraft local-excavation run failed' }
    & $python tools\run_logged.py "$resultDir\independent-cad-fine.log" $python tools\verify_local_excavation.py "$resultDir\cad-fine"
    if ($LASTEXITCODE -ne 0) { throw 'Refined-mesh independent comparison failed' }
    & $python tools\run_logged.py "$resultDir\spatial.log" $python tools\compare_local_meshes.py "$resultDir\cad" "$resultDir\cad-fine" "$resultDir\spatial-sensitivity.json"
    if ($LASTEXITCODE -ne 0) { throw 'Spatial sensitivity calculation failed' }
    & $python tools\run_logged.py "$resultDir\independent-plastic.log" $python tools\verify_stages.py --solver target\release\soil-fem-3d.exe --output-dir "$resultDir\independent-plastic" --small-strain
    if ($LASTEXITCODE -ne 0) { throw 'Independent plastic MCC regression failed' }
    Write-Host "Local excavation verification passed: $resultDir"
} finally { Stop-Transcript | Out-Null; Pop-Location }
