param(
    [string]$OutputDirectory,
    [string]$SolverPath,
    [string]$CliPath,
    [string]$PythonPath,
    [double]$CoarseMeshSizeM = 3.0,
    [double]$FineMeshSizeM = 2.0
)
$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
$python = if ($PythonPath) { (Resolve-Path -LiteralPath $PythonPath).Path } else { Join-Path $projectRoot '.venv/Scripts/python.exe' }
$cli = if ($CliPath) { (Resolve-Path -LiteralPath $CliPath).Path } else { Join-Path $projectRoot 'tools/solidcraft-0.20.0/solidcraft-0.20.0-windows-x64-portable/solidcraft-cli.exe' }
$runDir = if ($OutputDirectory) { [IO.Path]::GetFullPath($OutputDirectory) } else { Join-Path $projectRoot ('runs/solid-pile-' + (Get-Date -Format 'yyyyMMdd-HHmmss') + '-' + [guid]::NewGuid().ToString('N').Substring(0,8)) }
if (Test-Path -LiteralPath $runDir) { throw "Output exists: $runDir" }
Push-Location -LiteralPath $projectRoot
try {
    if ($SolverPath) { $SolverPath = (Resolve-Path -LiteralPath $SolverPath).Path }
    else {
        $env:PATH = "$env:USERPROFILE/.cargo/bin;$env:PATH"
        cargo build --release --locked
        if ($LASTEXITCODE -ne 0) { throw 'Solver build failed' }
        $SolverPath = Join-Path $projectRoot 'target/release/soil-fem-3d.exe'
    }
    & $python tools/prepare_solid_pile.py $runDir --cli $cli --coarse-mesh-size-m $CoarseMeshSizeM --fine-mesh-size-m $FineMeshSizeM
    if ($LASTEXITCODE -ne 0) { throw 'SolidCraft capture or conformal pile preparation failed' }
    $solverDirectory = Join-Path $runDir 'solver'
    New-Item -ItemType Directory -Path $solverDirectory | Out-Null
    $solver = Join-Path $solverDirectory 'soil-fem-3d.exe'
    $sourceHash = (Get-FileHash -LiteralPath $SolverPath -Algorithm SHA256).Hash.ToLower()
    Copy-Item -LiteralPath $SolverPath -Destination $solver
    if ((Get-FileHash -LiteralPath $solver -Algorithm SHA256).Hash.ToLower() -ne $sourceHash -or (Get-FileHash -LiteralPath $SolverPath -Algorithm SHA256).Hash.ToLower() -ne $sourceHash) { throw 'Solver changed while capturing the job' }
    @{solver_sha256=$sourceHash} | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $runDir 'solver-binding.json') -Encoding utf8NoBOM
    & $python tools/freeze_sources.py $runDir
    if ($LASTEXITCODE -ne 0) { throw 'Source capture failed' }
    foreach ($label in @('coarse','fine')) {
        $mesh = Join-Path $runDir "$label/mesh.json"
        $result = Join-Path $runDir "$label/stage-results.json"
        Copy-Item -LiteralPath (Join-Path $runDir 'stages.json') -Destination (Join-Path $runDir "$label/stages.json")
        & $python tools/run_logged.py (Join-Path $runDir "$label/solver.log") $solver solve-stages $mesh (Join-Path $runDir 'ground.json') (Join-Path $runDir 'stages.json') $result
        if ($LASTEXITCODE -ne 0) { throw "Solid pile $label solve failed" }
        & $python tools/export_stages.py $mesh $result (Join-Path $runDir "$label/results")
        if ($LASTEXITCODE -ne 0) { throw "Solid pile $label export failed" }
        & $python tools/run_logged.py (Join-Path $runDir "$label/package-verification.log") $python tools/verify_stage_package.py $mesh (Join-Path $runDir 'ground.json') $result (Join-Path $runDir "$label/results") --solver $solver
        if ($LASTEXITCODE -ne 0) { throw "Solid pile $label force, raw GP and VTU audit failed" }
    }
    & $python tools/run_logged.py (Join-Path $runDir 'independent.log') $python tools/verify_solid_pile.py $runDir
    if ($LASTEXITCODE -ne 0) { throw 'Independent solid pile comparison failed' }
    & $python tools/run_logged.py (Join-Path $runDir 'history.log') $python tools/verify_solid_pile_history.py $runDir
    if ($LASTEXITCODE -ne 0) { throw 'Solid pile checkpoint comparison failed' }
    Write-Host "Finite solid-pile elastic control completed: $runDir"
    Write-Host 'See independent-verification.json and history/verification.json. Spatial convergence, capacity, interface slip and pile-row equivalence remain unaccepted.'
} finally { Pop-Location }
