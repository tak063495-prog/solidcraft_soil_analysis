param(
    [string]$Design,
    [string]$GroundFile,
    [switch]$Live,
    [string]$Connect = '127.0.0.1:7878',
    [double]$MeshSizeM = 2.0,
    [ValidateSet(1,2)][int]$ElementOrder = 1,
    [ValidateSet(3,4,5)][int]$QuadratureOrder = 3,
    [string]$SolverPath,
    [switch]$PrepareOnly,
    [string]$InitializationCaseFile,
    [string]$OutputDirectory
)
$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
$env:PATH = "$env:USERPROFILE\.cargo\bin;$env:PATH"
$python = Join-Path $projectRoot '.venv\Scripts\python.exe'
$cli = Join-Path $projectRoot 'tools\solidcraft-0.20.0\solidcraft-0.20.0-windows-x64-portable\solidcraft-cli.exe'
if ($Live -and $Design) { throw 'Choose -Live or -Design' }
if (-not $GroundFile) { $GroundFile = Join-Path $projectRoot 'examples\ground.json' }
if (-not $Live -and -not $Design) { $Design = Join-Path $projectRoot 'examples\layers.json' }
if ($Design) { $Design = (Resolve-Path -LiteralPath $Design).Path }
$GroundFile = (Resolve-Path -LiteralPath $GroundFile).Path
if ($InitializationCaseFile) { $InitializationCaseFile=(Resolve-Path -LiteralPath $InitializationCaseFile).Path }
$generalGround = (Get-Content -LiteralPath $GroundFile -Raw | ConvertFrom-Json).initialization.method -eq 'equilibrate'
if ($generalGround -and -not $PrepareOnly -and -not $InitializationCaseFile) { throw 'General initialization requires -InitializationCaseFile (supports/Newton settings), or use Invoke-Stages with a general GroundFile.' }
if ($SolverPath) { $SolverPath = (Resolve-Path -LiteralPath $SolverPath).Path }
$runName = 'ground-' + (Get-Date -Format 'yyyyMMdd-HHmmss') + '-' + [guid]::NewGuid().ToString('N').Substring(0,8)
$runDir = if ($OutputDirectory) { [IO.Path]::GetFullPath($OutputDirectory) } else { Join-Path $projectRoot "runs\$runName" }
if (Test-Path -LiteralPath $runDir) { throw "Output directory already exists: $runDir" }
New-Item -ItemType Directory -Path $runDir | Out-Null
Push-Location -LiteralPath $projectRoot
try {
    if (-not $SolverPath) {
        cargo build --release --locked
        if ($LASTEXITCODE -ne 0) { throw 'Solver build failed' }
        $SolverPath = Join-Path $projectRoot 'target\release\soil-fem-3d.exe'
    }
    $sourceHash = (Get-FileHash -LiteralPath $SolverPath -Algorithm SHA256).Hash
    $binaryDirectory = Join-Path $runDir 'solver'
    New-Item -ItemType Directory -Path $binaryDirectory | Out-Null
    $jobSolver = Join-Path $binaryDirectory 'soil-fem-3d.exe'
    Copy-Item -LiteralPath $SolverPath -Destination $jobSolver
    if ((Get-FileHash -LiteralPath $jobSolver -Algorithm SHA256).Hash -ne $sourceHash -or (Get-FileHash -LiteralPath $SolverPath -Algorithm SHA256).Hash -ne $sourceHash) { throw 'Solver changed during job capture' }
    $SolverPath = $jobSolver
    Copy-Item -LiteralPath $GroundFile -Destination "$runDir\ground.json"
    $cadArgs = @('tools\cad_layers.py','--cli',$cli,'--output-dir',$runDir,'--ground',"$runDir\ground.json")
    if ($Live) { $cadArgs += @('--connect',$Connect) } else { $cadArgs += @('--design',$Design) }
    & $python @cadArgs
    if ($LASTEXITCODE -ne 0) { throw 'CAD body snapshot failed' }
    & $python tools\mesh_layers.py "$runDir\bodies.json" "$runDir\mesh.json" --mesh-size-m $MeshSizeM --order $ElementOrder --quadrature-order $QuadratureOrder
    if ($LASTEXITCODE -ne 0) { throw 'Conformal layered mesh generation failed' }
    if ($PrepareOnly) {
        & $python tools\prepare_model_record.py $runDir
        if ($LASTEXITCODE -ne 0) { throw 'Prepared model inventory failed' }
        Write-Host "Geometry prepared; equilibrium has not been evaluated: $runDir"
        return
    }
    if ($generalGround) {
        Copy-Item -LiteralPath $InitializationCaseFile -Destination "$runDir\initialization-case-source.json"
        $initialCase = Get-Content -LiteralPath "$runDir\initialization-case-source.json" -Raw | ConvertFrom-Json
        $initialCase.stages = @()
        $initialCase | Add-Member -NotePropertyName verification_reject_once_stages -NotePropertyValue @() -Force
        $initialCase | ConvertTo-Json -Depth 100 | Set-Content -LiteralPath "$runDir\stages.json" -Encoding utf8
        & $SolverPath solve-stages "$runDir\mesh.json" "$runDir\ground.json" "$runDir\stages.json" "$runDir\initial-equilibrium.json"
        if ($LASTEXITCODE -ne 0) { throw 'Initial equilibration failed; inspect last committed state; no accepted initial state created' }
        & $python tools\export_stages.py "$runDir\mesh.json" "$runDir\initial-equilibrium.json" "$runDir\initial-results"
        if ($LASTEXITCODE -ne 0) { throw 'Initial equilibrium export failed' }
        & $python tools\verify_stage_package.py "$runDir\mesh.json" "$runDir\ground.json" "$runDir\initial-equilibrium.json" "$runDir\initial-results" --solver $SolverPath
        if ($LASTEXITCODE -ne 0) { throw 'Independent initial equilibrium verification failed' }
        & $python tools\prepare_model_record.py $runDir --accepted
        if ($LASTEXITCODE -ne 0) { throw 'Accepted initial equilibrium inventory failed' }
        Write-Host "General initial equilibrium accepted: $runDir"
        return
    }
    & $SolverPath initialize-ground "$runDir\mesh.json" "$runDir\ground.json" "$runDir\initial-state.json"
    if ($LASTEXITCODE -ne 0) { throw 'K0/MCC initialization failed; no accepted initial state' }
    & $python tools\export_initial_vtu.py "$runDir\mesh.json" "$runDir\initial-state.json" "$runDir\initial-state.vtu"
    if ($LASTEXITCODE -ne 0) { throw 'Initial state VTU export failed' }
    $verifyArgs = @('tools\verify_ground.py',$runDir)
    if ($SolverPath) { $verifyArgs += @('--solver',$SolverPath) }
    & $python @verifyArgs
    if ($LASTEXITCODE -ne 0) { throw 'Independent ground verification failed' }
    Write-Host "Ground state accepted: $runDir"
} finally { Pop-Location }
