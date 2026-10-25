param(
    [string]$GroundRun,
    [string]$PreparedModel,
    [string]$StageFile,
    [string]$Design,
    [string]$GroundFile,
    [switch]$Live,
    [string]$Connect = '127.0.0.1:7878',
    [double]$MeshSizeM = 2.0,
    [ValidateSet(1,2)][int]$ElementOrder = 1,
    [ValidateSet(3,4,5)][int]$QuadratureOrder = 3,
    [string]$SolverPath,
    [string]$OutputDirectory,
    [string]$PythonPath,
    [switch]$PrepareOnly
)
$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
$env:PATH = "$env:USERPROFILE\.cargo\bin;$env:PATH"
$python = if($PythonPath){(Resolve-Path -LiteralPath $PythonPath).Path}else{Join-Path $projectRoot '.venv\Scripts\python.exe'}
if ($PreparedModel -and ($GroundRun -or $Design -or $Live -or -not $GroundFile -or -not $StageFile)) { throw '-PreparedModel requires explicit GroundFile/StageFile and cannot be combined with GroundRun/Design/Live' }
if ($PreparedModel) { $PreparedModel=(Resolve-Path -LiteralPath $PreparedModel).Path }
if ($GroundRun -and -not $StageFile) { throw 'Specify -StageFile with -GroundRun: stage pressures are absolute and must match the initial load' }
if (-not $StageFile) { $StageFile = Join-Path $projectRoot 'examples\stages-small.json' }
if (-not $GroundRun -and -not $GroundFile) { $GroundFile = Join-Path $projectRoot 'examples\ground-staged-small.json' }
$StageFile = (Resolve-Path -LiteralPath $StageFile).Path
if ($SolverPath) { $SolverPath = (Resolve-Path -LiteralPath $SolverPath).Path }
if ($Design) { $Design = (Resolve-Path -LiteralPath $Design).Path }
if ($GroundFile) { $GroundFile = (Resolve-Path -LiteralPath $GroundFile).Path }
if ($GroundRun -and ($Design -or $GroundFile -or $Live)) { throw 'Use -GroundRun or CAD input options' }
if ($GroundRun) { $GroundRun = (Resolve-Path -LiteralPath $GroundRun).Path }
$tag = 'stages-' + (Get-Date -Format 'yyyyMMdd-HHmmss') + '-' + [guid]::NewGuid().ToString('N').Substring(0,8)
$runDir = if ($OutputDirectory) { [IO.Path]::GetFullPath($OutputDirectory) } else { Join-Path $projectRoot "runs\$tag" }
if (Test-Path -LiteralPath $runDir) { throw "Output exists: $runDir" }
New-Item -ItemType Directory -Path $runDir | Out-Null
Push-Location -LiteralPath $projectRoot
Start-Transcript -Path "$runDir\run.log" | Out-Null
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
    Copy-Item -LiteralPath $StageFile -Destination "$runDir\stages.json"
    if ($PreparedModel) {
        & $python tools\copy_prepared_model.py $PreparedModel $GroundFile "$runDir\ground"
        if ($LASTEXITCODE -ne 0) { throw 'Prepared geometry/input capture failed' }
        $meshPath = "$runDir\ground\mesh.json"; $groundPath = "$runDir\ground\ground.json"
    } elseif ($GroundRun) {
        $manifest = Get-Content -LiteralPath "$GroundRun\accepted-run.json" -Raw | ConvertFrom-Json
        Copy-Item -LiteralPath "$GroundRun\accepted-run.json" -Destination "$runDir\initial-input-package.json"
        foreach ($name in @('mesh.json','ground.json')) {
            $actual = (Get-FileHash -LiteralPath "$GroundRun\$name" -Algorithm SHA256).Hash.ToLower()
            if ($manifest.files_sha256.$name -ne $actual) { throw "Initial package hash mismatch: $name" }
            Copy-Item -LiteralPath "$GroundRun\$name" -Destination "$runDir\$name"
            if ((Get-FileHash -LiteralPath "$runDir\$name" -Algorithm SHA256).Hash.ToLower() -ne $actual) { throw "Initial file changed during capture: $name" }
        }
        $meshPath = "$runDir\mesh.json"; $groundPath = "$runDir\ground.json"
    } else {
        $arguments = @{ OutputDirectory = "$runDir\ground"; MeshSizeM = $MeshSizeM; Connect = $Connect; ElementOrder=$ElementOrder; QuadratureOrder=$QuadratureOrder }
        if ($Design) { $arguments.Design = $Design }; if ($GroundFile) { $arguments.GroundFile = $GroundFile }; if ($Live) { $arguments.Live = $true }
        if ($SolverPath) { $arguments.SolverPath = $SolverPath }
        if ($GroundFile -and (Get-Content -LiteralPath $GroundFile -Raw | ConvertFrom-Json).initialization.method -eq 'equilibrate') { $arguments.PrepareOnly=$true }
        & "$PSScriptRoot\Invoke-Ground.ps1" @arguments
        $meshPath = "$runDir\ground\mesh.json"; $groundPath = "$runDir\ground\ground.json"
    }
    if ($PrepareOnly) {
        if (-not $PreparedModel) { throw '-PrepareOnly requires an explicit prepared model' }
        $preparedRecord = @{
            schema='soil-fem-prepared-stage-input/1'
            equilibrium_verified=$false
            mesh='ground/mesh.json'
            ground='ground/ground.json'
            stages='stages.json'
            solver='solver/soil-fem-3d.exe'
            files_sha256=@{
                'ground/mesh.json'=(Get-FileHash -LiteralPath $meshPath -Algorithm SHA256).Hash.ToLower()
                'ground/ground.json'=(Get-FileHash -LiteralPath $groundPath -Algorithm SHA256).Hash.ToLower()
                'stages.json'=(Get-FileHash -LiteralPath "$runDir\stages.json" -Algorithm SHA256).Hash.ToLower()
                'solver/soil-fem-3d.exe'=(Get-FileHash -LiteralPath $SolverPath -Algorithm SHA256).Hash.ToLower()
            }
        }
        $preparedRecord | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath "$runDir\prepared-stage-input.json" -Encoding utf8NoBOM
        Write-Host "Staged soil FEM inputs prepared; equilibrium not evaluated: $runDir"
        return
    }
    & $SolverPath solve-stages $meshPath $groundPath "$runDir\stages.json" "$runDir\stage-results.json"
    if ($LASTEXITCODE -ne 0) { throw 'Staged FEM incomplete: inspect failed result and last_committed; no accepted export was created' }
    & $python tools\export_stages.py $meshPath "$runDir\stage-results.json" "$runDir\results"
    if ($LASTEXITCODE -ne 0) { throw 'Stage result export/verification failed' }
    $verifyArgs = @('tools\verify_stage_package.py',$meshPath,$groundPath,"$runDir\stage-results.json","$runDir\results")
    if ($SolverPath) { $verifyArgs += @('--solver',$SolverPath) }
    & $python @verifyArgs
    if ($LASTEXITCODE -ne 0) { throw 'Independent staged force/load verification failed' }
    Write-Host "Staged soil FEM completed: $runDir"
} finally { Stop-Transcript | Out-Null; Pop-Location }
