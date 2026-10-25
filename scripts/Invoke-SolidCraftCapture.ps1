param(
    [Parameter(Mandatory=$true)][string]$CaptureFolder,
    [string]$SplitFile,
    [double]$MeshSizeM=2.0,
    [ValidateSet(1,2)][int]$ElementOrder=2,
    [ValidateSet(3,4,5)][int]$QuadratureOrder=3,
    [string]$SolverPath,
    [string]$OutputDirectory,
    [string]$PythonPath,
    [switch]$PrepareOnly
)
$ErrorActionPreference='Stop'
$projectRoot=Split-Path -Parent $PSScriptRoot
$python=if($PythonPath){(Resolve-Path -LiteralPath $PythonPath).Path}else{Join-Path $projectRoot '.venv\Scripts\python.exe'}
$CaptureFolder=(Resolve-Path -LiteralPath $CaptureFolder).Path
if(-not $SolverPath){$SolverPath=Join-Path $projectRoot 'target\release\soil-fem-3d.exe'}
$SolverPath=(Resolve-Path -LiteralPath $SolverPath).Path
if($SplitFile){$SplitFile=(Resolve-Path -LiteralPath $SplitFile).Path}
$tag='native-capture-'+(Get-Date -Format 'yyyyMMdd-HHmmss')+'-'+[guid]::NewGuid().ToString('N').Substring(0,8)
$runDir=if($OutputDirectory){[IO.Path]::GetFullPath($OutputDirectory)}else{Join-Path $projectRoot "runs\$tag"}
if(Test-Path -LiteralPath $runDir){throw "Output exists: $runDir"}
New-Item -ItemType Directory -Path $runDir | Out-Null
Push-Location -LiteralPath $projectRoot
Start-Transcript -Path "$runDir\run.log" | Out-Null
try {
    New-Item -ItemType Directory -Path "$runDir\solver" | Out-Null
    $digest=(Get-FileHash -LiteralPath $SolverPath -Algorithm SHA256).Hash
    Copy-Item -LiteralPath $SolverPath -Destination "$runDir\solver\soil-fem-3d.exe"
    if((Get-FileHash -LiteralPath $SolverPath -Algorithm SHA256).Hash -ne $digest -or (Get-FileHash -LiteralPath "$runDir\solver\soil-fem-3d.exe" -Algorithm SHA256).Hash -ne $digest){throw 'Solver changed during capture'}
    $SolverPath="$runDir\solver\soil-fem-3d.exe"
    & $python tools\import_solidcraft_capture.py $CaptureFolder "$runDir\conformal"
    if($LASTEXITCODE -ne 0){throw 'Native capture integrity/body import failed'}
    $stages=Get-Content -LiteralPath "$runDir\conformal\stages.json" -Raw | ConvertFrom-Json
    if($stages.joint_system -and -not $SplitFile){throw 'JOINT input requires an explicit SplitFile'}
    if($SplitFile -and -not $stages.joint_system){throw 'SplitFile requires a JOINT system in the captured stages'}
    & $python tools\mesh_layers.py "$runDir\conformal\bodies.json" "$runDir\conformal\mesh.json" --mesh-size-m $MeshSizeM --order $ElementOrder --quadrature-order $QuadratureOrder
    if($LASTEXITCODE -ne 0){throw 'Conformal meshing failed'}
    & $python tools\prepare_model_record.py "$runDir\conformal"
    if($LASTEXITCODE -ne 0){throw 'Prepared inventory failed'}
    $prepared="$runDir\conformal"
    if($SplitFile){
        & $python tools\split_joint_mesh.py --prepared $prepared --spec $SplitFile --output "$runDir\split"
        if($LASTEXITCODE -ne 0){throw 'JOINT displacement-node split failed'}
        $prepared="$runDir\split"
    }
    & "$PSScriptRoot\Invoke-Stages.ps1" -PreparedModel $prepared -GroundFile "$runDir\conformal\ground.json" -StageFile "$runDir\conformal\stages.json" -SolverPath $SolverPath -OutputDirectory "$runDir\analysis" -PythonPath $python -PrepareOnly:$PrepareOnly
    if($PrepareOnly){
        if(-not (Test-Path -LiteralPath "$runDir\analysis\prepared-stage-input.json")){throw 'Prepared-stage input receipt was not created'}
        Write-Host "Native capture prepared for owned staged analysis: $runDir"
        return
    }
    & $python tools\export_solidcraft_view.py $runDir "$runDir\view"
    if($LASTEXITCODE -ne 0){throw 'Native viewer data export failed; completed mechanical results remain in analysis'}
    & $python tools\verify_solidcraft_view.py $runDir "$runDir\view" "$runDir\view-verification.json"
    if($LASTEXITCODE -ne 0){throw 'Native viewer data/VTU comparison failed; inspect the retained job'}
    Write-Host "Native SolidCraft capture -> FEM and independent audit completed: $runDir"
} finally {Stop-Transcript | Out-Null;Pop-Location}
