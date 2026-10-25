param(
    [string]$Design,[string]$GroundFile,[string]$StageFile,[string]$SplitFile,
    [switch]$Live,[string]$Connect='127.0.0.1:7878',
    [switch]$Mcc,[double]$MeshSizeM=1.0,[ValidateSet(1,2)][int]$ElementOrder=2,
    [string]$SolverPath,[string]$OutputDirectory
)
$ErrorActionPreference='Stop'
$projectRoot=Split-Path -Parent $PSScriptRoot
$env:PATH="$env:USERPROFILE\.cargo\bin;$env:PATH"
$python=Join-Path $projectRoot '.venv\Scripts\python.exe'
if($Live -and $Design){throw 'Choose Live or Design'}
if (-not $Design -and -not $Live) {$Design=Join-Path $projectRoot 'examples\joint-column.json'}
if (-not $GroundFile) {$GroundFile=Join-Path $projectRoot $(if($Mcc){'examples\ground-joint-mcc.json'}else{'examples\ground-joint-mc.json'})}
if (-not $StageFile) {$StageFile=Join-Path $projectRoot 'examples\stages-joint-column.json'}
if (-not $SplitFile) {$SplitFile=Join-Path $projectRoot 'examples\joint-column-split.json'}
$tag=(Get-Date -Format 'yyyyMMdd-HHmmss')+'-'+[guid]::NewGuid().ToString('N').Substring(0,8)
$runDir=if($OutputDirectory){[IO.Path]::GetFullPath($OutputDirectory)}else{Join-Path $projectRoot "runs\joints-$tag"}
if(Test-Path -LiteralPath $runDir){throw "Output exists: $runDir"}
New-Item -ItemType Directory -Path $runDir | Out-Null
Push-Location -LiteralPath $projectRoot
Start-Transcript -Path "$runDir\run.log" | Out-Null
try {
    if(-not $SolverPath){cargo build --release --locked;if($LASTEXITCODE -ne 0){throw 'Build failed'};$SolverPath=Join-Path $projectRoot 'target\release\soil-fem-3d.exe'}
    $SolverPath=(Resolve-Path -LiteralPath $SolverPath).Path
    $digest=(Get-FileHash -LiteralPath $SolverPath -Algorithm SHA256).Hash
    New-Item -ItemType Directory -Path "$runDir\solver" | Out-Null
    Copy-Item -LiteralPath $SolverPath -Destination "$runDir\solver\soil-fem-3d.exe"
    if((Get-FileHash -LiteralPath $SolverPath -Algorithm SHA256).Hash -ne $digest -or (Get-FileHash -LiteralPath "$runDir\solver\soil-fem-3d.exe" -Algorithm SHA256).Hash -ne $digest){throw 'Solver changed during job capture'}
    $SolverPath="$runDir\solver\soil-fem-3d.exe"
    $geometry=@{GroundFile=$GroundFile;MeshSizeM=$MeshSizeM;ElementOrder=$ElementOrder;SolverPath=$SolverPath;PrepareOnly=$true;OutputDirectory="$runDir\conformal";Connect=$Connect}
    if($Live){$geometry.Live=$true}else{$geometry.Design=$Design}
    & "$PSScriptRoot\Invoke-Ground.ps1" @geometry
    & $python tools\split_joint_mesh.py --prepared "$runDir\conformal" --spec $SplitFile --output "$runDir\split"
    if($LASTEXITCODE -ne 0){throw 'JOINT displacement-node split failed'}
    & "$PSScriptRoot\Invoke-Stages.ps1" -PreparedModel "$runDir\split" -GroundFile $GroundFile -StageFile $StageFile -SolverPath $SolverPath -OutputDirectory "$runDir\analysis"
    Write-Host "JOINT analysis and independent mechanics audit completed: $runDir"
} finally {Stop-Transcript | Out-Null;Pop-Location}
