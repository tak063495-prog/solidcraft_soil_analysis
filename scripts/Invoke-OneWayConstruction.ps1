param(
    [string]$PreparedModel,[string]$Design,[switch]$Live,[string]$Connect='127.0.0.1:7878',
    [Parameter(Mandatory=$true)][string]$GroundFile,[Parameter(Mandatory=$true)][string]$CaseFile,
    [double]$MeshSizeM=2.,[ValidateSet(1,2)][int]$ElementOrder=2,[ValidateSet(3,4,5)][int]$QuadratureOrder=3,
    [string]$SolverPath,[string]$OutputDirectory
)
$ErrorActionPreference='Stop'
$projectRoot=Split-Path -Parent $PSScriptRoot
$python=Join-Path $projectRoot '.venv\Scripts\python.exe'
$env:PATH="$env:USERPROFILE\.cargo\bin;$env:PATH"
if($PreparedModel -and ($Live -or $Design)){throw 'Choose a prepared model or a fresh CAD input'}
$GroundFile=(Resolve-Path -LiteralPath $GroundFile).Path
$CaseFile=(Resolve-Path -LiteralPath $CaseFile).Path
$tag=(Get-Date -Format 'yyyyMMdd-HHmmss')+'-'+[guid]::NewGuid().ToString('N').Substring(0,8)
$resultDir=if($OutputDirectory){[IO.Path]::GetFullPath($OutputDirectory)}else{Join-Path $projectRoot "runs\one-way-construction-$tag"}
if(Test-Path -LiteralPath $resultDir){throw "Output exists: $resultDir"}
New-Item -ItemType Directory -Path $resultDir | Out-Null
Push-Location -LiteralPath $projectRoot
try {
    if(-not $SolverPath){cargo build --release --locked;if($LASTEXITCODE -ne 0){throw 'Build failed'};$SolverPath=Join-Path $projectRoot 'target\release\soil-fem-3d.exe'}
    $SolverPath=(Resolve-Path -LiteralPath $SolverPath).Path
    New-Item -ItemType Directory -Path "$resultDir\solver" | Out-Null
    $digest=(Get-FileHash -LiteralPath $SolverPath -Algorithm SHA256).Hash
    Copy-Item -LiteralPath $SolverPath -Destination "$resultDir\solver\soil-fem-3d.exe"
    $solver="$resultDir\solver\soil-fem-3d.exe"
    if((Get-FileHash -LiteralPath $solver -Algorithm SHA256).Hash -ne $digest){throw 'Solver changed during capture'}
    Copy-Item -LiteralPath $GroundFile -Destination "$resultDir\ground.json"
    Copy-Item -LiteralPath $CaseFile -Destination "$resultDir\case.json"
    if($PreparedModel){Copy-Item -LiteralPath (Join-Path (Resolve-Path -LiteralPath $PreparedModel).Path 'mesh.json') -Destination "$resultDir\mesh.json"}
    else {
        $cadArgs=@{GroundFile="$resultDir\ground.json";PrepareOnly=$true;ElementOrder=$ElementOrder;QuadratureOrder=$QuadratureOrder;MeshSizeM=$MeshSizeM;SolverPath=$solver;OutputDirectory="$resultDir\cad";Connect=$Connect}
        if($Live){$cadArgs.Live=$true}else {if(-not $Design){throw 'Specify Design or Live'};$cadArgs.Design=$Design}
        & "$PSScriptRoot\Invoke-Ground.ps1" @cadArgs
        if($LASTEXITCODE -ne 0){throw 'CAD preparation failed'}
        Copy-Item -LiteralPath "$resultDir\cad\mesh.json" -Destination "$resultDir\mesh.json"
    }
    & $solver solve-one-way-construction "$resultDir\mesh.json" "$resultDir\ground.json" "$resultDir\case.json" "$resultDir\result.json"
    if($LASTEXITCODE -ne 0){throw "Incomplete construction; inspect last committed state: $resultDir"}
    & $python tools\export_one_way_construction.py "$resultDir\mesh.json" "$resultDir\result.json" "$resultDir\results" --solver $solver
    if($LASTEXITCODE -ne 0){throw 'Construction export verification failed'}
    Write-Host "ONE_WAY construction complete: $resultDir"
} finally {Pop-Location}
