param([string]$SolverPath,[string]$OutputDirectory,[double]$MeshSizeM=4.0,[double]$K0=1.0,[double]$PcOverMean=4.0)
$ErrorActionPreference='Stop'
$projectRoot=Split-Path -Parent $PSScriptRoot
$python=Join-Path $projectRoot '.venv\Scripts\python.exe'
$tag=(Get-Date -Format 'yyyyMMdd-HHmmss')+'-'+[guid]::NewGuid().ToString('N').Substring(0,8)
$runDir=if ($OutputDirectory) {[IO.Path]::GetFullPath($OutputDirectory)} else {Join-Path $projectRoot "runs\geostatic-slope-pit-$tag"}
if (Test-Path -LiteralPath $runDir) {throw "Output exists: $runDir"}
New-Item -ItemType Directory -Path $runDir | Out-Null
Push-Location -LiteralPath $projectRoot
Start-Transcript -Path "$runDir\verification.log" | Out-Null
try {
    $template=Get-Content -LiteralPath 'examples\ground-geostatic-free-mcc-template.json' -Raw | ConvertFrom-Json
    $template | Add-Member -NotePropertyName excavation_bodies -NotePropertyValue @('Pit')
    $template | ConvertTo-Json -Depth 30 | Set-Content -LiteralPath "$runDir\template.json" -Encoding utf8
    $prepareArguments=@{Design='examples\geostatic-slope-pit.json';GroundFile="$runDir\template.json";PrepareOnly=$true;ElementOrder=2;MeshSizeM=$MeshSizeM;OutputDirectory="$runDir\prepared"}
    if ($SolverPath) {$prepareArguments.SolverPath=$SolverPath}
    & "$PSScriptRoot\Invoke-Ground.ps1" @prepareArguments
    $solver="$runDir\prepared\solver\soil-fem-3d.exe"
    & $solver mesh-identity "$runDir\prepared\mesh.json" "$runDir\mesh-identity.json"
    if ($LASTEXITCODE -ne 0) {throw 'Mesh identity failed'}
    & $python tools\make_planar_overburden_seed.py "$runDir\prepared\mesh.json" "$runDir\mesh-identity.json" "$runDir\template.json" "$runDir\ground.json" --top-at-origin-m 10 --gradient-x -0.3333333333333333 --k0-x $K0 --k0-y $K0 --pc-ratio $PcOverMean
    if ($LASTEXITCODE -ne 0) {throw 'Planar example seed generation failed'}
    & "$PSScriptRoot\Invoke-Stages.ps1" -PreparedModel "$runDir\prepared" -GroundFile "$runDir\ground.json" -StageFile 'examples\stages-geostatic-slope-pit.json' -SolverPath $solver -OutputDirectory "$runDir\analysis"
    Write-Host "General initial equilibrium and local excavation completed: $runDir"
} finally {Stop-Transcript | Out-Null; Pop-Location}
