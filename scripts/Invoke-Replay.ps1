param([Parameter(Mandatory)][string]$SourceRun,[string]$SpecificationFile,[string]$SolverPath,[string]$OutputDirectory)
$ErrorActionPreference='Stop'
$projectRoot=Split-Path -Parent $PSScriptRoot
$env:PATH="$env:USERPROFILE\.cargo\bin;$env:PATH"
$python=Join-Path $projectRoot '.venv\Scripts\python.exe'
$SourceRun=(Resolve-Path -LiteralPath $SourceRun).Path
if (-not $SpecificationFile) {$SpecificationFile=Join-Path $projectRoot 'examples\replay-mc.json'}
$SpecificationFile=(Resolve-Path -LiteralPath $SpecificationFile).Path
if ($SolverPath) {$SolverPath=(Resolve-Path -LiteralPath $SolverPath).Path}
$tag=(Get-Date -Format 'yyyyMMdd-HHmmss')+'-'+[guid]::NewGuid().ToString('N').Substring(0,8)
$runDir=if ($OutputDirectory) {[IO.Path]::GetFullPath($OutputDirectory)} else {Join-Path $projectRoot "runs\replay-$tag"}
if (Test-Path -LiteralPath $runDir) {throw "Output exists: $runDir"}
$manifestPath=Join-Path $SourceRun 'accepted-stage-run.json'
$manifest=Get-Content -LiteralPath $manifestPath -Raw | ConvertFrom-Json
if ($manifest.schema -ne 'soil-fem-accepted-staged-run/1') {throw 'An accepted source stage package is required'}
$sources=@{}
$sourceDigests=@{}
foreach ($item in $manifest.inputs_and_results_sha256.PSObject.Properties) {
    $path=(Resolve-Path -LiteralPath $item.Name).Path
    if (-not $path.StartsWith($SourceRun+[IO.Path]::DirectorySeparatorChar,[StringComparison]::OrdinalIgnoreCase)) {throw 'Source manifest points outside source run'}
    if ((Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash.ToLower() -ne $item.Value) {throw "Source changed: $path"}
    $name=[IO.Path]::GetFileName($path)
    if ($sources.ContainsKey($name)) {throw "Duplicate source file: $name"}
    $sources[$name]=$path
    $sourceDigests[$name]=$item.Value
}
foreach ($name in @('mesh.json','ground.json','stages.json','stage-results.json')) {if (-not $sources.ContainsKey($name)) {throw "Source package missing $name"}}
New-Item -ItemType Directory -Path $runDir | Out-Null
Push-Location -LiteralPath $projectRoot
Start-Transcript -Path "$runDir\run.log" | Out-Null
try {
    if (-not $SolverPath) {cargo build --release --locked;if ($LASTEXITCODE -ne 0) {throw 'Build failed'};$SolverPath=Join-Path $projectRoot 'target\release\soil-fem-3d.exe'}
    New-Item -ItemType Directory -Path "$runDir\solver" | Out-Null
    $sourceHash=(Get-FileHash -LiteralPath $SolverPath -Algorithm SHA256).Hash
    $jobSolver="$runDir\solver\soil-fem-3d.exe"
    Copy-Item -LiteralPath $SolverPath -Destination $jobSolver
    if ((Get-FileHash -LiteralPath $jobSolver -Algorithm SHA256).Hash -ne $sourceHash -or (Get-FileHash -LiteralPath $SolverPath -Algorithm SHA256).Hash -ne $sourceHash) {throw 'Solver changed during capture'}
    Copy-Item -LiteralPath $manifestPath -Destination "$runDir\source-input-package.json"
    foreach ($name in @('mesh.json','ground.json','stages.json','stage-results.json')) {
        $destination=if ($name -eq 'stage-results.json') {"$runDir\source-stage-results.json"} else {"$runDir\$name"}
        Copy-Item -LiteralPath $sources[$name] -Destination $destination
        if ((Get-FileHash -LiteralPath $destination -Algorithm SHA256).Hash.ToLower() -ne $sourceDigests[$name]) {throw "Source changed during capture: $name"}
    }
    Copy-Item -LiteralPath $SpecificationFile -Destination "$runDir\specification.json"
    & $jobSolver replay-stages "$runDir\mesh.json" "$runDir\ground.json" "$runDir\stages.json" "$runDir\specification.json" "$runDir\replay-results.json"
    $replayExit=$LASTEXITCODE
    if (-not (Test-Path -LiteralPath "$runDir\replay-results.json")) {throw 'Replay input rejected; inspect run log'}
    & $python tools\export_replays.py "$runDir\mesh.json" "$runDir\replay-results.json" "$runDir\trials" --solver $jobSolver
    if ($LASTEXITCODE -ne 0) {throw 'Replay export/force verification failed'}
    if ($replayExit -ne 0) {throw 'Some replay trials did not complete; diagnostics retained; no stability/Fs decision'}
    Write-Host "Independent replays captured: $runDir"
} finally {Stop-Transcript | Out-Null;Pop-Location}
