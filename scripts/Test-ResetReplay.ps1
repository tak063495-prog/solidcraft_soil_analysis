param([switch]$SkipChecks,[string]$MccBaseline)
$ErrorActionPreference='Stop'
$projectRoot=Split-Path -Parent $PSScriptRoot
$env:PATH="$env:USERPROFILE\.cargo\bin;$env:PATH"
$python=Join-Path $projectRoot '.venv\Scripts\python.exe'
$shellExecutable=(Get-Process -Id $PID).Path
if ($MccBaseline) {$MccBaseline=(Resolve-Path -LiteralPath $MccBaseline).Path}
$tag=(Get-Date -Format 'yyyyMMdd-HHmmss')+'-'+[guid]::NewGuid().ToString('N').Substring(0,8)
$resultDir=Join-Path $projectRoot "runs\reset-replay-test-$tag"
New-Item -ItemType Directory -Path $resultDir | Out-Null
Push-Location -LiteralPath $projectRoot
Start-Transcript -Path "$resultDir\verification.log" | Out-Null
function Invoke-Logged([string]$Name,[string[]]$Arguments) {
    & $python tools\run_logged.py "$resultDir\$Name.log" @Arguments
    if ($LASTEXITCODE -ne 0) {throw "Verification step failed: $Name"}
}
try {
    if (-not $SkipChecks) {
        Invoke-Logged format @('cargo','fmt','--check')
        Invoke-Logged rust-tests @('cargo','test','--release','--locked','--','--nocapture')
        Invoke-Logged lint @('cargo','clippy','--locked','--all-targets','--','-D','warnings')
        Invoke-Logged build @('cargo','build','--release','--locked')
    }
    New-Item -ItemType Directory -Path "$resultDir\solver" | Out-Null
    $solver="$resultDir\solver\soil-fem-3d.exe";Copy-Item -LiteralPath 'target\release\soil-fem-3d.exe' -Destination $solver
    $hash=(Get-FileHash -LiteralPath $solver -Algorithm SHA256).Hash.ToLower()
    Invoke-Logged cad @($shellExecutable,'-NoProfile','-File',"$PSScriptRoot\Invoke-Stages.ps1",'-Design','examples\construction-fill.json','-GroundFile','examples\ground-construction-fill.json','-StageFile','examples\stages-reset-replay.json','-ElementOrder','2','-MeshSizeM','4','-SolverPath',$solver,'-OutputDirectory',"$resultDir\cad")
    Invoke-Logged independent @($python,'tools\verify_construction_linear.py',"$resultDir\cad","$resultDir\independent.json")
    Invoke-Logged checkpoint @($python,'tools\verify_construction_checkpoint.py','--solver',$solver,'--run',"$resultDir\cad",'--output-dir',"$resultDir\checkpoint")
    Invoke-Logged replay @($shellExecutable,'-NoProfile','-File',"$PSScriptRoot\Invoke-Replay.ps1",'-SourceRun',"$resultDir\cad",'-SolverPath',$solver,'-OutputDirectory',"$resultDir\replay")
    Invoke-Logged replay-verify @($python,'tools\verify_replay_batch.py',"$resultDir\replay","$resultDir\replay-verification.json")
    Invoke-Logged failures @($python,'tools\verify_replay_failures.py','--solver',$solver,'--run',"$resultDir\cad",'--output-dir',"$resultDir\failures")
    $regression=@($shellExecutable,'-NoProfile','-File',"$PSScriptRoot\Test-Geostatic.ps1",'-SkipChecks','-SolverPath',$solver,'-OutputDirectory',"$resultDir\geostatic-regression")
    if ($MccBaseline) {$regression+=@('-MccBaseline',$MccBaseline)}
    Invoke-Logged regression $regression
    @{schema='soil-fem-reset-replay-suite/1';passed=$true;solver_sha256=$hash;scope='Display/state RESET and original-input MC construction replay; deformed-reference 2D RESET, Fs search and P6 remain incomplete';rust_format_lint_build_included=(-not $SkipChecks.IsPresent)} | ConvertTo-Json | Set-Content -LiteralPath "$resultDir\verification.json" -Encoding utf8
    Write-Host "Reset/replay verification passed: $resultDir"
} finally {Stop-Transcript | Out-Null;Pop-Location}
