param([switch]$SkipChecks,[string]$MccBaseline,[string]$OutputDirectory)
$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
$env:PATH = "$env:USERPROFILE\.cargo\bin;$env:PATH"
$python = Join-Path $projectRoot '.venv\Scripts\python.exe'
$shellExecutable = (Get-Process -Id $PID).Path
if ($MccBaseline) {$MccBaseline=(Resolve-Path -LiteralPath $MccBaseline).Path}
$tag = (Get-Date -Format 'yyyyMMdd-HHmmss') + '-' + [guid]::NewGuid().ToString('N').Substring(0,8)
$resultDir = if ($OutputDirectory) {[IO.Path]::GetFullPath($OutputDirectory)} else {Join-Path $projectRoot "runs\anchor-test-$tag"}
if (Test-Path -LiteralPath $resultDir) {throw "Output exists: $resultDir"}
New-Item -ItemType Directory -Path $resultDir | Out-Null
Push-Location -LiteralPath $projectRoot
Start-Transcript -Path "$resultDir\verification.log" | Out-Null
function Invoke-Logged([string]$Name,[string[]]$Arguments) {
    & $python tools\run_logged.py "$resultDir\$Name.log" @Arguments
    if ($LASTEXITCODE -ne 0) {throw "Verification step failed: $Name (retained $resultDir)"}
}
try {
    if (-not $SkipChecks) {
        Invoke-Logged format @('cargo','fmt','--check')
        Invoke-Logged rust-tests @('cargo','test','--release','--locked','--','--nocapture')
        Invoke-Logged lint @('cargo','clippy','--locked','--all-targets','--','-D','warnings')
        Invoke-Logged build @('cargo','build','--release','--locked')
    }
    New-Item -ItemType Directory -Path "$resultDir\solver" | Out-Null
    $solver = "$resultDir\solver\soil-fem-3d.exe"
    Copy-Item -LiteralPath 'target\release\soil-fem-3d.exe' -Destination $solver
    $digest = (Get-FileHash -LiteralPath $solver -Algorithm SHA256).Hash.ToLower()
    Invoke-Logged cad-tests @($python,'-m','unittest','discover','-s','tests','-p','test_cad_excavation.py','-v')
    Invoke-Logged cad @($shellExecutable,'-NoProfile','-File',"$PSScriptRoot\Invoke-Stages.ps1",'-Design','examples\geostatic-slope.json','-GroundFile','examples\ground-beam-slope.json','-StageFile','examples\stages-anchor-slope.json','-ElementOrder','2','-MeshSizeM','4','-SolverPath',$solver,'-OutputDirectory',"$resultDir\cad")
    Invoke-Logged independent @($python,'tools\verify_anchor_coupling.py',"$resultDir\cad","$resultDir\independent.json")
    Invoke-Logged history @($python,'tools\verify_beam_history.py','--solver',$solver,'--run',"$resultDir\cad",'--output-dir',"$resultDir\history")
    Invoke-Logged checkpoint @($python,'tools\verify_construction_checkpoint.py','--solver',$solver,'--run',"$resultDir\cad",'--output-dir',"$resultDir\checkpoint")
    Invoke-Logged mcc @($shellExecutable,'-NoProfile','-File',"$PSScriptRoot\Invoke-Stages.ps1",'-PreparedModel',"$resultDir\cad\ground",'-GroundFile','examples\ground-geostatic-mcc-slope.json','-StageFile','examples\stages-anchor-slope.json','-SolverPath',$solver,'-OutputDirectory',"$resultDir\mcc")
    Invoke-Logged mcc-history @($python,'tools\verify_beam_history.py','--solver',$solver,'--run',"$resultDir\mcc",'--output-dir',"$resultDir\mcc-history")
    Invoke-Logged mcc-checkpoint @($python,'tools\verify_construction_checkpoint.py','--solver',$solver,'--run',"$resultDir\mcc",'--output-dir',"$resultDir\mcc-checkpoint")
    $regression=@($shellExecutable,'-NoProfile','-File',"$PSScriptRoot\Test-Beams.ps1",'-SkipChecks','-OutputDirectory',"$resultDir\beam-regression")
    if ($MccBaseline) {$regression+=@('-MccBaseline',$MccBaseline)}
    Invoke-Logged beam-regression $regression
    if ((Get-FileHash -LiteralPath 'target\release\soil-fem-3d.exe' -Algorithm SHA256).Hash.ToLower() -ne $digest) {throw 'Build binary changed during suite'}
    @{
        schema='soil-fem-anchor-suite/1';passed=$true;solver_sha256=$digest
        rust_format_lint_build_included=(-not $SkipChecks.IsPresent)
        scope='Embedded tension-only anchors, free/bonded lengths, force-controlled jacking and lock-off, slack/reload, MC/MCC cutback, exact restart and beam/excavation regressions; perfect bond only'
    } | ConvertTo-Json | Set-Content -LiteralPath "$resultDir\verification.json" -Encoding utf8
    Write-Host "Anchor integration verification passed: $resultDir"
} finally {Stop-Transcript | Out-Null; Pop-Location}
if (-not $SkipChecks) {
    & $python "$projectRoot\tools\record_anchor_validation.py" --run $resultDir
    if ($LASTEXITCODE -ne 0) {throw "Anchor evidence binding failed: $resultDir"}
}
