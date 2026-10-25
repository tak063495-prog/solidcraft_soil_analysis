param([switch]$SkipChecks,[string]$MccBaseline,[string]$OutputDirectory,[string]$SolverPath)
$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
$env:PATH = "$env:USERPROFILE\.cargo\bin;$env:PATH"
$python = Join-Path $projectRoot '.venv\Scripts\python.exe'
$shellExecutable = (Get-Process -Id $PID).Path
if ($MccBaseline) { $MccBaseline=(Resolve-Path -LiteralPath $MccBaseline).Path }
$tag = (Get-Date -Format 'yyyyMMdd-HHmmss') + '-' + [guid]::NewGuid().ToString('N').Substring(0,8)
$resultDir = if ($OutputDirectory) {[IO.Path]::GetFullPath($OutputDirectory)} else {Join-Path $projectRoot "runs\geostatic-test-$tag"}
if (Test-Path -LiteralPath $resultDir) {throw "Output exists: $resultDir"}
if ($SolverPath) {$SolverPath=(Resolve-Path -LiteralPath $SolverPath).Path}
New-Item -ItemType Directory -Path $resultDir | Out-Null
Push-Location -LiteralPath $projectRoot
Start-Transcript -Path "$resultDir\verification.log" | Out-Null
function Invoke-Logged([string]$Name,[string[]]$Arguments) {
    & $python tools\run_logged.py "$resultDir\$Name.log" @Arguments
    if ($LASTEXITCODE -ne 0) { throw "Verification step failed: $Name" }
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
    if (-not $SolverPath) {$SolverPath=Join-Path $projectRoot 'target\release\soil-fem-3d.exe'}
    Copy-Item -LiteralPath $SolverPath -Destination $solver
    $hash = (Get-FileHash -LiteralPath $solver -Algorithm SHA256).Hash.ToLower()
    Invoke-Logged cad-tests @($python,'-m','unittest','discover','-s','tests','-p','test_cad_excavation.py','-v')
    Invoke-Logged mc-slope @($shellExecutable,'-NoProfile','-File',"$PSScriptRoot\Invoke-Stages.ps1",'-Design','examples\geostatic-slope.json','-GroundFile','examples\ground-geostatic-mc-slope.json','-StageFile','examples\stages-geostatic-slope.json','-ElementOrder','2','-MeshSizeM','4','-SolverPath',$solver,'-OutputDirectory',"$resultDir\mc-slope")
    Invoke-Logged independent-slope @($python,'tools\verify_geostatic_slope.py',"$resultDir\mc-slope","$resultDir\independent-slope.json")
    Invoke-Logged mcc-slope @($shellExecutable,'-NoProfile','-File',"$PSScriptRoot\Invoke-Stages.ps1",'-Design','examples\geostatic-slope.json','-GroundFile','examples\ground-geostatic-mcc-slope.json','-StageFile','examples\stages-geostatic-slope.json','-ElementOrder','2','-MeshSizeM','4','-SolverPath',$solver,'-OutputDirectory',"$resultDir\mcc-slope")
    Invoke-Logged identity @($solver,'mesh-identity',"$resultDir\mc-slope\ground\mesh.json","$resultDir\mesh-identity.json")
    Invoke-Logged seed @($python,'tools\make_planar_overburden_seed.py',"$resultDir\mc-slope\ground\mesh.json","$resultDir\mesh-identity.json",'examples\ground-geostatic-free-mcc-template.json',"$resultDir\free-mcc-input.json",'--top-at-origin-m','10','--gradient-x','-0.3333333333333333','--k0-x','0.5','--k0-y','0.5','--pc-ratio','1.5')
    Invoke-Logged free-mcc @($shellExecutable,'-NoProfile','-File',"$PSScriptRoot\Invoke-Stages.ps1",'-PreparedModel',"$resultDir\mc-slope\ground",'-GroundFile',"$resultDir\free-mcc-input.json",'-StageFile','examples\stages-geostatic-slope.json','-SolverPath',$solver,'-OutputDirectory',"$resultDir\free-mcc")
    Invoke-Logged checkpoint @($python,'tools\verify_geostatic_checkpoint.py','--solver',$solver,'--run',"$resultDir\free-mcc",'--output-dir',"$resultDir\checkpoint")
    Invoke-Logged slope-pit @($shellExecutable,'-NoProfile','-File',"$PSScriptRoot\Invoke-GeostaticSlopePit.ps1",'-SolverPath',$solver,'-OutputDirectory',"$resultDir\slope-pit")
    Invoke-Logged slope-pit-checkpoint @($python,'tools\verify_geostatic_checkpoint.py','--solver',$solver,'--run',"$resultDir\slope-pit\analysis",'--output-dir',"$resultDir\slope-pit-checkpoint")
    Invoke-Logged mc-pit @($shellExecutable,'-NoProfile','-File',"$PSScriptRoot\Invoke-Stages.ps1",'-PreparedModel',"$resultDir\slope-pit\prepared",'-GroundFile','examples\ground-geostatic-mc-slope-pit.json','-StageFile','examples\stages-geostatic-slope-pit.json','-SolverPath',$solver,'-OutputDirectory',"$resultDir\mc-pit")
    Invoke-Logged independent-pit @($python,'tools\verify_geostatic_slope.py',"$resultDir\mc-pit","$resultDir\independent-pit.json")
    Invoke-Logged initial-only @($shellExecutable,'-NoProfile','-File',"$PSScriptRoot\Invoke-Ground.ps1",'-Design','examples\geostatic-slope.json','-GroundFile','examples\ground-geostatic-mcc-slope.json','-InitializationCaseFile','examples\stages-geostatic-slope.json','-ElementOrder','2','-MeshSizeM','4','-SolverPath',$solver,'-OutputDirectory',"$resultDir\initial-only")
    Invoke-Logged mcc-regression-cad @($shellExecutable,'-NoProfile','-File',"$PSScriptRoot\Invoke-LocalExcavation.ps1",'-MeshSizeM','4','-SolverPath',$solver,'-OutputDirectory',"$resultDir\mcc-regression-cad")
    if ($MccBaseline) { Invoke-Logged mcc-compare @($python,'tools\verify_mcc_regression.py',$MccBaseline,"$resultDir\mcc-regression-cad\stage-results.json","$resultDir\mcc-regression.json") }
    $record = @{
        schema='soil-fem-geostatic-suite/1';passed=$true;solver_sha256=$hash;result_directory=$resultDir
        scope='General initial-state import, nonlinear equilibration/reference, slope comparison and exact multi-phase restart; P6 overall incomplete'
        rust_format_lint_build_included=(-not $SkipChecks.IsPresent)
        cargo_lock_sha256=(Get-FileHash -LiteralPath 'Cargo.lock' -Algorithm SHA256).Hash.ToLower()
    }
    $record | ConvertTo-Json | Set-Content -LiteralPath "$resultDir\verification.json" -Encoding utf8
    Write-Host "Geostatic integration verification passed: $resultDir"
} finally { Stop-Transcript | Out-Null; Pop-Location }
