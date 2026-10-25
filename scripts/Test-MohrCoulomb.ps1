param([switch]$SkipChecks,[string]$MccBaseline)
$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
$env:PATH = "$env:USERPROFILE\.cargo\bin;$env:PATH"
$python = Join-Path $projectRoot '.venv\Scripts\python.exe'
$shellExecutable = (Get-Process -Id $PID).Path
if ($MccBaseline) { $MccBaseline = (Resolve-Path -LiteralPath $MccBaseline).Path }
$tag = (Get-Date -Format 'yyyyMMdd-HHmmss') + '-' + [guid]::NewGuid().ToString('N').Substring(0,8)
$resultDir = Join-Path $projectRoot "runs\mc-test-$tag"
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
    $binaryDirectory = Join-Path $resultDir 'solver'
    New-Item -ItemType Directory -Path $binaryDirectory | Out-Null
    $solver = Join-Path $binaryDirectory 'soil-fem-3d.exe'
    Copy-Item -LiteralPath 'target\release\soil-fem-3d.exe' -Destination $solver
    $solverHash = (Get-FileHash -LiteralPath $solver -Algorithm SHA256).Hash.ToLower()
    Invoke-Logged cad-tests @($python,'-m','unittest','discover','-s','tests','-p','test_cad_excavation.py','-v')
    Invoke-Logged material-paths @($python,'tools\verify_mc.py','--solver',$solver,'--output-dir',"$resultDir\material-paths")
    Invoke-Logged columns @($python,'tools\verify_mc_fem.py','--solver',$solver,'--output-dir',"$resultDir\columns")
    Invoke-Logged mc-cad @($shellExecutable,'-NoProfile','-File',"$PSScriptRoot\Invoke-LocalExcavation.ps1",'-GroundFile','examples\ground-mc-local.json','-StageFile','examples\stages-mc-local.json','-MeshSizeM','4','-SolverPath',$solver,'-OutputDirectory',"$resultDir\mc-cad")
    Invoke-Logged mcc-cad @($shellExecutable,'-NoProfile','-File',"$PSScriptRoot\Invoke-LocalExcavation.ps1",'-MeshSizeM','4','-SolverPath',$solver,'-OutputDirectory',"$resultDir\mcc-cad")
    if ($MccBaseline) { Invoke-Logged mcc-regression @($python,'tools\verify_mcc_regression.py',$MccBaseline,"$resultDir\mcc-cad\stage-results.json","$resultDir\mcc-regression.json") }
    Invoke-Logged mcc-theory @($python,'tools\verify_stages.py','--solver',$solver,'--output-dir',"$resultDir\mcc-theory",'--small-strain')
    $record = @{
        schema = 'soil-fem-mc-suite/1'
        scope = '3D MC/Rankine/Davis material and global FEM, typed MCC/MC history, exact restart and MCC regression; P6 overall incomplete'
        solver_sha256 = $solverHash
        cargo_lock_sha256 = (Get-FileHash -LiteralPath 'Cargo.lock' -Algorithm SHA256).Hash.ToLower()
        rust_format_lint_build_included = (-not $SkipChecks.IsPresent)
        result_directory = $resultDir
        passed = $true
    }
    $record | ConvertTo-Json | Set-Content -LiteralPath "$resultDir\verification.json" -Encoding utf8
    Write-Host "MC integration verification passed: $resultDir"
} finally { Stop-Transcript | Out-Null; Pop-Location }
