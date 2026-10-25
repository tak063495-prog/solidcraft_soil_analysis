param([switch]$SkipChecks,[string]$OutputDirectory,[string]$MccBaseline)
$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
$env:PATH = "$env:USERPROFILE\.cargo\bin;$env:PATH"
$python = Join-Path $projectRoot '.venv\Scripts\python.exe'
$shellExecutable = (Get-Process -Id $PID).Path
if ($MccBaseline) {$MccBaseline=(Resolve-Path -LiteralPath $MccBaseline).Path}
$tag = (Get-Date -Format 'yyyyMMdd-HHmmss') + '-' + [guid]::NewGuid().ToString('N').Substring(0,8)
$resultDir = if ($OutputDirectory) {[IO.Path]::GetFullPath($OutputDirectory)} else {Join-Path $projectRoot "runs\local-acceptance-$tag"}
if (Test-Path -LiteralPath $resultDir) {throw "Output exists: $resultDir"}
New-Item -ItemType Directory -Path $resultDir | Out-Null
Push-Location -LiteralPath $projectRoot
Start-Transcript -Path "$resultDir\verification.log" | Out-Null
function Invoke-Logged([string]$Name,[string[]]$Arguments) {
    & $python tools\run_logged.py "$resultDir\$Name.log" @Arguments
    if ($LASTEXITCODE -ne 0) {throw "Verification step failed: $Name (retained $resultDir)"}
}
try {
    Invoke-Logged freeze @($python,'tools\freeze_sources.py',$resultDir)
    if (-not $SkipChecks) {
        Invoke-Logged format @('cargo','fmt','--check')
        Invoke-Logged rust-tests @('cargo','test','--release','--locked','--','--nocapture')
        Invoke-Logged lint @('cargo','clippy','--locked','--all-targets','--','-D','warnings')
        Invoke-Logged build @('cargo','build','--release','--locked')
    }
    New-Item -ItemType Directory -Path "$resultDir\solver" | Out-Null
    $solver = "$resultDir\solver\soil-fem-3d.exe"
    Copy-Item -LiteralPath 'target\release\soil-fem-3d.exe' -Destination $solver
    $solverHash = (Get-FileHash -LiteralPath $solver -Algorithm SHA256).Hash.ToLower()
    Invoke-Logged cad-tests @($python,'-m','unittest','discover','-s','tests','-p','test_cad_excavation.py','-v')
    Invoke-Logged flat @($shellExecutable,'-NoProfile','-File',"$PSScriptRoot\Invoke-LocalExcavation.ps1",'-MeshSizeM','4','-SolverPath',$solver,'-OutputDirectory',"$resultDir\flat")
    Invoke-Logged curved @($shellExecutable,'-NoProfile','-File',"$PSScriptRoot\Invoke-LocalExcavation.ps1",'-Design','examples\local-excavation-cylinder.json','-MeshSizeM','4','-SolverPath',$solver,'-OutputDirectory',"$resultDir\curved")
    if ($MccBaseline) {Invoke-Logged mcc-regression @($python,'tools\verify_mcc_regression.py',$MccBaseline,"$resultDir\flat\stage-results.json","$resultDir\mcc-regression.json")}
    Invoke-Logged independent-3d @($python,'tools\verify_quadratic_local.py','--solver',$solver,'--ground-run',"$resultDir\flat\ground",'--output-dir',"$resultDir\independent-3d")
    Invoke-Logged plastic-history @($python,'tools\verify_checkpoint.py','--solver',$solver,'--ground-run',"$resultDir\flat\ground",'--output-dir',"$resultDir\plastic-history")
    if ((Get-FileHash -LiteralPath 'target\release\soil-fem-3d.exe' -Algorithm SHA256).Hash.ToLower() -ne $solverHash) {throw 'Build binary changed during suite'}
    @{
        schema='soil-fem-local-excavation-acceptance/1';passed=$true;solver_sha256=$solverHash
        rust_format_lint_build_included=(-not $SkipChecks.IsPresent)
        scope='SolidCraft box/cylinder tools; Tet10 MCC local excavation; independent elastic-branch 3D reference; plastic cutback and fresh-process exact restart'
    } | ConvertTo-Json | Set-Content -LiteralPath "$resultDir\verification.json" -Encoding utf8
    Write-Host "Local excavation acceptance passed: $resultDir"
} finally {Stop-Transcript | Out-Null; Pop-Location}
if (-not $SkipChecks) {
    & $python "$projectRoot\tools\record_local_acceptance.py" --run $resultDir
    if ($LASTEXITCODE -ne 0) {throw "Acceptance evidence binding failed: $resultDir"}
    & $python "$projectRoot\tools\verify_acceptance_packet.py" $resultDir --output "$resultDir\integrity.json"
    if ($LASTEXITCODE -ne 0) {throw "Acceptance packet integrity failed: $resultDir"}
}
