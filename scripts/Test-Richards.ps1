param([string]$OutputDirectory)
$ErrorActionPreference='Stop'
$projectRoot=Split-Path -Parent $PSScriptRoot
$env:PATH="$env:USERPROFILE\.cargo\bin;$env:PATH"
$python=Join-Path $projectRoot '.venv\Scripts\python.exe'
$shellExecutable=(Get-Process -Id $PID).Path
$tag=(Get-Date -Format 'yyyyMMdd-HHmmss')+'-'+[guid]::NewGuid().ToString('N').Substring(0,8)
$resultDir=if($OutputDirectory){[IO.Path]::GetFullPath($OutputDirectory)}else{Join-Path $projectRoot "runs\richards-acceptance-$tag"}
if(Test-Path -LiteralPath $resultDir){throw "Output exists: $resultDir"}
New-Item -ItemType Directory -Path $resultDir | Out-Null
Push-Location -LiteralPath $projectRoot
Start-Transcript -Path "$resultDir\verification.log" | Out-Null
function Invoke-Logged([string]$Name,[string[]]$Arguments) {
    & $python tools\run_logged.py "$resultDir\$Name.log" @Arguments
    if($LASTEXITCODE -ne 0){throw "Verification failed: $Name (retained $resultDir)"}
}
try {
    Invoke-Logged freeze @($python,'tools\freeze_sources.py',$resultDir)
    Invoke-Logged format @('cargo','fmt','--check')
    Invoke-Logged rust-tests @('cargo','test','--release','--locked')
    Invoke-Logged lint @('cargo','clippy','--locked','--all-targets','--','-D','warnings')
    Invoke-Logged build @('cargo','build','--release','--locked')
    New-Item -ItemType Directory -Path "$resultDir\solver" | Out-Null
    $solver="$resultDir\solver\soil-fem-3d.exe"
    Copy-Item -LiteralPath target\release\soil-fem-3d.exe -Destination $solver
    $digest=(Get-FileHash -LiteralPath $solver -Algorithm SHA256).Hash.ToLower()
    Invoke-Logged cad-tests @($python,'-m','unittest','discover','-s','tests','-p','test_cad_excavation.py','-v')
    Invoke-Logged split-tests @($python,'-m','unittest','discover','-s','tests','-p','test_joint_split.py','-v')
    Invoke-Logged decimal @($python,'tools\retention_decimal_reference.py','--verify','tests\data\retention-decimal-reference.csv','--output-dir',"$resultDir\decimal")
    Invoke-Logged cad-capture @($shellExecutable,'-NoProfile','-File',"$PSScriptRoot\Invoke-Ground.ps1",'-Design','examples\local-excavation-cylinder.json','-GroundFile','examples\ground-local.json','-PrepareOnly','-MeshSizeM','4','-ElementOrder','2','-SolverPath',$solver,'-OutputDirectory',"$resultDir\cad")
    Invoke-Logged cases @($python,'tools\verify_richards_cases.py','--solver',$solver,'--output-dir',"$resultDir\cases",'--cad-mesh',"$resultDir\cad\mesh.json")
    Invoke-Logged convergence @($python,'tools\verify_richards_convergence.py','--solver',$solver,'--output-dir',"$resultDir\convergence")
    Invoke-Logged legacy @($python,'tools\verify_seepage_legacy.py','--solver',$solver,'--output-dir',"$resultDir\legacy")
    if((Get-FileHash -LiteralPath target\release\soil-fem-3d.exe -Algorithm SHA256).Hash.ToLower() -ne $digest){throw 'Build binary changed during suite'}
    @{schema='soil-fem-richards-suite/1';passed=$true;solver_sha256=$digest} | ConvertTo-Json | Set-Content -LiteralPath "$resultDir\verification.json" -Encoding utf8
} finally {Stop-Transcript | Out-Null;Pop-Location}
& $python "$projectRoot\tools\record_richards_validation.py" --run $resultDir
if($LASTEXITCODE -ne 0){throw "Hydraulic evidence binding failed: $resultDir"}
& $python "$projectRoot\tools\verify_acceptance_packet.py" $resultDir --output "$resultDir\integrity.json"
if($LASTEXITCODE -ne 0){throw "Hydraulic packet integrity failed: $resultDir"}
