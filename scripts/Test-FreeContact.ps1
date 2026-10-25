param([string]$OutputDirectory,[string]$ReuseRustChecks)
$ErrorActionPreference='Stop'
$projectRoot=Split-Path -Parent $PSScriptRoot
$env:PATH="$env:USERPROFILE\.cargo\bin;$env:PATH"
$python=Join-Path $projectRoot '.venv\Scripts\python.exe'
$shellExecutable=(Get-Process -Id $PID).Path
$tag=(Get-Date -Format 'yyyyMMdd-HHmmss')+'-'+[guid]::NewGuid().ToString('N').Substring(0,8)
$resultDir=if($OutputDirectory){[IO.Path]::GetFullPath($OutputDirectory)}else{Join-Path $projectRoot "runs\free-contact-acceptance-$tag"}
if(Test-Path -LiteralPath $resultDir){throw "Output exists: $resultDir"}
New-Item -ItemType Directory -Path $resultDir | Out-Null
Push-Location -LiteralPath $projectRoot
Start-Transcript -Path "$resultDir\verification.log" | Out-Null
function Invoke-Logged([string]$Name,[string[]]$Arguments) {
    & $python tools\run_logged.py "$resultDir\$Name.log" @Arguments
    if($LASTEXITCODE -ne 0){throw "Verification step failed: $Name (retained $resultDir)"}
}
try {
    Invoke-Logged freeze @($python,'tools\freeze_sources.py',$resultDir)
    Invoke-Logged format @('cargo','fmt','--check')
    if($ReuseRustChecks){
        Invoke-Logged reuse-rust-tests @($python,'tools\reuse_rust_test_evidence.py',$ReuseRustChecks,$resultDir)
    }else{
        Invoke-Logged rust-tests @('cargo','test','--release','--locked','--','--nocapture')
    }
    Invoke-Logged lint @('cargo','clippy','--locked','--all-targets','--','-D','warnings')
    Invoke-Logged build @('cargo','build','--release','--locked')
    New-Item -ItemType Directory -Path "$resultDir\solver" | Out-Null
    $solver="$resultDir\solver\soil-fem-3d.exe"
    Copy-Item -LiteralPath target\release\soil-fem-3d.exe -Destination $solver
    $digest=(Get-FileHash -LiteralPath $solver -Algorithm SHA256).Hash.ToLower()
    if($ReuseRustChecks -and (Get-Content -Raw -LiteralPath "$resultDir\rust-test-reuse.json" | ConvertFrom-Json).expected_solver_sha256 -ne $digest){throw 'Reused Rust evidence does not match the newly built binary'}
    Invoke-Logged cad-tests @($python,'-m','unittest','discover','-s','tests','-p','test_cad_excavation.py','-v')
    Invoke-Logged split-tests @($python,'-m','unittest','discover','-s','tests','-p','test_joint_split.py','-v')
    Invoke-Logged cad @($shellExecutable,'-NoProfile','-File',"$PSScriptRoot\Invoke-Joints.ps1",'-GroundFile','examples\ground-free-contact-mc.json','-StageFile','examples\stages-free-contact.json','-SolverPath',$solver,'-OutputDirectory',"$resultDir\cad")
    Invoke-Logged free-contact @($python,'-u','tools\verify_free_contact.py','--solver',$solver,'--output-dir',"$resultDir\free",'--cad-run',"$resultDir\cad\analysis")
    foreach($material in @('mc','mcc')) {
        Invoke-Logged $material @($python,'tools\run_joint_fixture.py','--solver',$solver,'--ground',"examples\ground-joint-$material.json",'--stages','examples\stages-joint-verification.json','--output',"$resultDir\$material")
    }
    Invoke-Logged compatibility @($python,'tools\verify_joint_compatibility.py','runs\joint-test-20261024-045026-efdb1ae7',$resultDir,"$resultDir\contact-state-compatibility.json")
    foreach($name in @('flat','curved')) {
        $design=if($name -eq 'flat'){'examples\local-excavation.json'}else{'examples\local-excavation-cylinder.json'}
        Invoke-Logged "local-$name" @($shellExecutable,'-NoProfile','-File',"$PSScriptRoot\Invoke-LocalExcavation.ps1",'-Design',$design,'-MeshSizeM','4','-SolverPath',$solver,'-OutputDirectory',"$resultDir\local-$name")
        Invoke-Logged "local-$name-regression" @($python,'tools\verify_local_same_input.py','--baseline','runs\local-acceptance-20261024-033221-8ada5d9c','--current',"$resultDir\local-$name",'--output-dir',"$resultDir\local-$name-same-input",'--comparison',"$resultDir\local-$name-regression.json")
    }
    if((Get-FileHash -LiteralPath target\release\soil-fem-3d.exe -Algorithm SHA256).Hash.ToLower() -ne $digest){throw 'Build binary changed during suite'}
    @{schema='soil-fem-free-contact-suite/1';passed=$true;solver_sha256=$digest} | ConvertTo-Json | Set-Content -LiteralPath "$resultDir\verification.json" -Encoding utf8
} finally {Stop-Transcript | Out-Null;Pop-Location}
& $python "$projectRoot\tools\record_free_contact.py" --run $resultDir
if($LASTEXITCODE -ne 0){throw "Free contact evidence binding failed: $resultDir"}
& $python "$projectRoot\tools\verify_acceptance_packet.py" $resultDir --output "$resultDir\integrity.json"
if($LASTEXITCODE -ne 0){throw "Acceptance packet integrity failed: $resultDir"}
