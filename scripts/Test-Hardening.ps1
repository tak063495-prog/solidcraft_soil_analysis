param([string]$OutputDirectory,[string]$ReuseRustChecks)
$ErrorActionPreference='Stop'
$projectRoot=Split-Path -Parent $PSScriptRoot
$env:PATH="$env:USERPROFILE\.cargo\bin;$env:PATH"
$python=Join-Path $projectRoot '.venv\Scripts\python.exe'
$shellExecutable=(Get-Process -Id $PID).Path
$tag=(Get-Date -Format 'yyyyMMdd-HHmmss')+'-'+[guid]::NewGuid().ToString('N').Substring(0,8)
$resultDir=if($OutputDirectory){[IO.Path]::GetFullPath($OutputDirectory)}else{Join-Path $projectRoot "runs\hardening-acceptance-$tag"}
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
        Invoke-Logged rust-tests @('cargo','test','--release','--locked')
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
    Invoke-Logged elastic-cycle @($python,'tools\verify_hardening_elastic_cycle.py','--solver',$solver,'--rejected','runs\hardening-elastic-cycle-counterexample-20261024','--output-dir',"$resultDir\elastic-cycle")
    Invoke-Logged material @($python,'-u','tools\verify_hardening_material.py','--solver',$solver,'--output-dir',"$resultDir\material")
    Invoke-Logged triaxial @($python,'-u','tools\verify_hardening_triaxial.py','--solver',$solver,'--output-dir',"$resultDir\triaxial")
    foreach($order in @('tet4','tet10')) {
        $mesh=if($order -eq 'tet4'){'examples\mc-column-mesh.json'}else{'examples\hardening-column-tet10-mesh.json'}
        Invoke-Logged $order @($python,'tools\run_hardening_fixture.py','--solver',$solver,'--mesh',$mesh,'--ground','examples\ground-hardening-column.json','--stages','examples\stages-hardening-column.json','--output',"$resultDir\$order")
        Invoke-Logged "$order-history" @($python,'tools\verify_hardening_history.py','--solver',$solver,'--run',"$resultDir\$order",'--output-dir',"$resultDir\$order-history")
        Invoke-Logged "$order-checkpoint" @($python,'tools\verify_construction_checkpoint.py','--solver',$solver,'--run',"$resultDir\$order",'--output-dir',"$resultDir\$order-checkpoint")
    }
    Invoke-Logged column-reference @($python,'-u','tools\verify_hardening_column.py','--solver',$solver,'--run',"$resultDir\tet4",'--output-dir',"$resultDir\column-reference")
    Invoke-Logged local @($shellExecutable,'-NoProfile','-File',"$PSScriptRoot\Invoke-LocalExcavation.ps1",'-GroundFile','examples\ground-hardening-local.json','-StageFile','examples\stages-hardening-local.json','-MeshSizeM','4','-SolverPath',$solver,'-OutputDirectory',"$resultDir\local")
    Invoke-Logged local-history @($python,'tools\verify_hardening_history.py','--solver',$solver,'--run',"$resultDir\local",'--output-dir',"$resultDir\local-history")
    Invoke-Logged local-checkpoint @($python,'tools\verify_construction_checkpoint.py','--solver',$solver,'--run',"$resultDir\local",'--output-dir',"$resultDir\local-checkpoint")
    Invoke-Logged legacy @($python,'tools\verify_hardening_legacy.py','--solver',$solver,'--baseline','runs\local-acceptance-20261024-033221-8ada5d9c','--output-dir',"$resultDir\legacy")
    if((Get-FileHash -LiteralPath target\release\soil-fem-3d.exe -Algorithm SHA256).Hash.ToLower() -ne $digest){throw 'Build binary changed during suite'}
    @{schema='soil-fem-hardening-suite/1';passed=$true;solver_sha256=$digest} | ConvertTo-Json | Set-Content -LiteralPath "$resultDir\verification.json" -Encoding utf8
} finally {Stop-Transcript | Out-Null;Pop-Location}
& $python "$projectRoot\tools\record_hardening_validation.py" --run $resultDir
if($LASTEXITCODE -ne 0){throw "Hardening evidence binding failed: $resultDir"}
& $python "$projectRoot\tools\verify_acceptance_packet.py" $resultDir --output "$resultDir\integrity.json"
if($LASTEXITCODE -ne 0){throw "Acceptance packet integrity failed: $resultDir"}
