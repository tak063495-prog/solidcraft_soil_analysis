param([string]$OutputDirectory)
$ErrorActionPreference='Stop'
$projectRoot=Split-Path -Parent $PSScriptRoot
$env:PATH="$env:USERPROFILE\.cargo\bin;$env:PATH"
$python=Join-Path $projectRoot '.venv\Scripts\python.exe'
$shellExecutable=(Get-Process -Id $PID).Path
$tag=(Get-Date -Format 'yyyyMMdd-HHmmss')+'-'+[guid]::NewGuid().ToString('N').Substring(0,8)
$resultDir=if($OutputDirectory){[IO.Path]::GetFullPath($OutputDirectory)}else{Join-Path $projectRoot "runs\water-acceptance-$tag"}
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
    Invoke-Logged rust-tests @('cargo','test','--release','--locked')
    Invoke-Logged lint @('cargo','clippy','--locked','--all-targets','--','-D','warnings')
    Invoke-Logged build @('cargo','build','--release','--locked')
    New-Item -ItemType Directory -Path "$resultDir\solver" | Out-Null
    $solver="$resultDir\solver\soil-fem-3d.exe"
    Copy-Item -LiteralPath target\release\soil-fem-3d.exe -Destination $solver
    $digest=(Get-FileHash -LiteralPath $solver -Algorithm SHA256).Hash.ToLower()
    Invoke-Logged cad-tests @($python,'-m','unittest','discover','-s','tests','-p','test_cad_excavation.py','-v')
    Invoke-Logged split-tests @($python,'-m','unittest','discover','-s','tests','-p','test_joint_split.py','-v')
    Invoke-Logged wet-faces @($python,'tools\verify_water_faces.py','--solver',$solver,'--output-dir',"$resultDir\wet-faces")
    Invoke-Logged column-reference @($python,'-u','tools\verify_water_column.py','--solver',$solver,'--output-dir',"$resultDir\column-reference")
    foreach($profile in @('mc-tet4','mc-tet10','mcc-tet10')) {
        $material=if($profile.StartsWith('mcc')){'mcc'}else{'mc'}
        $mesh=if($profile -eq 'mc-tet4'){'examples\mc-column-mesh.json'}else{'examples\hardening-column-tet10-mesh.json'}
        Invoke-Logged $profile @($python,'tools\run_water_fixture.py','--solver',$solver,'--mesh',$mesh,'--ground',"examples\ground-water-column-$material.json",'--stages','examples\stages-water-column.json','--output',"$resultDir\$profile")
        $history=@($python,'tools\verify_water_history.py','--solver',$solver,'--run',"$resultDir\$profile",'--output-dir',"$resultDir\$profile-history")
        if($material -eq 'mcc'){$history+='--require-plastic'}
        Invoke-Logged "$profile-history" $history
        Invoke-Logged "$profile-checkpoint" @($python,'tools\verify_construction_checkpoint.py','--solver',$solver,'--run',"$resultDir\$profile",'--output-dir',"$resultDir\$profile-checkpoint")
    }
    Invoke-Logged cad-mc @($shellExecutable,'-NoProfile','-File',"$PSScriptRoot\Invoke-LocalExcavation.ps1",'-GroundFile','examples\ground-water-local-mc.json','-StageFile','examples\stages-water-local.json','-MeshSizeM','4','-SolverPath',$solver,'-OutputDirectory',"$resultDir\cad-mc")
    Invoke-Logged independent-3d @($python,'-u','tools\verify_water_linear_3d.py','--solver',$solver,'--run',"$resultDir\cad-mc",'--output-dir',"$resultDir\independent-3d")
    Invoke-Logged cad-mcc @($shellExecutable,'-NoProfile','-File',"$PSScriptRoot\Invoke-LocalExcavation.ps1",'-GroundFile','examples\ground-water-local.json','-StageFile','examples\stages-water-local-mcc.json','-MeshSizeM','4','-SolverPath',$solver,'-OutputDirectory',"$resultDir\cad-mcc")
    Invoke-Logged cad-mcc-history @($python,'tools\verify_water_history.py','--solver',$solver,'--run',"$resultDir\cad-mcc",'--output-dir',"$resultDir\cad-mcc-history",'--require-plastic')
    Invoke-Logged cad-mcc-checkpoint @($python,'tools\verify_construction_checkpoint.py','--solver',$solver,'--run',"$resultDir\cad-mcc",'--output-dir',"$resultDir\cad-mcc-checkpoint")
    Invoke-Logged legacy @($python,'tools\verify_water_legacy.py','--solver',$solver,'--output-dir',"$resultDir\legacy")
    if((Get-FileHash -LiteralPath target\release\soil-fem-3d.exe -Algorithm SHA256).Hash.ToLower() -ne $digest){throw 'Build binary changed during suite'}
    @{schema='soil-fem-water-suite/1';passed=$true;solver_sha256=$digest} | ConvertTo-Json | Set-Content -LiteralPath "$resultDir\verification.json" -Encoding utf8
} finally {Stop-Transcript | Out-Null;Pop-Location}
& $python "$projectRoot\tools\record_water_validation.py" --run $resultDir
if($LASTEXITCODE -ne 0){throw "Water evidence binding failed: $resultDir"}
& $python "$projectRoot\tools\verify_acceptance_packet.py" $resultDir --output "$resultDir\integrity.json"
if($LASTEXITCODE -ne 0){throw "Acceptance packet integrity failed: $resultDir"}
