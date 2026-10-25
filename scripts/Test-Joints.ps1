param([switch]$SkipChecks,[string]$OutputDirectory,[string]$MccBaseline)
$ErrorActionPreference='Stop'
$projectRoot=Split-Path -Parent $PSScriptRoot
$env:PATH="$env:USERPROFILE\.cargo\bin;$env:PATH"
$python=Join-Path $projectRoot '.venv\Scripts\python.exe'
$shellExecutable=(Get-Process -Id $PID).Path
if($MccBaseline){$MccBaseline=(Resolve-Path -LiteralPath $MccBaseline).Path}
$tag=(Get-Date -Format 'yyyyMMdd-HHmmss')+'-'+[guid]::NewGuid().ToString('N').Substring(0,8)
$resultDir=if($OutputDirectory){[IO.Path]::GetFullPath($OutputDirectory)}else{Join-Path $projectRoot "runs\joint-test-$tag"}
if(Test-Path -LiteralPath $resultDir){throw "Output exists: $resultDir"}
New-Item -ItemType Directory -Path $resultDir | Out-Null
Push-Location -LiteralPath $projectRoot
Start-Transcript -Path "$resultDir\verification.log" | Out-Null
function Invoke-Logged([string]$Name,[string[]]$Arguments) {
    & $python tools\run_logged.py "$resultDir\$Name.log" @Arguments
    if($LASTEXITCODE -ne 0){throw "Verification step failed: $Name (retained $resultDir)"}
}
try {
    if(-not $SkipChecks){
        Invoke-Logged format @('cargo','fmt','--check')
        Invoke-Logged rust-tests @('cargo','test','--release','--locked','--','--nocapture')
        Invoke-Logged lint @('cargo','clippy','--locked','--all-targets','--','-D','warnings')
        Invoke-Logged build @('cargo','build','--release','--locked')
    }
    New-Item -ItemType Directory -Path "$resultDir\solver" | Out-Null
    $solver="$resultDir\solver\soil-fem-3d.exe"
    Copy-Item -LiteralPath 'target\release\soil-fem-3d.exe' -Destination $solver
    $digest=(Get-FileHash -LiteralPath $solver -Algorithm SHA256).Hash.ToLower()
    Invoke-Logged cad-tests @($python,'-m','unittest','discover','-s','tests','-p','test_cad_excavation.py','-v')
    Invoke-Logged split-tests @($python,'-m','unittest','discover','-s','tests','-p','test_joint_split.py','-v')
    foreach($material in @('mc','mcc')) {
        Invoke-Logged $material @($python,'tools\run_joint_fixture.py','--solver',$solver,'--ground',"examples\ground-joint-$material.json",'--stages','examples\stages-joint-verification.json','--output',"$resultDir\$material")
        Invoke-Logged "$material-history" @($python,'tools\verify_joint_history.py','--solver',$solver,'--run',"$resultDir\$material",'--output-dir',"$resultDir\$material-history")
        Invoke-Logged "$material-checkpoint" @($python,'tools\verify_construction_checkpoint.py','--solver',$solver,'--run',"$resultDir\$material",'--output-dir',"$resultDir\$material-checkpoint")
    }
    Invoke-Logged independent @($python,'tools\verify_joint_coupling.py',"$resultDir\mc","$resultDir\independent.json")
    Invoke-Logged local-integration-refinement @($python,'tools\verify_joint_coupling.py',"$resultDir\mc","$resultDir\refinement.json",'--solver',$solver)
    Invoke-Logged excavation @($python,'tools\run_joint_fixture.py','--solver',$solver,'--ground','examples\ground-joint-mc.json','--stages','examples\stages-joint-excavation.json','--output',"$resultDir\excavation")
    Invoke-Logged excavation-checkpoint @($python,'tools\verify_construction_checkpoint.py','--solver',$solver,'--run',"$resultDir\excavation",'--output-dir',"$resultDir\excavation-checkpoint")
    Invoke-Logged normal-column @($python,'tools\run_joint_fixture.py','--solver',$solver,'--ground','examples\ground-joint-mc.json','--stages','examples\stages-joint-normal-column.json','--output',"$resultDir\normal-column")
    Invoke-Logged normal-column-theory @($python,'tools\verify_joint_column_theory.py',"$resultDir\normal-column","$resultDir\normal-column-theory.json")
    Invoke-Logged normal-column-checkpoint @($python,'tools\verify_construction_checkpoint.py','--solver',$solver,'--run',"$resultDir\normal-column",'--output-dir',"$resultDir\normal-column-checkpoint")
    Invoke-Logged explicit-normal-column @($python,'tools\run_joint_fixture.py','--solver',$solver,'--ground','examples\ground-joint-mc.json','--stages','examples\stages-joint-explicit-normal-column.json','--output',"$resultDir\explicit-normal-column")
    Invoke-Logged explicit-normal-column-theory @($python,'tools\verify_joint_column_theory.py',"$resultDir\explicit-normal-column","$resultDir\explicit-normal-column-theory.json")
    Invoke-Logged explicit-normal-column-checkpoint @($python,'tools\verify_construction_checkpoint.py','--solver',$solver,'--run',"$resultDir\explicit-normal-column",'--output-dir',"$resultDir\explicit-normal-column-checkpoint")
    Invoke-Logged unseating @($python,'tools\verify_joint_support.py','--solver',$solver,'--output-dir',"$resultDir\unseating")
    Invoke-Logged local-excavation @($python,'tools\run_joint_fixture.py','--solver',$solver,'--ground','examples\ground-joint-local-mc.json','--stages','examples\stages-joint-local-excavation.json','--mesh','examples\joint-column-local-mesh.json','--output',"$resultDir\local-excavation")
    Invoke-Logged local-excavation-checkpoint @($python,'tools\verify_construction_checkpoint.py','--solver',$solver,'--run',"$resultDir\local-excavation",'--output-dir',"$resultDir\local-excavation-checkpoint")
    Invoke-Logged cad @($shellExecutable,'-NoProfile','-File',"$PSScriptRoot\Invoke-Joints.ps1",'-SolverPath',$solver,'-OutputDirectory',"$resultDir\cad")
    Invoke-Logged cad-independent @($python,'tools\verify_joint_coupling.py',"$resultDir\cad\analysis","$resultDir\cad-independent.json")
    if($MccBaseline){
        Invoke-Logged legacy-cad @($shellExecutable,'-NoProfile','-File',"$PSScriptRoot\Invoke-LocalExcavation.ps1",'-MeshSizeM','4','-SolverPath',$solver,'-OutputDirectory',"$resultDir\legacy-cad")
        Invoke-Logged mcc-regression @($python,'tools\verify_mcc_regression.py',$MccBaseline,"$resultDir\legacy-cad\stage-results.json","$resultDir\mcc-regression.json")
    }
    if((Get-FileHash -LiteralPath 'target\release\soil-fem-3d.exe' -Algorithm SHA256).Hash.ToLower() -ne $digest){throw 'Build binary changed during suite'}
    @{schema='soil-fem-joint-suite/1';passed=$true;solver_sha256=$digest;rust_format_lint_build_included=(-not $SkipChecks.IsPresent)} | ConvertTo-Json | Set-Content -LiteralPath "$resultDir\verification.json" -Encoding utf8
    Write-Host "JOINT verification passed: $resultDir"
} finally {Stop-Transcript | Out-Null;Pop-Location}
if(-not $SkipChecks){
    & $python "$projectRoot\tools\record_joint_validation.py" --run $resultDir
    if($LASTEXITCODE -ne 0){throw "JOINT evidence binding failed: $resultDir"}
}
