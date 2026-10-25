param([string]$OutputDirectory)
$ErrorActionPreference='Stop'
$projectRoot=Split-Path -Parent $PSScriptRoot
$python=Join-Path $projectRoot '.venv\Scripts\python.exe'
$shellExecutable=(Get-Process -Id $PID).Path
$tag=(Get-Date -Format 'yyyyMMdd-HHmmss')+'-'+[guid]::NewGuid().ToString('N').Substring(0,8)
$resultDir=if($OutputDirectory){[IO.Path]::GetFullPath($OutputDirectory)}else{Join-Path $projectRoot "runs\construction-acceptance-$tag"}
if(Test-Path -LiteralPath $resultDir){throw "Output exists: $resultDir"}
New-Item -ItemType Directory -Path $resultDir | Out-Null
function Invoke-Logged([string]$Name,[string[]]$Arguments) {
    & $python "$projectRoot\tools\run_logged.py" "$resultDir\$Name.log" @Arguments
    if($LASTEXITCODE -ne 0){throw "Verification failed: $Name (retained $resultDir)"}
}
Push-Location -LiteralPath $projectRoot
try {
    Invoke-Logged baseline @($shellExecutable,'-NoProfile','-File',"$PSScriptRoot\Test-OneWay.ps1",'-OutputDirectory',"$resultDir\baseline")
    foreach($name in @('source-snapshot.json','source-snapshot.zip')){Copy-Item -LiteralPath "$resultDir\baseline\$name" -Destination "$resultDir\$name"}
    New-Item -ItemType Directory -Path "$resultDir\solver" | Out-Null
    $solver="$resultDir\solver\soil-fem-3d.exe"
    Copy-Item -LiteralPath "$resultDir\baseline\solver\soil-fem-3d.exe" -Destination $solver
    Invoke-Logged projection-reference @($python,'tools\test_surface_projection_reference.py')
    Invoke-Logged cases @($python,'tools\verify_one_way_construction.py','--solver',$solver,'--output-dir',"$resultDir\cases",'--cad',"$resultDir\baseline\baseline\core\cad")
    & $python tools\record_one_way_construction_validation.py --run $resultDir
    if($LASTEXITCODE -ne 0){throw "Construction evidence binding failed: $resultDir"}
    & $python tools\verify_acceptance_packet.py $resultDir --output "$resultDir\integrity.json"
    if($LASTEXITCODE -ne 0){throw "Construction packet integrity failed: $resultDir"}
} finally {Pop-Location}
