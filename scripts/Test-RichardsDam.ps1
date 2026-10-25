param([string]$OutputDirectory)
$ErrorActionPreference='Stop'
$projectRoot=Split-Path -Parent $PSScriptRoot
$python=Join-Path $projectRoot '.venv\Scripts\python.exe'
$shellExecutable=(Get-Process -Id $PID).Path
$tag=(Get-Date -Format 'yyyyMMdd-HHmmss')+'-'+[guid]::NewGuid().ToString('N').Substring(0,8)
$resultDir=if($OutputDirectory){[IO.Path]::GetFullPath($OutputDirectory)}else{Join-Path $projectRoot "runs\richards-dam-acceptance-$tag"}
if(Test-Path -LiteralPath $resultDir){throw "Output exists: $resultDir"}
New-Item -ItemType Directory -Path $resultDir | Out-Null
function Invoke-Logged([string]$Name,[string[]]$Arguments) {
    & $python "$projectRoot\tools\run_logged.py" "$resultDir\$Name.log" @Arguments
    if($LASTEXITCODE -ne 0){throw "Verification failed: $Name (retained $resultDir)"}
}
Push-Location -LiteralPath $projectRoot
try {
    # The hydraulic core freezes *all* sources before its fresh compilation.
    Invoke-Logged core @($shellExecutable,'-NoProfile','-File',"$PSScriptRoot\Test-Richards.ps1",'-OutputDirectory',"$resultDir\core")
    foreach($name in @('source-snapshot.json','source-snapshot.zip')){Copy-Item -LiteralPath "$resultDir\core\$name" -Destination "$resultDir\$name"}
    New-Item -ItemType Directory -Path "$resultDir\solver" | Out-Null
    $solver="$resultDir\solver\soil-fem-3d.exe"
    Copy-Item -LiteralPath "$resultDir\core\solver\soil-fem-3d.exe" -Destination $solver
    Invoke-Logged reference-tests @($python,'-m','unittest','discover','-s','tests','-p','test_richards_dam_reference.py','-v')
    Invoke-Logged dam @($python,'tools\verify_richards_dam.py','--solver',$solver,'--output-dir',"$resultDir\dam")
    Invoke-Logged counterexample @($python,'tools\verify_richards_dam_counterexample.py','--solver',$solver,'--output-dir',"$resultDir\counterexample")
    & $python tools\record_richards_dam_validation.py --run $resultDir
    if($LASTEXITCODE -ne 0){throw "Dam evidence binding failed: $resultDir"}
    & $python tools\verify_acceptance_packet.py $resultDir --output "$resultDir\integrity.json"
    if($LASTEXITCODE -ne 0){throw "Dam packet integrity failed: $resultDir"}
} finally {Pop-Location}
