param([switch]$SkipChecks,[string]$ReuseCadDirectory)
$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
$env:PATH = "$env:USERPROFILE\.cargo\bin;$env:PATH"
$python = Join-Path $projectRoot '.venv\Scripts\python.exe'
$shellExecutable = (Get-Process -Id $PID).Path
$tag = (Get-Date -Format 'yyyyMMdd-HHmmss') + '-' + [guid]::NewGuid().ToString('N').Substring(0,8)
$resultDir = Join-Path $projectRoot "runs\tet10-boundary-test-$tag"
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
    $solverHash = (Get-FileHash -LiteralPath 'target\release\soil-fem-3d.exe' -Algorithm SHA256).Hash.ToLower()
    $cadRoot = $resultDir
    if ($ReuseCadDirectory) {
        $cadRoot = (Resolve-Path -LiteralPath $ReuseCadDirectory).Path
        foreach ($name in @('flat','curved','curved-fine')) {
            $accepted = Get-Content -LiteralPath "$cadRoot\$name\accepted-stage-run.json" -Raw | ConvertFrom-Json
            if ($accepted.solver_binary_sha256 -ne $solverHash) { throw "Reused CAD result binary mismatch: $name" }
            foreach ($entry in $accepted.inputs_and_results_sha256.PSObject.Properties) {
                if ((Get-FileHash -LiteralPath $entry.Name -Algorithm SHA256).Hash.ToLower() -ne $entry.Value) { throw "Reused CAD input/result mismatch: $($entry.Name)" }
            }
        }
    } else {
        Invoke-Logged cad-tests @($python,'-m','unittest','discover','-s','tests','-p','test_cad_excavation.py','-v')
        Invoke-Logged flat @($shellExecutable,'-NoProfile','-File',"$PSScriptRoot\Invoke-LocalExcavation.ps1",'-MeshSizeM','4','-OutputDirectory',"$cadRoot\flat")
        Invoke-Logged curved @($shellExecutable,'-NoProfile','-File',"$PSScriptRoot\Invoke-LocalExcavation.ps1",'-Design','examples\local-excavation-cylinder.json','-MeshSizeM','4','-OutputDirectory',"$cadRoot\curved")
        Invoke-Logged curved-fine @($shellExecutable,'-NoProfile','-File',"$PSScriptRoot\Invoke-LocalExcavation.ps1",'-Design','examples\local-excavation-cylinder.json','-MeshSizeM','3','-OutputDirectory',"$cadRoot\curved-fine")
    }
    Invoke-Logged high-order @($python,'tools\verify_tet10_package.py','--curved-run',"$cadRoot\curved",'--fine-curved-run',"$cadRoot\curved-fine",'--output-dir',"$resultDir\high-order")
    Invoke-Logged spatial @($python,'tools\compare_quadratic_meshes.py',"$cadRoot\curved","$cadRoot\curved-fine","$resultDir\spatial-sensitivity.json")
    Invoke-Logged local-reference @($python,'tools\verify_quadratic_local.py','--ground-run',"$cadRoot\flat\ground",'--output-dir',"$resultDir\local-reference")
    Invoke-Logged prescribed @($shellExecutable,'-NoProfile','-File',"$PSScriptRoot\Invoke-Stages.ps1",'-GroundRun',"$cadRoot\flat\ground",'-StageFile','examples\stages-prescribed.json','-OutputDirectory',"$resultDir\prescribed")
    Invoke-Logged checkpoint @($python,'tools\verify_checkpoint.py','--ground-run',"$cadRoot\flat\ground",'--output-dir',"$resultDir\checkpoint")
    Invoke-Logged independent-plastic @($python,'tools\verify_stages.py','--solver','target\release\soil-fem-3d.exe','--output-dir',"$resultDir\independent-plastic",'--small-strain')
    if ((Get-FileHash -LiteralPath 'target\release\soil-fem-3d.exe' -Algorithm SHA256).Hash.ToLower() -ne $solverHash) { throw 'Solver changed during suite' }
    $record = @{
        schema = 'soil-fem-tet10-boundary-suite/1'
        scope = 'Tet10 CAD local excavation, general boundary conditions, independent calculations and exact history restart; P6 overall incomplete'
        solver_sha256 = $solverHash
        cargo_lock_sha256 = (Get-FileHash -LiteralPath 'Cargo.lock' -Algorithm SHA256).Hash.ToLower()
        rust_format_lint_build_included = (-not $SkipChecks.IsPresent)
        result_directory = $resultDir
        cad_directory = $cadRoot
        passed = $true
    }
    $record | ConvertTo-Json | Set-Content -LiteralPath "$resultDir\verification.json" -Encoding utf8
    Write-Host "Tet10 and boundary verification passed: $resultDir"
} finally { Stop-Transcript | Out-Null; Pop-Location }
