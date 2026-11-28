param(
    [Parameter(Mandatory=$true)][string]$InputDirectory,
    [Parameter(Mandatory=$true)][string]$Trace,
    [ValidateSet('baseline','refined','repeated_refinement','extended_residual','ruiz','gmres_polish')]
    [string[]]$Modes = @('baseline','refined','repeated_refinement','extended_residual','ruiz','gmres_polish')
)
$ErrorActionPreference = 'Stop'
$directory = (Resolve-Path $InputDirectory).Path
$binary = (Resolve-Path $Trace).Path
$hash = (Get-FileHash -Algorithm SHA256 $binary).Hash.ToLowerInvariant()
$mesh = Get-Content (Join-Path $directory 'mesh.json') -Raw | ConvertFrom-Json
if ($mesh.tetrahedra10.Count -ne 96 -or $mesh.tetrahedra.Count -ne 96) { throw 'Expected 96 Tet10' }
foreach ($mode in $Modes) {
    $stage = if ($mode -eq 'baseline') { 'stages.json' } else { "stages-$mode.json" }
    $case = Get-Content (Join-Path $directory $stage) -Raw | ConvertFrom-Json
    if ($case.settings.relative_linear_tolerance -ne 1e-10) { throw 'Unexpected linear tolerance' }
    foreach ($name in @($mode,"$mode-process.json","$mode.stdout","$mode.stderr")) {
        if (Test-Path (Join-Path $directory $name)) { throw "Existing output $name" }
    }
    $clock = [Diagnostics.Stopwatch]::StartNew()
    $p = Start-Process -FilePath $binary -WorkingDirectory $directory -PassThru -NoNewWindow `
        -ArgumentList @('direct','mesh.json','ground.json',$stage,$mode) `
        -RedirectStandardOutput (Join-Path $directory "$mode.stdout") `
        -RedirectStandardError (Join-Path $directory "$mode.stderr")
    $peak = 0L
    $budget = $false
    while (-not $p.WaitForExit(250)) {
        $p.Refresh()
        if ($p.HasExited) { break }
        $peak = [Math]::Max($peak,$p.PeakWorkingSet64)
        if ($clock.Elapsed.TotalHours -ge 6 -or $peak -ge 32000000000 -or $p.PrivateMemorySize64 -ge 32000000000) {
            $budget = $true
            $p.Kill()
            break
        }
    }
    $p.WaitForExit()
    $record = [ordered]@{ mode=$mode; binary_sha256=$hash; seconds=$clock.Elapsed.TotalSeconds
        peak_hwm_bytes=$peak; budget_stop=$budget; exit_code=$p.ExitCode }
    $record | ConvertTo-Json | Out-File (Join-Path $directory "$mode-process.json") -Encoding utf8
    if ($budget -or $p.ExitCode -ne 0) { throw "$mode stopped. Inspect the process record and retained output." }
}
