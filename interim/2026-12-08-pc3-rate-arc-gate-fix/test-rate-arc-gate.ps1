param(
    [string]$Python = 'python',
    [Parameter(Mandatory=$true)][string]$Output
)
$ErrorActionPreference = 'Stop'
New-Item -ItemType Directory -Force -Path $Output | Out-Null
$Output = (Resolve-Path $Output).Path
function Read-Json([string]$Path) {
    Get-Content -LiteralPath $Path -Raw -Encoding utf8 | ConvertFrom-Json
}
$Samples = @(0..127 | ForEach-Object { [ordered]@{step=$_} })
$Cases = @('pc3-values','count127','residual','stop-string','empty','missing','residual-string','passed-string','nested','count-string')
$Results = @()
foreach ($Name in $Cases) {
    $Arc = Join-Path $Output $Name
    New-Item -ItemType Directory -Path $Arc -ErrorAction Stop | Out-Null
    $Report = [ordered]@{equilibrium_samples=128; stop_reason=$null; maximum_bordered_true_linear_residual=5.4291880423103136e-11; passed=$false}
    $Rows = $Samples
    switch ($Name) {
        'count127' { $Rows = $Samples[0..126] }
        'residual' { $Report.maximum_bordered_true_linear_residual=1.1e-10 }
        'stop-string' { $Report.stop_reason='null' }
        'empty' { $Rows=@() }
        'missing' { $Report.Remove('stop_reason') }
        'residual-string' { $Report.maximum_bordered_true_linear_residual='5.4291880423103136e-11' }
        'passed-string' { $Report.passed='false' }
        'nested' { $Rows=, $Samples }
        'count-string' { $Report.equilibrium_samples='128' }
    }
    $Report | ConvertTo-Json -Depth 8 | Out-File "$Arc/report.json" -Encoding utf8
    ConvertTo-Json -InputObject $Rows -Depth 8 | Out-File "$Arc/samples.json" -Encoding utf8
    & $Python tools/nc_pc_gate.py --arc $Arc --output "$Arc/gate.json" > "$Arc/gate.log"
    $Rc = $LASTEXITCODE
    $Result = Read-Json "$Arc/gate.json"
    $Expected = $Name -eq 'pc3-values'
    if (($Expected -and $Rc -ne 0) -or (-not $Expected -and $Rc -ne 1)) { throw "$Name unexpected rc=$Rc" }
    if ($Result.accepted -ne $Expected -or $Result.passed -ne $false -or $Result.nc01_complete -ne $false -or $Result.p6_complete -ne $false) { throw "$Name unexpected flags" }
    if (-not $Expected -and @($Result.failures).Count -lt 1) { throw "$Name missing diagnostic" }
    $Results += [ordered]@{name=$Name; rc=$Rc; accepted=$Result.accepted; failures=$Result.failures}
}
# Execute the exact old reader/count expression on this real shell.
$OldCount = @(Read-Json "$Output/pc3-values/samples.json").Count
$Direct = Read-Json "$Output/pc3-values/samples.json"
if ($Direct.Count -ne 128) { throw 'fixture is not a top-level 128-item array' }
if ($PSVersionTable.PSVersion.Major -eq 5 -and $OldCount -ne 1) { throw "unexpected Windows PowerShell 5.1 old count=$OldCount" }
if ($PSVersionTable.PSVersion.Major -ge 7 -and $OldCount -ne 128) { throw "unexpected PowerShell 7 old count=$OldCount" }
[ordered]@{ps_version=$PSVersionTable.PSVersion.ToString(); old_wrapped_count=$OldCount; direct_count=$Direct.Count; passed=$Cases.Count; failed=0; skipped=0; cases=$Results} |
    ConvertTo-Json -Depth 10 | Out-File "$Output/shell-result.json" -Encoding utf8
Get-Content "$Output/shell-result.json" -Raw -Encoding utf8
