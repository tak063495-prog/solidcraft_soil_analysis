param(
    [string]$Python = '.venv/Scripts/python.exe',
    [string]$Output = 'runs/NEW-nc-fold-confirmation-PC-3',
    [string]$Pc2 = 'runs/NEW-nc-fold-confirmation-PC-2',
    [string]$ResumeFrom = ''
)
$ErrorActionPreference = 'Stop'
$Repo = (Get-Location).Path
if (-not (Test-Path 'tools/nc_fold_path.py')) { throw '対象worktreeのルートから実行してください' }
if (Test-Path $Output) { throw '新しい出力先を指定してください。既存出力を削除/上書きしません' }
$ResumeRoot = $null
if ($ResumeFrom) {
    $ResumeRoot = (Resolve-Path -LiteralPath $ResumeFrom).Path
    $NewRoot = $ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath($Output)
    if ($NewRoot -eq $ResumeRoot -or $NewRoot.StartsWith($ResumeRoot + [System.IO.Path]::DirectorySeparatorChar, [System.StringComparison]::OrdinalIgnoreCase)) {
        throw 'PC-3参照元の内部には出力できません。読み取り専用で保持します'
    }
}
$Py = (Resolve-Path $Python).Path
$Pc2Root = (Resolve-Path $Pc2).Path
New-Item -ItemType Directory -Path $Output -ErrorAction Stop | Out-Null
$TaskRun = (Resolve-Path $Output).Path
$env:OPENBLAS_NUM_THREADS='1'; $env:MKL_NUM_THREADS='1'; $env:OMP_NUM_THREADS='1'
$env:RUST_TEST_THREADS='1'; $env:PYTHONUTF8='1'; $env:PYTHONIOENCODING='utf-8'
[Console]::OutputEncoding = [System.Text.UTF8Encoding]::new($false)
$OutputEncoding = [Console]::OutputEncoding
$Head = (git rev-parse HEAD).Trim()
if ($LASTEXITCODE -ne 0) { throw 'HEAD取得失敗' }
$Head | Out-File "$TaskRun/head.txt" -Encoding utf8
Get-Date -Format o | Out-File "$TaskRun/start.txt" -Encoding utf8
$PSVersionTable | ConvertTo-Json -Depth 8 | Out-File "$TaskRun/powershell-version.json" -Encoding utf8
$CurrentStep = 'source-manifest'
$GateDiagnostic = $null
$Commands = @()
function Invoke-Step([string]$Name, [string]$Exe, [string[]]$Arguments) {
    $script:CurrentStep = $Name
    $Start = Get-Date
    $Record = [ordered]@{ name=$Name; executable=$Exe; arguments=$Arguments; head=$Head; start=$Start.ToString('o'); return_code=$null }
    $Record | ConvertTo-Json -Depth 8 | Out-File "$TaskRun/$Name-command.json" -Encoding utf8
    Get-Command $Exe -CommandType Application -ErrorAction Stop | Out-Null
    # Windows PowerShell 5.1 represents native stderr as ErrorRecord. Retain it
    # in the log and judge the native exit code, including on a failed process.
    $Previous = $ErrorActionPreference
    $ErrorActionPreference = 'Continue'
    try {
        & $Exe @Arguments 2>&1 | Tee-Object -FilePath "$TaskRun/$Name.log"
        $Rc = $LASTEXITCODE
    } finally { $ErrorActionPreference = $Previous }
    $Record.return_code=$Rc; $Record.end=(Get-Date).ToString('o')
    $Record.elapsed_seconds=((Get-Date)-$Start).TotalSeconds
    $Record | ConvertTo-Json -Depth 8 | Out-File "$TaskRun/$Name-command.json" -Encoding utf8
    if ($Rc -ne 0) { throw "$Name rc=$Rc。report/log/途中出力を保持して停止" }
}
function Read-Json([string]$Path) {
    Get-Content -LiteralPath $Path -Raw -Encoding utf8 | ConvertFrom-Json
}
function Write-NewJson([string]$Path, $Value) {
    $Bytes = [System.Text.UTF8Encoding]::new($false).GetBytes(($Value | ConvertTo-Json -Depth 16))
    $Stream = [System.IO.File]::Open($Path, [System.IO.FileMode]::CreateNew, [System.IO.FileAccess]::Write)
    try { $Stream.Write($Bytes, 0, $Bytes.Length) } finally { $Stream.Dispose() }
}
function Assert-PathCount([string]$Path, [int]$Count) {
    $Report = Read-Json "$Path/report.json"
    if ($Report.equilibrium_samples -ne $Count -or $Count -le 0) { throw "$Path の samples不足。reportを保持" }
    if ($Report.stop_reason -ne 'bounded step count; inconclusive') { throw "$Path の別の停止理由: $($Report.stop_reason)" }
    if ($Report.passed -ne $false) { throw '診断を本番受入にしてはいけません' }
}
$NativeReports = @()
try {
    foreach ($Denom in @(128,256)) {
        foreach ($Size in @('coarse','fine','finer')) {
            $Path = Join-Path $Pc2Root "native$Denom-$Size-audit/fold-results.json"
            if (-not (Test-Path $Path)) { throw "PC-2の未変更native監査が必要: $Path" }
            $NativeReports += $Path
        }
    }
    $Manifest = @($NativeReports | ForEach-Object { [ordered]@{path=$_; sha256=(Get-FileHash -Algorithm SHA256 -LiteralPath $_).Hash.ToLowerInvariant(); read_only_source='PC-2; no recalculation or overwrite'} })
    $Manifest | ConvertTo-Json -Depth 8 | Out-File "$TaskRun/pc2-native-source.json" -Encoding utf8
    Invoke-Step 'runtime' $Py @('-c','import sys,numpy,scipy,platform; print(sys.version); print(numpy.__version__,scipy.__version__); print(platform.platform())')
    Invoke-Step 'build' 'cargo' @('build','--locked','--release')
    Invoke-Step 'focused-tests' $Py @('-m','unittest','discover','-s','tests','-p','test_nc_independent_arc_residual.py','-v')
    Invoke-Step 'gate-tests' $Py @('-m','unittest','discover','-s','tests','-p','test_nc_pc_gate.py','-v')
    if ($ResumeRoot) {
        $GateDiagnostic = "$TaskRun/pc3-reuse-validation.json"
        Invoke-Step 'pc3-reuse-validation' $Py @('tools/nc_pc_gate.py','--resume-from',$ResumeRoot,'--repo',$Repo,'--output',$GateDiagnostic)
        $GateDiagnostic = $null
    }
    $Reports = @($NativeReports)
    foreach ($Steps in @(400,800)) {
        $Ref = "$TaskRun/reference$Steps"; $Arc = "$TaskRun/rate-arc$Steps"
        if ($ResumeRoot -and $Steps -eq 400) {
            $Ref = "$ResumeRoot/reference400"; $Arc = "$ResumeRoot/rate-arc400"
        } else {
            Invoke-Step "reference$Steps" $Py @('tools/low_preload_nc_reference.py','--output',$Ref,'--steps',"$Steps",'--method','rk4','--budget-seconds','7200')
        }
        $Reference = Read-Json "$Ref/reference.json"
        # The old reference stiffness gate normally stops the input-only history.
        # Its last accepted excavate_upper state, not rc=0, is the required seed.
        if ($Reference.last_accepted.name -ne 'excavate_upper') { throw "reference$Steps に必要な撤去起点なし" }
        if (-not ($ResumeRoot -and $Steps -eq 400)) {
            Invoke-Step "rate-arc$Steps" $Py @('interim/2026-11-09-nc-pseudo-arclength/independent_arc.py','--reference',$Ref,'--output',$Arc,'--step-size','0.00002','--count','128')
        }
        $GateDiagnostic = "$TaskRun/rate-arc$Steps-gate.json"
        Invoke-Step "rate-arc$Steps-gate" $Py @('tools/nc_pc_gate.py','--arc',$Arc,'--output',$GateDiagnostic)
        $GateDiagnostic = $null
        foreach ($Spec in @(@('coarse','0.00002',64,4),@('fine','0.00001',128,4),@('material8','0.00002',64,8))) {
            $Name = "independent$Steps-$($Spec[0])"; $PathOut = "$TaskRun/$Name"
            Invoke-Step $Name $Py @('tools/nc_fold_path.py','--map','independent','--reference',$Ref,'--arc',$Arc,'--step-size',"$($Spec[1])",'--count',"$($Spec[2])",'--material-steps',"$($Spec[3])",'--budget-seconds','1800','--output',$PathOut)
            Assert-PathCount $PathOut $Spec[2]
            Invoke-Step "$Name-audit" $Py @('tools/analyze_nc_fold_path.py','--path',$PathOut,'--output',"$PathOut-audit")
            $Reports += "$PathOut-audit/fold-results.json"
        }
        # §4.1: same new source/runtime/seed/ds, 32+32 vs the direct 64 above.
        $First = "$TaskRun/independent$Steps-restart-first32"
        $Second = "$TaskRun/independent$Steps-restart-second32"
        Invoke-Step "independent$Steps-first32" $Py @('tools/nc_fold_path.py','--map','independent','--reference',$Ref,'--arc',$Arc,'--step-size','0.00002','--count','32','--material-steps','4','--budget-seconds','1800','--output',$First)
        Assert-PathCount $First 32
        Invoke-Step "independent$Steps-second32" $Py @('tools/nc_fold_path.py','--map','independent','--reference',$Ref,'--arc',$Arc,'--step-size','0.00002','--count','32','--material-steps','4','--budget-seconds','1800','--resume',"$First/restart.json",'--output',$Second)
        Assert-PathCount $Second 32
        foreach ($PathOut in @($First,$Second)) {
            $Name = Split-Path $PathOut -Leaf
            Invoke-Step "$Name-audit" $Py @('tools/analyze_nc_fold_path.py','--path',$PathOut,'--output',"$PathOut-audit")
        }
        $Check = @'
import json,sys
from pathlib import Path
root,steps=Path(sys.argv[1]),sys.argv[2]
read=lambda p:[json.loads(x) for x in p.read_text(encoding='utf-8').splitlines()]
base=f'independent{steps}'
names=('samples.jsonl','gp-tangents.jsonl','linear-systems.jsonl')
pairs=[(n,read(root/f'{base}-coarse'/n),read(root/f'{base}-restart-first32'/n)+read(root/f'{base}-restart-second32'/n)) for n in names]
assert len(pairs[0][1])==len(pairs[0][2])==64
assert all(a==b for n,a,b in pairs),{n:a==b for n,a,b in pairs}
assert json.loads((root/f'{base}-coarse/restart.json').read_text())==json.loads((root/f'{base}-restart-second32/restart.json').read_text())
print(json.dumps(dict(all_saved_fields_gp_tangents_actual_systems_and_terminal_checkpoint_exact=True,rows={n:len(a) for n,a,b in pairs}),indent=2))
'@
        Invoke-Step "independent$Steps-restart-comparison" $Py @('-c',$Check,$TaskRun,"$Steps")
    }
    # §5: retain PC-2 native provenance; these are not PC-3 native calculations.
    $CompareArgs = @('tools/compare_nc_fold_reports.py','--reports') + $Reports + @('--output',"$TaskRun/comparison.json")
    Invoke-Step 'comparison' $Py $CompareArgs
    $Comparison = Read-Json "$TaskRun/comparison.json"
    if ($Comparison.passed -ne $false -or $Comparison.nc01_complete -ne $false) { throw '比較の未達フラグを維持してください' }
    $ManifestAfter = @($NativeReports | ForEach-Object { (Get-FileHash -Algorithm SHA256 -LiteralPath $_).Hash.ToLowerInvariant() })
    for ($I=0; $I -lt $Manifest.Count; $I++) { if ($Manifest[$I].sha256 -ne $ManifestAfter[$I]) { throw 'PC-2元監査が変化' } }
    if ($ResumeRoot) {
        $GateDiagnostic = "$TaskRun/pc3-reuse-after.json"
        Invoke-Step 'pc3-reuse-after' $Py @('tools/nc_pc_gate.py','--verify-reuse',"$TaskRun/pc3-reuse-validation.json",'--output',$GateDiagnostic)
        $GateDiagnostic = $null
    }
    'Bounded procedures finished; passed=false, NC-01 incomplete. Review every verdict and field comparison.' | Out-File "$TaskRun/procedure-ended.txt" -Encoding utf8
} catch {
    $Failure = [ordered]@{passed=$false; nc01_complete=$false; p6_complete=$false; stop_reason=$_.ToString(); failure_stage=$CurrentStep; head=$Head; time=(Get-Date).ToString('o')}
    if ($GateDiagnostic -and (Test-Path -LiteralPath $GateDiagnostic)) {
        $Failure.gate_diagnostic_path = $GateDiagnostic
        $Failure.gate_diagnostic = Read-Json $GateDiagnostic
        $Failure.failed_conditions = $Failure.gate_diagnostic.failures
    }
    $FailureName = if ($ResumeRoot) { 'pc4-failure.json' } else { 'pc3-failure.json' }
    Write-NewJson "$TaskRun/$FailureName" $Failure
    throw
} finally {
    Get-Date -Format o | Out-File "$TaskRun/end.txt" -Encoding utf8
}
