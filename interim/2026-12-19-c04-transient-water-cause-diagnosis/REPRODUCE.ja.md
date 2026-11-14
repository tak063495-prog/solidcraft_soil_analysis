# C-04原因診断の再現（Windows / Linux）

C-04は未達、passed=false。PCは未実行。GPD WIN Max 2・使用可能約48GBを前提に、以下は12要素だけの手順とする。元の力/線形/参照残差・全場定義・カテゴリを緩めない。出力先は毎回新しい名前にし、失敗後のフォルダーを消して使い回さない。

## 準備と旧封印

Rust 1.99.1、Python 3.13、固定依存。既存PR #13の同じブランチを使い、マージしない。PowerShell例：

```powershell
git switch codex/plastic-spatial-convergence
git pull --ff-only
# この局所再現用の疎取得。全件試験用の全アーカイブ配置とは別。
git sparse-checkout set --no-cone '/*' '!/**/*.zip*' '!/**/*.tar*' '!/**/*.gz' '!/**/*.npz' '/interim/2026-11-22-wet-joint-transient-water/evidence.zip.part*' '/interim/2026-12-19-c04-transient-water-cause-diagnosis/evidence.zip.part*'
py -3.13 -m venv .venv-c04-cause
.\.venv-c04-cause\Scripts\Activate.ps1
python -m pip install -r requirements-reference.lock.txt
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
$env:CARGO_BUILD_JOBS='2'
$env:RUST_TEST_THREADS='2'
$Study='interim/2026-12-19-c04-transient-water-cause-diagnosis'
$Old='interim/2026-11-22-wet-joint-transient-water'
$Saved='work/NEW-C04-cause-old'
python "$Old/archive.py" --mode replay --archive $Old --output $Saved
if ($LASTEXITCODE -ne 0) { throw 'old SHA or exact isolated measurement failed' }
$Base="$Saved/unpacked/data/final"
```

旧封印の5part・1,462エントリー・ZIPのSHAと隔離再計測を検査する。これは新規求解ではない。RAM 2GB程度、出力約0.6GB、PC初期予算2〜15分（未実測）。Windowsのlongdoubleによる旧RT0参照の丸め差があれば失敗を保持し、同版WSL/Linuxで厳密再現する。許容値を追加して成功へ変更しない。

## クラウドと同じ小型全体平衡対照

```powershell
cargo build --locked --release
if ($LASTEXITCODE -ne 0) { throw 'CLI build failed' }
cargo build --locked --release --example wet_joint_transient_probe
if ($LASTEXITCODE -ne 0) { throw 'observer build failed' }
$Run='work/NEW-C04-cause-global'
python "$Study/c04cause_driver.py" --output $Run --initialize --probe target/release/examples/wet_joint_transient_probe.exe --solver target/release/soil-fem-3d.exe
python "$Study/c04cause_binding.py" --controls $Run --capture
python "$Study/c04cause_batch.py" --output $Run --timeout 1800
if ($LASTEXITCODE -ne 0) { throw 'global matrix incomplete; retain partial outputs' }
foreach ($Kind in @('native','adaptive')) {
  python "$Study/c04cause_driver.py" --output $Run --label "$Kind-m128-w32" --mechanical-steps 128 --water-steps 32 --timeout 1800
  if ($LASTEXITCODE -ne 0) { throw '128 mechanical control incomplete' }
}
python "$Old/study.py" --mode refusals --output $Run --folder scope-refusals --timeout 900
if ($LASTEXITCODE -ne 0) { throw 'current Richards JOINT refusal changed' }
python "$Study/c04cause_binding.py" --baseline $Base --controls $Run --report "$Run/old-baseline-binding.json"
if ($LASTEXITCODE -ne 0) { throw 'old baseline or captured source differs' }
python "$Study/c04cause_observe.py" --baseline $Base --controls $Run --output work/NEW-C04-cause-observation
if ($LASTEXITCODE -ne 0) { throw 'full-field/global audit incomplete' }
python "$Study/c04cause_reference_assessment.py" --baseline $Base --controls $Run --report work/NEW-C04-cause-reference.json
```

旧Linux封印と新規Windowsバイナリの場のビット一致は未保証。bindingで差が出た場合は失敗として保存し、旧値を直さず、同版WSL/Linuxで厳密照合する。Windows差を新しい許容帯で受け入れない。最初の10経路＋追加2経路、全7,008受理状態。切戻しなしの入力列と元の平衡基準を検査する。クラウド実測は各プロセスJSONを参照。小型求解のRSSは約50〜55MiB（Python親からのLinuxサンプリングで子プロセス全体の厳密最大ではない）。Windowsは**1経路1〜15分、行列30〜180分、求解RAM 2GB・後処理RAM 4GB**を保守的な初期予算とする。未実測であり速度・完走保証ではない。旧封印の再監査・ビルド・全件試験は別時間。各子プロセスの時間切れ・未収束を未達として残す。

`c04cause_driver.py` の単独独立経路には外側の時間制限がないため、長い実験は下記PowerShellの予算付き実行を使う。batchは各子を監視し、タイムアウト時もログ・部分場・失敗JSONを保持する。Linuxでは `.exe` を外し、同じコマンドを使う。

## 出力と見る値

|出力|照合|
|---|---|
|`old-baseline-binding.json`|本体/独立の元の入力・416状態全場が旧封印と完全一致。Linuxから別Windows版のbinary SHAを付け替えない。|
|各`fields.jsonl`、本体`raw/accepted.jsonl`・`constitutive.jsonl`|全受理時刻・GPの全場、実分割、力残差。局所再生と区別する。|
|`mismatch-w32.json`|133件すべて保持。試行/戻り後余裕、最終塑性増分、元のカテゴリ。境界帯のラベルで一致へ置換しない。|
|`dissipation-identity-w*.json`|端点traction項・塑性増分項・増分内部残差。元の散逸差をこの3ベクトルで説明し、ノルムを足して原因割合にしない。|
|`local-replay-w*-c*.json`|共通接触分割での保存jump再生。`global_equilibrium_solved=false`を保持。|
|`measurements.json`|元の全履歴/13端点L2、絶対差、GP対応、全経路の別再組立平衡、工程別の散逸増分差。|
|`reference-assessment.json`|全3水準を同じ粗い座標で比較した差・次数・cos。水位方向だけの減少を全体連続解の保証にしない。|
|`scope-refusals/verification.json`|3現行拒否。Richards水理単独の旧RT0/水収支と湿潤JOINT連成を分ける。|
|各`process.json`・ログ・`failure*.json`|未完了・時間切れ・不収束を削除せず区別する。|

本体力1e-8・真の線形1e-10、独立力1e-10を維持。旧RT0基準も変更しない。未定義の全場受入条件に新しい値を入れない。`--boost-nontrivial` で散逸差が縮小しても元のC-04を合格にしない。

## PC：上限内の長系列と12 Tet10（今回は未実行）

追加の目的は全体経路・高次要素で同じ寄与が残るかの確認であり、元の合格条件を変えることではない。768/6144要素や実大規模連成はこの手順に含めない。

予算付きPowerShell関数（1経路最大4時間、メモリ監視の目安8GB。強制終了も未完了としてログ/部分出力を残す）：

```powershell
function Invoke-C04Bounded([string[]]$CmdArgs, [string]$LogBase, [int]$Seconds=14400) {
  $Py=(Get-Command python).Source
  $P=Start-Process -FilePath $Py -ArgumentList $CmdArgs -PassThru -NoNewWindow -RedirectStandardOutput "$LogBase.stdout.log" -RedirectStandardError "$LogBase.stderr.log"
  $Started=Get-Date
  while (-not $P.HasExited) {
    $P.Refresh()
    if ($P.WorkingSet64 -gt 8GB -or ((Get-Date)-$Started).TotalSeconds -gt $Seconds) {
      taskkill /PID $P.Id /T /F | Out-Null
      $P.WaitForExit()
      throw 'C-04 resource budget exhausted: incomplete; retain output and logs'
    }
    Start-Sleep -Seconds 2
  }
  $P.WaitForExit()
  if ($P.ExitCode -ne 0) { throw "C-04 failed: $($P.ExitCode); retain every file" }
}
```

パスに空白を含まないチェックアウト/出力名を使う。この関数は単独driver用で、batchを中に入れない。時間切れ時はnative driverが起動したCLI子も含めtaskkillで終了する。メモリ監視は親driverのworking setで、全子孫の厳密な合計測定ではない。小型CLIの追加分を含む予算は上の求解RAMに含める。

Tet4長系列（同じ13工程、水位32固定、全体256/512。新しい長時間物理を作るものではない）：

```powershell
$Long='work/NEW-PC-C04-cause-long'
python "$Study/c04cause_driver.py" --output $Long --initialize --probe target/release/examples/wet_joint_transient_probe.exe --solver target/release/soil-fem-3d.exe
python "$Study/c04cause_binding.py" --controls $Long --capture
foreach ($N in @(256,512)) {
  foreach ($Kind in @('native','adaptive')) {
    Invoke-C04Bounded @("$Study/c04cause_driver.py",'--output',$Long,'--label',"$Kind-m$N-w32",'--mechanical-steps',"$N",'--water-steps','32','--timeout','14400') "$Long/$Kind-m$N-w32"
  }
}
```

PC初期予算は**1経路10〜90分、求解RAM 2GB・一括後処理8GB以下、出力1〜5GB**。12 Tet4・48変位成分で、今回の128までの実測からの参考見積り。256/512は未実測。新たな刻み細分化で失敗を消したり、その値を新しい受入参照にしたりしない。

12 Tet10（54節点・324材料GP・32接触GP）の追加対照：

```powershell
$T10='work/NEW-PC-C04-cause-tet10'
python "$Study/c04cause_driver.py" --output $T10 --initialize --element tet10 --probe target/release/examples/wet_joint_transient_probe.exe --solver target/release/soil-fem-3d.exe
python "$Study/c04cause_binding.py" --controls $T10 --capture
foreach ($N in @(32,64,128)) {
  foreach ($Kind in @('native','adaptive')) {
    Invoke-C04Bounded @("$Study/c04cause_driver.py",'--output',$T10,'--label',"$Kind-m$N-w32",'--mechanical-steps',"$N",'--water-steps','32','--timeout','14400') "$T10/$Kind-m$N-w32"
  }
}
Invoke-C04Bounded @("$Study/c04cause_driver.py",'--output',$T10,'--label','boosted-m32-w32','--mechanical-steps','32','--water-steps','32','--boost-nontrivial') "$T10/boosted-m32-w32"
```

PC初期予算は**1経路0.5〜4時間、求解/逐次後処理RAM 8GB以内、出力2〜10GB**。未実測。GP数27倍・変位成分162という小さい密行列の規模からの保守的予算で、4時間や8GBを超えた経路は未完了として終了する。メモリ40GB超や1経路1日超をPC前提にしない。後処理は一度に全ケースをロードせず2経路ずつ行う。

この2つのPC系列は旧Tet4封印と同じケース数/時刻列ではないため、旧baseline-bindingや固定名の12経路用reference-assessmentをそのまま適用しない。各対の比較例：

```powershell
python -c "import sys;from pathlib import Path;sys.path.insert(0,'$Study');from c04cause_metrics import load,compare;from c04cause_observe import common_pair,audit_global;from c04cause_driver import write;p=Path('$T10');a=p/'adaptive-m32-w32';b=p/'adaptive-m64-w32';write(p/'pair-m32-m64.json',dict(comparison=compare(*common_pair(load(a/'fields.jsonl'),load(b/'fields.jsonl'))),left_audit=audit_global(a,p),right_audit=audit_global(b,p)))"
```

32→64、64→128とnative/boostedについて同様に比較する。Tet4長系列は256→512。本体と参照の各元力基準、共通座標数、散逸・すべり・カテゴリの全場/端点差を見る。Richards湿潤JOINT、一般3Dの連続極限・規模上限を超える検証は**規模上限のため確認不可・未達**（現行未対応の接続は資源量とは別の理由でも未達）。

## 封印の軽い検証と再計測

```powershell
python "$Study/c04cause_evidence.py" --mode verify --source $Study --output work/NEW-C04-cause-sealed --extract
if ($LASTEXITCODE -ne 0) { throw 'new archive SHA/inventory failed' }
python -m unittest discover -s tests -p test_c04_transient_water_cause_diagnosis.py -v
```

part・ZIP・全エントリーSHAを照合し、新規展開する。`raw/data/global-final` が保存求解、`raw/data/observation-final` が全原因集計。旧baselineは最初に照合した `$Base` を使い、新しい出力先でc04cause_observeを再実行できる。再監査・再計測を新規全体求解として数えない。アーカイブ検査の時間/容量はarchive-index.json、全件検証ログを参照する。

## ソフトウェア全件

局所再現用の疎取得では過去NC/C-01等のアーカイブが不足する。全件試験にはCIと同様に、旧postfoldだけを省いた全証拠を配置する。既存NCコード/設定・試験は変更しない。

```powershell
git sparse-checkout set --no-cone '/*' '!/interim/2026-11-10-nc-postfold-branches/evidence.tar.gz.part*'
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release
cargo test --locked --release
foreach ($G in @('python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence','python-c01-switch','python-c01-path','python-c01-cause')) {
  python tools/run_ci_python_evidence_shards.py --group $G --output "work/NEW-C04-software/$G"
  if ($LASTEXITCODE -ne 0) { throw "full Python group failed: $G" }
}
```

Git履歴の取得量・ディスクには別途余裕を取る。全件試験は数十分〜数時間の初期予算で、数値検証の1経路予算とは分ける。全群の発見IDと実行IDの和集合を照合し、欠落・重複・failure/error/未実行/skipを正確に記録する。P6台帳登録後の監査と関連回帰も実行する。push後の新先端Windows CIは未確定。
