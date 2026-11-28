# A-03原因診断の再現と規模上限

PR #13に積む診断で、マージしない。2026-12-22は依頼ラベル。
起点9a5e9d182299b6fce64cc18ea61a65be4702c4e1の子孫を使う。
A-03の元基準と定義は[固定契約](SPEC.ja.md)および旧11-26記録のまま。
未存在の出力先を毎回指定し、以前の失敗/生証拠/封印は上書きしない。

## Windowsの準備

GPD WIN Max 2・使用可能約48GBを前提とする。以下は**PC未実行の手順**。
追加の解析は同じ54節点・12 Tet10だけ。768/6144要素、見込み40GB超または
1経路1日超の新しい手順は作らない。必要な確認がこの上限を超える場合は
**「規模上限のため確認不可・未達」**と記録し、合格条件を変えない。

```powershell
git status --short
git fetch origin codex/plastic-spatial-convergence
git merge-base --is-ancestor 9a5e9d182299b6fce64cc18ea61a65be4702c4e1 HEAD
if ($LASTEXITCODE -ne 0) { throw '指定起点の子孫ではありません' }
$Py = '.venv/Scripts/python.exe'
if (-not (Test-Path $Py)) { py -3.13 -m venv .venv }
& $Py -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw '固定依存の配置失敗' }
$env:CARGO_BUILD_JOBS='2'
$env:RUST_TEST_THREADS='2'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
if ($LASTEXITCODE -ne 0) { throw 'toolchainの配置失敗' }
cargo build --locked --release --all-targets
if ($LASTEXITCODE -ne 0) { throw '本体と既存観測器のビルド失敗' }
$Diag = 'interim/2026-12-22-a03-wet-soil-birth-cause-diagnosis'
```

ビルドは10〜60分・メモリ概算2〜8GBの予算（Windows未実測）。全アーカイブを含む
checkoutと作業出力には十分なディスク空きを確保する。通常core群では既存CIのsparse規則が
新しい証拠partを除外しない。アーカイブ欠落なら元Git blob/SHAで戻し、skipで回避しない。

## 小型の全体平衡と原因観測

```powershell
$env:A03_PROCESS_BUDGET_SECONDS='3600'
& $Py "$Diag/a03cause_driver.py" --counts 2 4 8 --output runs/NEW-a03-cause-PC
if ($LASTEXITCODE -ne 0) { throw '診断未完了。途中状態・失敗ログを保持' }
& $Py "$Diag/a03cause_report.py" --run runs/NEW-a03-cause-PC --output runs/NEW-a03-cause-PC-report.json
if ($LASTEXITCODE -ne 0) { throw '集計失敗。求解成功に読み替えない' }
Remove-Item Env:A03_PROCESS_BUDGET_SECONDS
```

12 Tet10のh=1.5（外水面指定なし）/2.0（exterior明示）、N=2/4/8。
各経路は初期平衡→土BIRTH→JOINT BIRTH/B→上載圧20 kPa。
N=4では同一実行版で8実checkpointの別PID復旧・12再hash改変、旧24拒否と
全GPが水下でも一部節点が水上の追加拒否、GL5/GL7全体平衡、JSONキー順対照を含む。
正式クラウド行列は65.826秒。PCの保守的な計画枠は全体5〜30分、各プロセス上限1時間、ピーク概算1〜4GB。
これは予算であってWindows実測の速度・メモリ保証ではない。実測は各process.jsonで記録する。

|出力|見る値・意味|
|---|---|
|`source-binding.json`, `source-snapshot/`, `binaries.json`, `bin/`|求解前ソースと固定バイナリSHA。解析中に実行版を差し替えない|
|`h*-n*/native/run.json`, `accepted.jsonl`, `constitutive.jsonl`, `preparations.json`|全Run、全受理点、力の分解、準備状態。公開API観測を独立求解に数えない|
|`independent-fields.jsonl`, `independent-prepared.jsonl`, `independent-result.json`|入力だけからの独立全体平衡、初期状態、力残差。GL7は343体積点/要素で出力GPは元の27点|
|`legacy-comparison.json`, `legacy-prepared-comparison.json`|旧比較器をそのまま適用した結果。旧保存値を書き換えない|
|`driver-result.json`, `*-process.json`, `*.log`|実コマンド/PID/終了コード/時間切れ、実際の完走経路数。復旧/監査/拒否は新規経路に加算しない|
|`restarts/`, `tamper/`, `refusals/`|実復旧と再hash内容改変拒否。有効な先行Darcyの完了prefixとBIRTH停止を分ける|
|`NEW-a03-cause-PC-report.json`|固定順の全場差、絶対差/分母、材料状態/順序/ULP、水圧力の各寄与、共通格子、GL感度。passed=falseとnullを保持|

元の変位/応力/traction相対1e-7、全節点力1e-10 kN、拘束1e-10 m、接触仕事恒等式1e-10、
入力の力1e-11/線形1e-12を保持する。全54場を新しいしきい値で包括合格にしない。
比較のGPは所有要素/local index、接触はJOINT/face/固定面内順。最近接点へ付け替えない。
同じRunの比較ではelapsed_secondsだけを除く。標準CLIと観測exampleの版SHA差だけは
別対照として区別し、指紋を付け替えない。全場の符号付き0、null、活動/世代も保持する。

## 任意の長系列：同じ12 Tet10だけ

```powershell
$env:A03_PROCESS_BUDGET_SECONDS='21600'
try {
  & $Py "$Diag/a03cause_driver.py" --counts 16 32 --output runs/NEW-a03-cause-long-PC
  if ($LASTEXITCODE -ne 0) { throw '長系列未完了。ログ・途中出力を保持' }
  & $Py "$Diag/a03cause_report.py" --run runs/NEW-a03-cause-long-PC --output runs/NEW-a03-cause-long-PC-report.json
  if ($LASTEXITCODE -ne 0) { throw '長系列の集計未完了' }
} finally { Remove-Item Env:A03_PROCESS_BUDGET_SECONDS -ErrorAction SilentlyContinue }
```

N=16は各64状態、N=32は各128状態。観測なしと標準CLI対照も各々別プロセスで記録する。
クラウドの観測本体N2/N4/N8は約0.82/2.09/5.27秒で、seal/replayを含む。
PCの保守的な計画枠はN16で各経路1〜15分、N32で3〜60分、メモリ概算1〜6GB。
これは小型の実時間からの参考外挿に広い予算を置いたもので、Windows実測ではない。
seal/resumeによる再演算があるので単純なN比例の外挿ではない。各プロセス6時間で停止したら
数値未完了・未達。1経路1日超の見込みへ手順を延長しない。メモリ見込み40GB超となる変更も行わない。
全体時間が長くても、並列化して48GBを使い切る計画にはしない。

この追加ではN=4専用の復旧/拒否/GL対照を反復しない。N=2/4/8結果と同じ工程の
共通進行度だけを比較する。現在のreportは2/4/8の次数表を作るため、16/32だけの出力では
その表を作らない。16/32の全場ファイルは保存される。必要なら以下で補間なしの差を記録する。

```powershell
@'
import sys,json
from pathlib import Path
p=Path('interim/2026-12-22-a03-wet-soil-birth-cause-diagnosis');sys.path.insert(0,str(p))
from a03cause_report import load
from a03cause_metrics import compare,common
from a03cause_binding import fields
for h in ('h1_5','h2_0'):
 for side in ('native','independent'):
  a=common(load(Path('runs/NEW-a03-cause-long-PC')/(h+'-n16'),side),16)
  b=common(load(Path('runs/NEW-a03-cause-long-PC')/(h+'-n32'),side),16)
  result=compare(a,b,fields.CATEGORIES)
  target=Path('runs')/('NEW-a03-'+h+'-'+side+'-16-32.json')
  with target.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2,allow_nan=False)
'@ | & $Py -
if ($LASTEXITCODE -ne 0) { throw '共通格子の対応/全場差の確認失敗' }
```

見る量は変位/全有効応力/全水圧/力分解/接触履歴/元定義の仕事・数値散逸。
差の減少や観測次数は参考値であり、外挿/新定義だけで合格を作らない。
独立参照は有限離散化の対照で、最細を厳密解へ昇格させない。
Tet10は既に今回の小型系列で使うため、追加の大きなTet10メッシュを要求しない。
一般3D空間収束が上限内で証明できない場合は未達を保持する。

## 全件検証と保存証拠

```powershell
cargo fmt --all --check
if ($LASTEXITCODE -ne 0) { throw 'fmt失敗' }
cargo clippy --locked --all-targets -- -D warnings
if ($LASTEXITCODE -ne 0) { throw 'clippy失敗' }
cargo test --locked --release
Write-Output "Rust main exit=$LASTEXITCODE"
cargo test --locked --release --examples
Write-Output "Rust examples exit=$LASTEXITCODE"
foreach ($group in @('python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence','python-c01-switch','python-c01-path','python-c01-cause')) {
  & $Py tools/run_ci_python_evidence_shards.py --group $group --output "runs/NEW-a03-validation/$group"
  Write-Output "$group exit=$LASTEXITCODE"
}
& $Py "$Diag/a03cause_evidence.py" --verify $Diag --extract runs/NEW-a03-archive-check
if ($LASTEXITCODE -ne 0) { throw 'part/結合ZIP/entryのSHA不一致' }
& $Py interim/2026-12-06-p6-01-status-reconciliation/audit_status.py --help
```

全件は既存の分割と基準のまま。想定はRust全体10〜90分、Python全群30〜120分、
メモリ概算2〜8GB（Windows未実測の予算）。発見/実行/成功/失敗/error/未実行/skipを正確に保存し、
初回失敗と再実行を重複算入しない。文書はuniversal newlineで読み、CRLFの同内容でも一致させる。
ノルムの集計順はC-order/math.fsumに固定し、Windowsの順序差を許容値で吸収しない。
証拠SHA検査/保存場再集計は新規求解に数えない。

Linuxでは `source /workspace/soil-environment/activate.sh` 後に同じPythonコマンドを使える。
P6台帳の最終監査は固定コミットとSHAを確認する。詳細オプションは上記helpと今回保存ログを参照。
当該push後のWindows CIは未確定。Linux成功や旧CI結果を転用しない。
A-01/02/04、B、C-01〜04、NC-01、P6も未達のまま。
