# B-02原因診断の再現（ユーザーPCでは未実行）

対象はGPD WIN Max 2、使用可能約48GB。新規系列は12/96 Tet10と各工程4/8/16刻みだけ。
768/6144、40GB超、1経路1日超の新しいPC計画は作らない。
元の支持・荷重・接続・許容値・受入条件を維持する。
元の12-09手順に書かれた6144追試は今回のPC上限では実行対象にしない。
それが必要な確認は「規模上限のため確認不可・未達」。

## Windowsの準備

PR #13の本作業を含む固定先端をcheckoutし、マージしない。
PowerShell 7、Rust 1.99.1、Python 3.13とrequirements-reference.lock.txtを使用する。
既存のプロジェクト環境がある場合はそれを使う。以下の出力先は未存在の新しい名前とする。

```powershell
$Diag='interim/2026-12-24-b02-offset-composite-cause-diagnosis'
$Py='.venv/Scripts/python.exe'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:CARGO_BUILD_JOBS='2'
$env:RUST_TEST_THREADS='2'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
cargo build --locked --release --all-targets
if ($LASTEXITCODE -ne 0) { throw 'build未完了、ログを保持' }
& $Py interim/2026-12-01-tri6-offset-composite-b02/build_observer.py work/NEW-b02cause-observer
if ($LASTEXITCODE -ne 0) { throw '既存公開API観測器のbuild未完了' }
```

初回buildの計画枠は10〜60分、メモリ概算2〜8GB。Windowsでは未実測。
生産srcやCargoターゲットを変更する手順ではなく、旧観測器の隔離manifestをwork/に作る。

## 12/96・4/8/16の全体系（PC未実行）

```powershell
& $Py "$Diag/b02cause_driver.py" --output runs/NEW-b02cause-PC `
  --observer target/release/b02-observer.exe --cli target/release/soil-fem-3d.exe
if ($LASTEXITCODE -ne 0) { throw '診断未完了。途中場・失敗ログを保持' }
& $Py "$Diag/b02cause_report.py" --run runs/NEW-b02cause-PC --output runs/NEW-b02cause-PC-report.json
if ($LASTEXITCODE -ne 0) { throw '再集計失敗。保存場を変更しない' }
```

本体6経路＋独立6経路、観測なし6、土GL7/積分点逆順各2、同入力再実行/キー順本体・独立各2を実行する。
再開、内容改変、範囲拒否は別記録。刻み8/16でも同じ12/96格子を保つ。
1子プロセス上限7200秒、全体計画枠10〜120分、1経路1秒〜10分、メモリ概算1〜6GB。
これは小規模の計画見積りで、PC実測ではない。時間切れを完走/破壊/受入にせず、上限延長で1日超へ回さない。
元参照のDOF上限やガードを解除するコマンドはない。

見る値は、plan-results.jsonの各経路完走/既存ゲート・受理数・実PID、
レポートのinput_only_vs_native（u/σ'/回転/設置参照）、reference_sensitivity（GL5→7/加算逆順）、
load_step_sensitivity（状態/Q/瞬間仕事/元の累積仕事）、spatial_12_to_96、
maximum（膜の連成恒等式/荷重仮想仕事/全節点再組立/明示構造自由度）、
scope・mutations・restartsである。
force-work-energy-raw.jsonには土・板・梁の符号付き全自由度ベクトル、仕事源、各板の交差/腕二乗項を残す。
work-raw.jsonは旧端点台形規約の値。新しい回復法・仕事定義・ゲートとして採用しない。

元ゲートは再組立1e-10、力1e-8、真線形1e-10、土相対1e-7、構造絶対1e-9、
板生値atol1e-8/rtol1e-7、付着1e-10 m。
診断量の一致や負次数の改善だけでB-02全体をpassed=trueにしない。
拒否先の全体平衡はnull。局所エネルギー/局所再生だけで代用しない。
2空間水準から漸近次数を算出しない。最細の独立参照も厳密連続解ではない。

## 封印証拠の読込と再集計（新規求解なし）

```powershell
& $Py "$Diag/b02cause_evidence.py" --verify $Diag --extract runs/NEW-b02cause-unsealed
if ($LASTEXITCODE -ne 0) { throw 'part/ZIP/member SHA不一致' }
& $Py "$Diag/b02cause_report.py" --run runs/NEW-b02cause-unsealed/formal --output runs/NEW-b02cause-remeasured.json
if ($LASTEXITCODE -ne 0) { throw '再集計失敗' }
```

想定1〜10分、メモリ概算0.5〜3GB。固定順math.fsumの観測値を再現する。
旧numpy集計を上書きせず、Windowsの差が出た場合も許容値・参照値を緩めず記録する。
実行版・入力・全状態のSHAを確認し、異なるバイナリのSHA拒否を内容改変拒否の代わりにしない。

## 全件検証と状態監査

sparse checkoutを使う場合も、新しいevidence.zip.part*と全既存試験が要求するアーカイブを含める。
Linuxクラウドで使った復元記録はvalidation/に保存する。元Git blobとの不一致を無視しない。

```powershell
cargo fmt --all --check
Write-Output "fmt exit=$LASTEXITCODE"
cargo clippy --locked --all-targets -- -D warnings
Write-Output "clippy exit=$LASTEXITCODE"
cargo test --locked --release
Write-Output "Rust main exit=$LASTEXITCODE"
cargo test --locked --release --examples
Write-Output "Rust examples exit=$LASTEXITCODE"
foreach ($group in @('python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence','python-c01-switch','python-c01-path','python-c01-cause')) {
  & $Py tools/run_ci_python_evidence_shards.py --group $group --output "runs/NEW-b02cause-validation/$group"
  Write-Output "$group exit=$LASTEXITCODE"
}
& $Py interim/2026-12-06-p6-01-status-reconciliation/audit_status.py
if ($LASTEXITCODE -ne 0) { throw 'P6台帳integrity不合格' }
```

Rust10〜90分、Python全群30〜120分、メモリ概算2〜8GBの計画枠。Windows未実測。
全発見ID・割当・実行を照合し、成功/失敗/error/未実行/skip/ignore/filterを別記する。
初回失敗は保持し、再実行を同じ固有IDへ重複加算しない。
文書末尾検査はuniversal newlineで読み、実CRLFで試験する。
ノルム順序の違いを許容値で吸収しない。

Linuxはsource /workspace/soil-environment/activate.sh後に同じPythonコマンドを使う。
push後Windows CIは未確定。B-02・他の未達の判定を全件試験の成功へ付け替えない。
