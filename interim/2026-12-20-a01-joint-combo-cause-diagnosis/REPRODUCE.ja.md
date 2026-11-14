# A-01原因診断の再現とPCでの追加確認

判定は **未達、passed=false**。計測・SHA監査・ソフトウェア試験の終了0を、A-01全体の合格へ読み替えない。PR #13の同じブランチを使い、マージしない。2026-12-20は保存先ラベル、実施セッションは2026-11-14（日本時間）。PC実行はしていない。

## 規模と停止条件

ユーザーPCはGPD WIN Max 2、利用可能メモリ約48GB。新しいPC手順は下表の小型ケースだけに限定する。時間・メモリはWindowsで未実測の予算で、最初の1ケースで確認する。元の力1e-8、真の線形1e-10、独立力1e-10、独立弾性相対1e-7、構造並進/回転絶対1e-9、板生値絶対1e-8/相対1e-7、局所/再組立/仕事の旧基準を維持する。

|処理|メッシュ/系列|想定時間（ビルド/取得を除く）|想定メモリ|出力と見る値|
|---|---|---|---|---|
|封印監査・派生集計|旧7経路と今回の小型生証拠|1–10分|1–4GB|全part/entry SHA、元仕事L2、恒等式の残り、同一線形系の残差|
|同じTet4の全体32分割|12 Tet4・16節点・6接触GP、辺共有B。必要なら旧7経路へ展開|1ケース2–15分、旧7経路は1–4時間の予算|2–8GB|`fields.jsonl` / `native-fields.jsonl`、変位/応力/traction/塑性すべり/仕事、全受理座標|
|独立の固定接触1/2/4/8、全体8|同じ12 Tet4|各1–5分|1–4GB|独立の力残差と分割、元仕事差。局所再生と区別|
|既存の小型Tet10対照|12 Tet10・54土節点・324材料GP・32接触GP、辺共有Bまたは部材1種、全体16|1ケース0.5–2時間、最大4ケース2–8時間の予算|8–32GB|全GP/全節点・構造回転/力・仕事・再開と拒否|

768/6144要素、見込み40GB超、1経路1日超の計算は新たにPC手順へ追加しない。これらが必要な一般化は **「規模上限のため確認不可・未達」**。表の予算を超える場合は途中ファイルと停止理由を保持し、上限超過を要求したり、許容値を緩めたりして通さない。支持・拘束の変更は行わない。数値未完了を入力拒否や物理破壊へ変換しない。

## Windowsの準備

Windows 11 / PowerShell、Python 3.13、Rust 1.99.1。パスを毎回新しくする。`work/`を計算出力先に使い、旧interimへ上書きしない。

```powershell
git switch codex/plastic-spatial-convergence
git pull --ff-only
# 全試験に必要な旧証拠を保持する。C-01専用CI群の除外規則はここへ転用しない。
git sparse-checkout set --no-cone '/*' '!/interim/2026-11-10-nc-postfold-branches/evidence.tar.gz.part*'
py -3.13 -m venv .venv-a01-cause
.\.venv-a01-cause\Scripts\Activate.ps1
python -m pip install -r requirements-reference.lock.txt
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
$env:CARGO_BUILD_JOBS='2'
$env:RUST_TEST_THREADS='2'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
$Diag='interim/2026-12-20-a01-joint-combo-cause-diagnosis'
$Old='interim/2026-11-23-joint-combo-extension'
cargo build --locked --release --bin soil-fem-3d --example joint_combo_probe --example a01_joint_cause_probe --example a01_joint_linear_replay
if ($LASTEXITCODE -ne 0) { throw 'A-01 build failed' }
```

取得・全件検査には旧証拠を含む相応のディスクが必要。必要partが欠けた場合は元Git blobを取得し元SHAで照合する。欠落をskip、manifestの更新、別の参照値への置換で回避しない。新しいpartは通常のPython core群で除外されない。全件検査の前にP6の固定コミットも取得する（shallow cloneなら `git fetch --unshallow`）。

## 封印の監査と旧経路の再集計

```powershell
python "$Diag/a01cause_evidence.py" --verify $Diag
if ($LASTEXITCODE -ne 0) { throw 'new evidence SHA failed' }
python "$Diag/a01cause_evidence.py" --extract-old $Old --output work/NEW-A01-old-verified
if ($LASTEXITCODE -ne 0) { throw 'old evidence SHA failed' }
```

旧5 part・2,620エントリーを全検査し、`selected/`の入力・全場・実分割と平行2面の旧途中結果を新規出力先へ展開する。旧checkpoint/restartを消去・変更する処理ではない。新アーカイブの `formal/` が正式求解、`retained-preparation/` が準備・失敗、`old-verified-subset/` が監査済み旧保存値。新アーカイブの全entryは `archive-index.json` に固定してある。実行版は同梱せず、SHAとソースを保持する。

分割ZIPを新規フォルダーへ展開する場合は先に上の全SHA検査を通す。Pythonの`zipfile`でpartを番号順に一時ファイルへ結合して展開する。任意の保存パスを変更せず、`formal/linear/replay-input.json`、`replay.json`、`stacked/native/trace.jsonl`と入力3JSONを一組として読む。

Linuxで実施した正式行列をそのまま新規求解する一括入口は次のとおり。`a01cause_batch.py`の固定バイナリ名はLinux用なので、Windowsでは次節の個別入口を使う。

```bash
source /workspace/soil-environment/activate.sh  # このCloud環境の設定。別環境では同じ固定版/変数を設定
cargo build --locked --release --bin soil-fem-3d --example joint_combo_probe \
  --example a01_joint_cause_probe --example a01_joint_linear_replay
python interim/2026-12-20-a01-joint-combo-cause-diagnosis/a01cause_batch.py \
  --output work/NEW-A01-cause-matrix --baseline /PATH/TO/5908506/soil-fem-3d
```

`--baseline`には起点5908506から別の小型ソースcheckoutでビルドした旧CLIを指定する。古い実行版のSHAを新しいものへ付け替えない。比較時に除くのは計測時間とその実行版固有SHAだけで、入力SHA、数値、履歴は残す。

一括行列は本体5＋観測なし5、独立14の24完走経路と、平行2面4未完了経路を記録する。各子プロセスの時間上限900秒は未完了を検出するためのもので、タイムアウトを成功にしない。Cloudの正式本体選択経路は約4.6–25.2秒、各独立対照は約0.4–0.9秒だったが、Windowsの実測ではない。

## Windows PC：小型長系列と局所分割対照

同じ入力から全体32分割を独立に求める。本体の受理列・場・実分割数を独立側の入力へ渡さない。

```powershell
python "$Diag/a01cause_driver.py" --mode native --case shared-stress_free --steps 32 --output work/NEW-PC-A01-native-g32 --probe target/release/examples/joint_combo_probe.exe
if ($LASTEXITCODE -ne 0) { throw 'native incomplete; retain partials' }
python "$Diag/a01cause_driver.py" --mode controls --case shared-stress_free --steps 32 --output work/NEW-PC-A01-independent-g32
if ($LASTEXITCODE -ne 0) { throw 'independent incomplete; retain partials' }
foreach ($parts in @(1,2,4,8)) {
  python "$Diag/a01cause_driver.py" --mode controls --case shared-stress_free --steps 8 --parts $parts --output "work/NEW-PC-A01-fixed$parts-g8"
  if ($LASTEXITCODE -ne 0) { throw 'fixed-part global control incomplete' }
}
```

本体個別入口の子プロセスは900秒上限。上限到達は失敗/未完了として保存する。より長い旧7経路の比較は既存入口で、時間上限を明示して行う。

```powershell
python "$Old/study.py" --output work/NEW-PC-A01-tet4-g32 --steps 32 --probe target/release/examples/joint_combo_probe.exe --solver target/release/soil-fem-3d.exe --timeout 14400
if ($LASTEXITCODE -ne 0) { throw 'Tet4 matrix incomplete; retain every partial file' }
python "$Old/diagnose.py" --output work/NEW-PC-A01-tet4-g32 --steps 32
if ($LASTEXITCODE -ne 0) { throw 'saved-state audit incomplete' }
```

出力で確認する値は `native/run.json.complete`、全受理座標 `i/32`、元の力/線形残差、`independent-result.json.complete`、全GPの仕事・すべり・traction、世代の0開始と凍結、別PID再開。元の比較器の基準を保持する。今回の固定順ノルムで追加比較する時は `a01cause_metrics.compare` を使い、共通入力座標へ抽出し、成分の順序を固定する。積分定義・正規化を変えて合格にしない。

## Windows PC：12 Tet10の限定した範囲確認

まず辺共有Bの1ケースだけを実測する。その後も表の時間/メモリ上限内なら、同じ既存12 Tet10に梁/板/アンカーを1種ずつ追加した最大4ケースへ展開する。2柱複製のTet10や、大型メッシュへ拡大しない。

```powershell
python "$Old/study.py" --output work/NEW-PC-A01-tet10-shared-g16 --steps 16 --quadratic --user-pc --case shared-stress_free --probe target/release/examples/joint_combo_probe.exe --solver target/release/soil-fem-3d.exe --timeout 14400
if ($LASTEXITCODE -ne 0) { throw 'small Tet10 incomplete; retain partials' }
python "$Old/diagnose.py" --output work/NEW-PC-A01-tet10-shared-g16 --steps 16
if ($LASTEXITCODE -ne 0) { throw 'Tet10 saved-state audit incomplete' }
# 必要なら別出力先で、--case joint-beam / joint-plate / joint-anchor を各1ケースずつ実行する。
```

Tet10でも適合小メッシュの対照に留まる。先端付き湿潤BIRTH、複数先端、先端＋構造の拒否を回避する入力は作らない。これらは独立平衡・投影・世代・仕事の仕様/受入が揃うまで未対応・未達。

## 平行2面の観測を再現する

```powershell
python "$Diag/a01cause_driver.py" --mode stacked --steps 8 --output work/NEW-PC-A01-observations/stacked --probe target/release/examples/a01_joint_cause_probe.exe
if ($LASTEXITCODE -ne 0) { throw 'observation did not finish; retain partials' }
python "$Diag/a01cause_linear.py" --input work/NEW-PC-A01-observations/stacked --output work/NEW-PC-A01-linear
if ($LASTEXITCODE -ne 0) { throw 'same-system analysis failed' }
target/release/examples/a01_joint_linear_replay.exe work/NEW-PC-A01-linear/replay-input.json work/NEW-PC-A01-linear/replay.json
if ($LASTEXITCODE -ne 0) { throw 'linear replay recording failed' }
python "$Diag/a01cause_modes.py" --formal work/NEW-PC-A01-observations --output work/NEW-PC-A01-mode-localization.json
```

`a01cause_modes.py`は同じ `work/NEW-PC-A01-observations/stacked/` 配置の実行列を読む。`a01cause_nonlinear.py`も同じ配置と監査済み旧subsetを引数に取り、全Newtonの開閉/剛体支持診断と実分割の局所監査を行う。

見る値は、観測/非観測のRun一致、`complete=false`と停止理由、予測子前/後のA,b、最小特異値、中央土塊のモード局在、直接LUと元GMRESの**同じA,bの真残差**、全line search候補。新しい解や支持を本体へ戻さない。独立BIRTHの失敗を元本体openの失敗と同一段階の失敗にはしない。小型で残る原因分離はCloud残件であり、PCへ必須移管したわけではない。

## 全件試験とP6の再監査

```powershell
cargo fmt --all --check
if ($LASTEXITCODE -ne 0) { throw 'fmt failed' }
cargo clippy --locked --all-targets -- -D warnings
if ($LASTEXITCODE -ne 0) { throw 'clippy failed' }
cargo build --locked --release
if ($LASTEXITCODE -ne 0) { throw 'build failed' }
foreach ($group in @('rust-core','rust-integration-0','rust-integration-1')) {
  python tools/run_ci_rust_tests.py --group $group --output "work/NEW-PC-A01-validation/$group"
  if ($LASTEXITCODE -ne 0) { throw 'Rust group failed' }
}
cargo test --locked --release --examples
if ($LASTEXITCODE -ne 0) { throw 'example tests failed' }
foreach ($group in @('python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence','python-c01-switch','python-c01-path','python-c01-cause')) {
  python tools/run_ci_python_evidence_shards.py --group $group --output "work/NEW-PC-A01-validation/$group"
  if ($LASTEXITCODE -ne 0) { throw 'Python group failed or incomplete' }
}
python interim/2026-12-06-p6-01-status-reconciliation/audit_status.py --root . > work/NEW-PC-A01-validation/p601.json
if ($LASTEXITCODE -ne 0) { throw 'P6 inventory audit failed' }
```

P6監査は同入口の `--root` と標準出力のリダイレクトを使い、新しい出力先へ保存する。固定コミットのGit blob、台帳SHA、型付きbool/nullを検査し、`integrity_ok=true`とSOURCE_NOT_INVENTORIED=0を確認する。試験件数は各inventoryと全発見IDから成功/失敗/未実行/skipを集計し、単独スモークや再確認を二重加算しない。旧NC封印の厳密比較で失敗した場合も、旧保存値・許容値・NCコードを変えず失敗として残す。

Windows全件は旧証拠の再監査を含み、ビルド/取得込みで数時間になり得る。科学系列の1経路1日超を要求する手順とは分ける。push後の新先端Windows CIは未確定で、Linux結果や旧CI成功を転用しない。
