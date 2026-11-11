# A-01の再現：Linux／Windows PC

判定は **未達、passed=false**。再現スクリプトの終了コード0は計測・比較・再開・拒否確認を完了した意味であり、全履歴の数値受入を意味しない。保存ラベル2026-11-23、実施2026-11-12（日本時間）。PR #13の同じブランチを使い、マージしない。旧証拠・失敗・参照値・許容値・NC設定は変更しない。

## Windowsの準備

Windows 11 / PowerShell、Python 3.13、Rust 1.99.1。小型行列はCPU 4〜8コア、RAM 16 GBで開始できる。全件試験では旧封印証拠も読み、取得と一時展開用に空き80〜100 GB程度を計画する。初回ビルドと全Rust試験のコンパイルは30〜90分程度、全件試験の旧証拠取得は数GB以上で回線に依存する。Windowsで今回の新先端を実行した時間・CI結果は未確認である。

```powershell
git switch codex/plastic-spatial-convergence
git pull --ff-only
# 既存CIと同じ。試験に使わない旧postfold生アーカイブだけを省く。
git sparse-checkout set --no-cone '/*' '!/interim/2026-11-10-nc-postfold-branches/evidence.tar.gz.part*'
py -3.13 -m venv .venv-a01
.\.venv-a01\Scripts\Activate.ps1
python -m pip install -r requirements-reference.lock.txt
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
$env:CARGO_BUILD_JOBS='2'
$env:RUST_TEST_THREADS='2'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
cargo build --locked --release --bin soil-fem-3d --example joint_combo_probe
if ($LASTEXITCODE -ne 0) { throw 'A-01 build failed' }
cargo build --locked --release --example joint_combo_structure_loads
if ($LASTEXITCODE -ne 0) { throw 'structure load observer build failed' }
$Study='interim/2026-11-23-joint-combo-extension'
$Run='work/NEW-A01-small'
```

`git sparse-checkout`で旧証拠を取得する。Python全件を開始する前に必要な全partが存在することを確認する。今回Linuxの全件検査に必要な4アーカイブ（NC mean/spatial、湿潤convergence、step）の合計は約7.47 GB。取得欠落、0バイトの途中ファイル、SHA不一致をskipやmanifest更新で回避しない。必要ならGitから同じパスを復元し、元のSHAで確認してから全群を実行する。

## 小型の入力生成・両側新規求解・全点比較

```powershell
python "$Study/study.py" --output $Run --steps 8 --probe target/release/examples/joint_combo_probe.exe --solver target/release/soil-fem-3d.exe --timeout 900
if ($LASTEXITCODE -ne 0) { throw 'matrix incomplete; retain every partial file' }
python "$Study/diagnose.py" --output $Run --steps 8
if ($LASTEXITCODE -ne 0) { throw 'local/nodal diagnostic incomplete' }
python "$Study/geometry_audit.py" --matrix $Run --report work/NEW-A01-geometry.json
if ($LASTEXITCODE -ne 0) { throw 'geometry observation differs' }
python "$Study/structure_work.py" --matrix $Run --output work/NEW-A01-work --producer target/release/examples/joint_combo_structure_loads.exe
if ($LASTEXITCODE -ne 0) { throw 'whole-system external work observation incomplete' }
python "$Study/stacked_attempt.py" --output work/NEW-A01-stacked --steps 8 --solver target/release/soil-fem-3d.exe --timeout 900
if ($LASTEXITCODE -ne 0) { throw 'stacked attempt was not recorded' }
```

出力先は毎回未存在の名前にする。途中結果を消して同じ名前で再試行しない。各行列は実行版を`_bin/`に固定コピーし、ソース指紋と新スナップショットを保存する。ソースが途中で変わった場合も失敗とする。

クラウドの小型入力は各工程8増分・全7経路・672受理対。Windowsの初期予算は行列15〜45分程度（ビルド・取得・全件試験を除く）。各CLI900秒の上限に達した場合は未完了として保存する。小型の実測時間は生の各`*-process.json`に記録し、PCの予算を実測値に置き換える。意図的な範囲拒否18件は別扱いであり、時間切れや未知のエラーを既知の拒否へ変換しない。

相互作用するstackedケースは現状の入力で本体・独立とも未完了となる。`native.json`の`complete=false`、`native.log`の真の残差、`independent-result.json`の失敗を確認する。記録スクリプトの正常終了を解析合格としない。中間土塊に支持を追加したり、許容値を変えたりして同じケースの成功にしない。

Linuxでは仮想環境を有効にし、実行版の`.exe`を省く。

```bash
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
export CARGO_BUILD_JOBS=2 RUST_TEST_THREADS=2 PYTHONUTF8=1 PYTHONIOENCODING=utf-8
cargo build --locked --release --bin soil-fem-3d --example joint_combo_probe
cargo build --locked --release --example joint_combo_structure_loads
python interim/2026-11-23-joint-combo-extension/study.py \
  --output work/NEW-A01-linux --steps 8 --timeout 900
python interim/2026-11-23-joint-combo-extension/diagnose.py \
  --output work/NEW-A01-linux --steps 8
python interim/2026-11-23-joint-combo-extension/structure_work.py \
  --matrix work/NEW-A01-linux --output work/NEW-A01-work-linux
python interim/2026-11-23-joint-combo-extension/stacked_attempt.py \
  --output work/NEW-A01-stacked-linux --steps 8
```

## 出力と合否基準

|ファイル|確認すること|
|---|---|
|`identity.json`、`source-fingerprints.json`、`source-snapshot.zip`|指定起点、固定コピーした本体/観測CLIのSHA、入力格子、実行時ソース。旧hash/snapshotと別物|
|各経路の`mesh.json`、`ground.json`、`case.json`|既存入力機能・材料値、BIRTH A/B、物理境界条件。`case.json`だけから受理格子`i/8`が定まる|
|`native/accepted.jsonl`、`constitutive.jsonl`、`structure-reports.jsonl`|本体の全確定状態、実内部力/target/分割、構造の生レポート。各96状態、切戻し0、欠落0|
|`independent-fields.jsonl`、`independent-result.json`|入力だけから新規求解。全96状態、`reference_uses_native_results=false`。保存JSON監査とは区別|
|`native-fields.jsonl`、`comparison.json`|全受理点・全場、全成分の最大絶対差/プールL2、各工程、カテゴリ不一致数。0誤差や適当な許容値を補完しない|
|`native/run.json`、`whole/run.json`|両方complete=trueで、全Runがelapsed_secondsだけを除き完全一致。符号付き0、全入力指紋、全世代を保持|
|`native/checkpoint-0012.json`、`0028`、`0076`|birth/open/rebirthのprogress=.5に実ファイルがあること。別PIDが同じSHAの実行版で読む|
|`restarts/*/accepted.jsonl`、`restarts/verification.json`|各84/68/20後続全状態、全7経路合計1,204状態と全Runが時間以外完全一致。許容差0|
|`cause-diagnostics.json`|本体jump/本体分割での局所診断、全節点再組立、仕事収支、世代の厳密0開始。新規3D求解数へ数えない|
|`refusals/verification.json`、各`cli.log`|11未対応組合せと7改変の既知メッセージ、終了コード1、結果ファイルなし|
|`failure-*.json`、`completion.json`|未知エラー、欠落、途中、時間切れは失敗。`measurement_complete=true`でもA01/passedは未達/false|
|stackedの`verification.json`|スキーマ拒否・数値未完了・物理失敗を混同しない。両側の途中状態と理由を保存|

既存の非線形平衡1e-8、真の線形残差1e-10、独立平衡1e-10、拘束1e-10 m、全節点再組立1e-10 kN、接触仕事収支1e-10、独立弾性相対1e-7、構造並進/回転絶対1e-9、板生値絶対1e-8/相対1e-7を維持する。新規全場数値ゲートを観測値から作らない。今回の変位/応力/構造差と接触仕事差が残る限り、**A-01の各全履歴は未達**。世代・再開・カテゴリの一致だけで限定受入へ変えない。

`package.py`は今回の96状態×7経路を封印するためのツールである。既存成果物へ再出力せず、新しい保存先に対して使う。長いPC系列や異なるケース数の受入器ではない。

## PC：長い履歴とTet10

次の系列はユーザーPC用で、今回クラウドで完了したとはしない。全体格子だけを長くし、材料/接触の入力上限・物理・許容値は変えない。C-01〜C-04の収束課題の代行でもない。

```powershell
python "$Study/study.py" --output work/NEW-PC-A01-tet4-g32 --steps 32 --probe target/release/examples/joint_combo_probe.exe --solver target/release/soil-fem-3d.exe --timeout 14400
if ($LASTEXITCODE -ne 0) { throw 'Tet4 long history incomplete' }
python "$Study/diagnose.py" --output work/NEW-PC-A01-tet4-g32 --steps 32
if ($LASTEXITCODE -ne 0) { throw 'Tet4 diagnosis incomplete' }

# 12 Tet10・54土節点・324材料GP・32接触GPの既存メッシュ。
# 2柱を複製するseparatedはこの全場観測のGP上限を超えるので選ばない。
python "$Study/study.py" --output work/NEW-PC-A01-tet10-g16 --steps 16 --quadratic --user-pc --case shared-stress_free --case joint-beam --case joint-plate --case joint-anchor --probe target/release/examples/joint_combo_probe.exe --solver target/release/soil-fem-3d.exe --timeout 14400
if ($LASTEXITCODE -ne 0) { throw 'Tet10 composite history incomplete' }
python "$Study/diagnose.py" --output work/NEW-PC-A01-tet10-g16 --steps 16
if ($LASTEXITCODE -ne 0) { throw 'Tet10 diagnosis incomplete' }
```

Tet4の32増分全7経路は初期予算1〜4時間、Tet10の16増分4経路は2〜8時間、RAM 16〜32 GB、出力数百MB〜数GBを見込む（未実測のWindows予算）。観測CLIは受理点ごとに本体checkpointを認証しながら継続するため、増分数に比例するだけの時間ではない。各子プロセス4時間の上限は未完了を検出するためで、超過を成功にしない。まず各1ケースで実測し、必要なら別の未存在出力先と明示した時間予算で継続する。

Tet10も既存の適合小メッシュであり、一般の大規模メッシュ、GUI/GPU、実CAD、交差するJOINT、先端BIRTHの湿潤/構造の新しい入力仕様は用意していない。本体の実内部力診断は512材料GP/256縮約DOFに制限され、この観測exampleは物理節点DOFにも小型上限を置く。より大きい実CADについてはユーザーPCで有効なメッシュ・JOINT宣言・独立参照・観測方法を別途揃え、同じ全場/再開/拒否基準で検証するまで未達とする。未定義のCADやGPUオプションのコマンドは創作しない。

## Linux／Windows全件試験

```powershell
cargo fmt --all --check
if ($LASTEXITCODE -ne 0) { throw 'fmt failed' }
cargo clippy --locked --all-targets -- -D warnings
if ($LASTEXITCODE -ne 0) { throw 'clippy failed' }
cargo test --locked --release 2>&1 | Tee-Object work/NEW-A01-rust-tests.log
if ($LASTEXITCODE -ne 0) { throw 'Rust tests failed' }
cargo build --locked --release --bin soil-fem-3d --example joint_combo_probe
if ($LASTEXITCODE -ne 0) { throw 'CLI build failed' }
foreach ($Group in @('python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence')) {
  python tools/run_ci_python_tests.py --group $Group --output "work/NEW-A01-check/$Group"
  if ($LASTEXITCODE -ne 0) { throw "Python group failed: $Group" }
}
```

Linuxでは同じcargoコマンド、同じ4つのPython群を実行する。全群のinventoryで発見集合と実行集合の一致・重複0を確認し、testsだけでなくtoolsの試験も含める。failure/error/skip/ignored/filtered-out/expected-failure/unexpected-successは0を要求する。途中の単独スモーク、失敗、キャンセル、時間切れは最終成功件数へ加算しない。今回の確定Linux件数・ログ・実行環境は`validation.json`を参照し、push後のWindows CI成功を推測で補わない。

今回のLinuxではNC mean/spatialの旧封印証拠2群が厳密比較で失敗しており、全件成功は未達である。`--allow-platform-float-differences`をLinuxに付けて成功へ変えたり、旧JSONを書き換えたり、当該群を収集から外したりしない。再現時も失敗なら、元SHA・実行環境・差を保存し、NCの実装変更なしで報告する。詳細はREADME末尾とvalidation記録を参照。
