# C-03 原因診断の再現と PC 手順

C-03 は未達、passed=false。保存先の2026-12-18は作業系列ラベル。クラウドで小型を実行し、以下の PC 計算は実行していない。元の許容値・入力拒否・旧結果を変更しない。途中失敗・時間切れ・自動切戻しを保存し、新しい出力名で続ける。

## 環境と取得

PR #13 の `codex/plastic-spatial-convergence` を使う。マージしない。Python 3.13、Rust 1.99.1、固定依存版。全件試験は旧アーカイブと履歴が必要。空き80 GB程度を計画する。通常CIと同じく未使用postfold本体だけを除外し、`reference400-fixture.zip` を含む他の証拠は取得する。

```powershell
git switch codex/plastic-spatial-convergence
git pull --ff-only
git fetch --unshallow # shallow checkout の場合だけ。通常 clone では不要。
git sparse-checkout set --no-cone '/*' '!/interim/2026-11-10-nc-postfold-branches/evidence.tar.gz.part*'
py -3.13 -m venv .venv-c03
.\.venv-c03\Scripts\Activate.ps1
python -m pip install -r requirements-reference.lock.txt
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
$env:CARGO_BUILD_JOBS='2'
$env:RUST_TEST_THREADS='2'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
$Study='interim/2026-12-18-c03-plastic-bulk-cause-diagnosis'
cargo build --locked --release --example c03_bulk_cause_probe --example c03_bulk_return_trace
if ($LASTEXITCODE -ne 0) { throw 'probe build failed' }
cargo build --locked --release
```

上記 fetch 条件以外のコマンドは各非0終了時に止める。既存出力を削除・上書きしない。Linux は `.exe` を除く。新 example の `--diagnose` が Newton 観測の明示スイッチで、既定はオフ。生産CLI・旧example・入力形式には新しいスイッチを追加していない。

## クラウドと同じ小型行列

```powershell
$Base='work/NEW-C03'
python "$Study/c03diag_driver.py" --mode matrix --output "$Base/matrix" --solver target/release/examples/c03_bulk_cause_probe.exe --steps 32 --element tet4 --timeout 900
if ($LASTEXITCODE -ne 0) { throw 'matrix driver failed; keep evidence' }
```

12 Tet4 / 16節点 / 12材料GP / 6接触GP、水位.75、MC c=30、open=.12、全7工程。材料上限 `1e-4,5e-5,2.5e-5` と接触上限 `1e-4,5e-5,2.5e-5` の9組合せを本体・独立で新規求解する。材料m0/m2では接触1.25e-5も測り、m2c2だけ全体64を測る。分割1対照も保持する。独立固定分割は材料8/32×接触4/16の4経路。各端点で必ず全体平衡を解く。入力のみの独立系へ本体の変位・状態・実分割数を渡さない。

本体13、独立17（固定4含む）。通常経路224点、全体64は448点。`matrix-terminal.json` の全終了コードと各 `*-terminal.json` を必ず見る。ドライバの終了コードだけで完了を判断しない。基準の各7工程途中から、同じ実行版SHAの実checkpointを別PIDで復旧し、全Runと784後続状態を時計以外完全一致で照合する。観測なしwholeも同じ全Runで照合する。再開は新規独立求解に数えない。

固定分割は弾性GPにも8/32分割を課すため、上限方式より大幅に遅い。小型の全行列は15〜45分を初期予算とする（PC未実測）。`--timeout` は失敗保存用の上限であり、成功保証ではない。

```powershell
python "$Study/c03diag_scope.py" --matrix "$Base/matrix" --output "$Base/scope" --cli target/release/soil-fem-3d.exe
python "$Study/c03diag_mcc_driver.py" --case "$Base/scope/mcc-pc200" --mesh "$Base/matrix/mesh.json"
python "$Study/c03diag_mcc_driver.py" --case "$Base/scope/mcc-pc105" --mesh "$Base/matrix/mesh.json"
python "$Study/c03diag_mcc_open.py" --matrix "$Base/matrix" --source "$Base/scope/mcc-pc105" --output "$Base/scope/mcc-pc105-open1" --cli target/release/soil-fem-3d.exe
python "$Study/c03diag_mcc_driver.py" --case "$Base/scope/mcc-pc105-open1" --mesh "$Base/matrix/mesh.json"
python "$Study/c03diag_stops.py" --case "$Base/scope/mc-open1" --output work/NEW-C03-open1-stops.json
python "$Study/c03diag_measure.py" --matrix "$Base/matrix" --scope "$Base/scope" --output work/NEW-C03-measurements.json
python "$Study/c03diag_replay.py" --matrix "$Base/matrix" --output work/NEW-C03-local-replay
New-Item -ItemType Directory -Path "$Base/reports"
python "$Study/c03diag_branchtrace.py" --matrix "$Base/matrix" --solver target/release/examples/c03_bulk_return_trace.exe --output "$Base/return-trace" --publish "$Base/reports"
```

`c03diag_stops.py` は同じ新しい親フォルダ下の `matrix/` と `scope/` の配置を使う。入力を移動する場合も、この関係を保持する。

MCCはCLI受理・固定列・観測raw完了を先に確かめる。標準MCCのpc=200/105（open=.12）、追加pc=105（open=1.0）、Backward Euler、同じ12 Tet4・固定水位・CONTINUOUSのみ。混合MCCの拒否は期待する保存証拠であり、受理へ変更しない。独立MCCは既存Backward Eulerの局所離散方程式を2未知数へ消去して解く別実装。塑性仕事・相当塑性ひずみは現行MCC状態にないので捏造しない。pcを含む実在する全場を比較する。

見る値：元のL2と分母・絶対L2、全224/448点と7端点、工程別・GP別、差ベクトルの方向、実分割ヒストグラム、MC return種別、MCC塑性ひずみ/pc、水圧、反力、全体残差、切戻し、実再開完全一致。`local-replay/` は保存経路に沿う局所再生であり、全体平衡の代用ではない。局所再生の比率を加算可能な因果割合にしない。

## PC：Tet10、長系列、96/768/6144規模

新しいPC専用ドライバは既存Tet10細分化を利用し、初期状態を入力から作り直す。大規模ではNewton観測の既存512材料GP上限を外さず、観測コピーをオフにして全受理確定状態だけ保存する。dense独立系の上限は明示 `--user-pc` のときだけ15000変位DOF。本体の元許容値や接触則は同じ。

まず12 Tet10（324材料GP / 32接触GP）を別フォルダで実行する。RAM32 GB・空き20 GBを初期予算にする。

```powershell
$P='work/NEW-PC-C03-tet10'
python "$Study/c03diag_pc.py" --user-pc --mode init --output $P --solver target/release/examples/c03_bulk_cause_probe.exe --refine 0
foreach ($n in 32,64,128,256) {
  foreach ($c in '0.0001','0.00005','0.000025','0.0000125') {
    $Tag="g$n-m0-c$c"
    python "$Study/c03diag_pc.py" --user-pc --mode native --output $P --folder "$Tag-native" --steps $n --material 0.0001 --contact $c --timeout 14400
    if ($LASTEXITCODE -ne 0) { throw 'native failed; retain partial evidence' }
    python "$Study/c03diag_pc.py" --user-pc --mode independent --output $P --folder "$Tag-independent" --steps $n --material 0.0001 --contact $c --timeout 14400
    if ($LASTEXITCODE -ne 0) { throw 'independent failed; retain partial evidence' }
    python "$Study/c03diag_pc.py" --user-pc --mode compare --output $P --folder "$Tag-native" --other "$Tag-independent"
  }
}
python "$Study/c03diag_pc.py" --user-pc --mode restart --output $P --folder 'g32-m0-c0.0001-native' --steps 32 --timeout 14400
```

12 Tet10の独立1経路0.5〜4時間、全系列は数日を予算とする（未実測）。材料上限2.5e-5に固定して同じ接触系列を別タグで実行し、材料×接触相互作用を確認する。長系列だけを小型Tet4で先に測るなら `c03diag_driver.py --mode matrix --element tet4 --steps 64`、128、256を新しい出力へ実行できる。こちらも元基準の代用にはしない。

96 / 768 / 6144 Tet10はrefine=1/2/3で準備する。初めに96・全体32・同一材料上限/接触上限で本体と独立を一組解く。空間水準と時間刻みを同時に変更しない。

```powershell
foreach ($h in 1,2,3) {
  $P="work/NEW-PC-C03-h$h"
  python "$Study/c03diag_pc.py" --user-pc --mode init --output $P --solver target/release/examples/c03_bulk_cause_probe.exe --refine $h
  python "$Study/c03diag_pc.py" --user-pc --mode native --output $P --folder native-g32 --steps 32 --material 0.0001 --contact 0.000025 --timeout 86400
  if ($LASTEXITCODE -ne 0) { throw 'native incomplete; do not advance automatically' }
  python "$Study/c03diag_pc.py" --user-pc --mode independent --output $P --folder independent-g32 --steps 32 --material 0.0001 --contact 0.000025 --timeout 86400
  if ($LASTEXITCODE -ne 0) { throw 'reference incomplete/budget rejected; keep failure' }
  python "$Study/c03diag_pc.py" --user-pc --mode compare --output $P --folder native-g32 --other independent-g32
}
```

96ではRAM64 GB・1経路4〜24時間、768はRAM128 GB以上・1日〜数日、6144はRAM256 GB以上・数日以上を計画上の予算とする（全て未実測、完走保証なし）。dense参照は予算拒否・時間切れになり得る。15000DOFを超えたら現行参照は拒否し、その規模で独立全体求解済みとはしない。PC専用コードを置いたことをPC実行済みに数えない。

`native-g32--independent-g32-comparison.json`、`schedule.json`、`process.json`、全場 `fields.jsonl` と全実行時間を見る。同一メッシュの本体対独立を最初に評価し、単純なGP配列引き算で異なるメッシュを比較しない。異なる空間水準の場比較にはC-02既存の親要素/界面対応と演算子を用いる。最細条件だけを真値にせず、同じ規模で時間と局所分割を別々に調べる。外挿・Richardson補正・新閾値を受入へ使わない。

## 封印確認と全件ソフトウェア検証

```powershell
python "$Study/c03diag_evidence.py" --mode verify --source $Study --output work/NEW-C03-archive-check
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release
cargo test --locked --release
foreach ($g in 'python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence','python-c01-switch','python-c01-path','python-c01-cause') {
  python tools/run_ci_python_evidence_shards.py --group $g --output "work/NEW-C03-validation/$g"
  if ($LASTEXITCODE -ne 0) { throw 'full test group failed' }
}
```

全Python発見IDの和集合と実行IDが重複・欠落なく一致し、skip/expected failure/unexpected success/未実行を成功へ加えない。重い求解は通常CIで再実行しない。新しい軽い試験は通常core群に属し、C-01専用のアーカイブ除外群へ割り当てない。文書検査はテキスト読取りの改行正規化を用い、実CRLFコピーでも検査する。ノルムは固定順C-orderの `math.fsum` で計算し、Windows集計順差を許容値で吸収しない。

新規JSONは固定コミットへ保存してから既存 `register_sources.py --commit FULL_SHA --id UNIQUE_ID PATH` でP6補助台帳に登録する。旧台帳原本・既存登録を再生成しない。最終Linux件数は `validation.json`。push後のWindows CIは未確定。
