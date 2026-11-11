# 再現・追加検証（範囲拡張(c)）

依頼の保存先ラベルは2026-11-15。実施日は2026-11-11（日本時間）。PR #13、`codex/plastic-spatial-convergence`、起点 `bc6eb1c39749559db91e20d2a4cba2d20679649c`。同じブランチを使い、PRをマージしない。出力先は必ず未存在のパスにする。

```sh
git fetch origin codex/plastic-spatial-convergence
git switch codex/plastic-spatial-convergence
git pull --ff-only origin codex/plastic-spatial-convergence
git merge-base --is-ancestor bc6eb1c39749559db91e20d2a4cba2d20679649c HEAD
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
python -m pip install -r requirements-reference.lock.txt
```

Linux/WSLは `OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 CARGO_BUILD_JOBS=2 RUST_TEST_THREADS=2 PYTHONUTF8=1 PYTHONIOENCODING=utf-8` をexportする。WindowsはPowerShellの `$env:OPENBLAS_NUM_THREADS='1'` 等を使う。以下の `python` はロック依存を入れた環境を指す。

sparse checkoutを使用している場合は、新しい証拠partとコードを取得してから再監査する。

```sh
git sparse-checkout add '/interim/2026-11-15-wet-joint-convergence/' '/interim/2026-11-15-wet-joint-convergence/*.part*' '/interim/2026-11-15-wet-joint-convergence/retained-failures/*.part*'
```

## 保存証拠を再計測する（新規求解なし）

```sh
python interim/2026-11-15-wet-joint-convergence/evidence.py --output work/NEW-wet-c-replay
```

`archive-sha256.json` の全part・ZIP本体・全エントリーSHAを確認してから展開する。全受理場から比較・Richardson観測次数・切替イベント・イベント整列診断を再計算し、保存JSONとの完全一致を要求する。`software_replay_passed=true` と `numerical_acceptance=false` を区別する。欠落・改変・不正な合格フラグは停止条件である。これは保存場の再計測であり、新規求解や別プロセスcheckpoint再開ではない。封印ソースを隔離したプロセスで再計測する場合は、上の展開先に対して次を実行する。

```sh
python -I -X utf8 work/NEW-wet-c-replay/sources/interim/2026-11-15-wet-joint-convergence/measure.py work/NEW-wet-c-replay/paths work/NEW-wet-c-sealed-measure
```

その `measurements.json`・`switch-events.json` を主記録とJSONオブジェクトで完全一致比較する。キーの表示順は数値差ではない。不一致を検出したら保存し、許容値や参照値を変更して通さない。

## 小型系列を入力から新規求解する

```sh
cargo build --locked --release
cargo build --locked --release --example wet_joint_history_probe
python interim/2026-11-15-wet-joint-convergence/study.py --output work/NEW-wet-c-solve
python interim/2026-11-15-wet-joint-convergence/measure.py work/NEW-wet-c-solve work/NEW-wet-c-measure
```

Windowsの本体保存プローブは `--solver target/release/examples/wet_joint_history_probe.exe` を `study.py` に指定する。物理モデルを変更するプローブではなく、既存 `solve_until` → `Checkpoint::seal/resume(Some(1))` を繰り返し、全受理点の確定状態を保存する。最後に同じ入力の中断なし求解と、全Runの非時間値を完全一致で照合する。完了フラグだけが変わる再開は受理点に数えない。同一プロセスのAPI再開であり、実ファイルからの別プロセス完了再開は未認定。

独立計算は既存 `reference_history.py` の入力のみの弾性MC＋Coulomb解を観測用に複製した `reference_probe.py` を使う。固定接触分割の既存演算は保持する。増分上限対照では点ごとの `ceil(norm(Δjump)/maximum)` を使い、Newton中は分割数を固定し、収束候補で不足する場合に増やして同じ旧状態から再評価する。本体の予測子・分割倍増戦略とは同一アルゴリズムではなく、上限と全体刻み・メッシュが同じ対照である。分割数は独立経路の各受理行に保存する。本体場・本体受理係数の注入はしない。

空間系列は既存12 Tet10→既存 `refine_tet10.refine` による96 Tet10。形状・材料・水位・面集合・上下の拘束変位・7段階を保持する。固定した全体刻みでhだけを変える対照と、h/2かつ全体刻み1/2の対照を分離する。異なる接触GPを同じ点とは扱わず、空間比較の新指標は親面ごとの積分平均と元の共通節点に限定する。面内分布を隠す平均値から点ごとの全履歴収束を認定しない。

未収束・時間切れ・キャンセル・欠落では `failure.json` と途中場を残す。終了コードや7段階完了、整合監査の通過を数値受入へ読み替えない。新規反復が失敗した場合は独立・本体のどちらか、進行度、残差、分割数、終了理由を保存する。

## 256増分以上・768要素以上はユーザーPCだけ

次はクラウドでは実行していない。空間3水準と全体刻み3水準を交差させる。資源予算を確保し、停止時も結果を消さない。

```sh
cargo build --locked --release --example wet_joint_history_probe
python interim/2026-11-15-wet-joint-convergence/user_pc.py \
  --solver target/release/examples/wet_joint_history_probe \
  --output work/NEW-PC-wet-c-h012-g64-128-256 \
  --levels 0 1 2 --steps 64 128 256 --maximum 0.00005 --native-timeout-seconds 14400
```

PowerShellでは実行ファイル名に `.exe` を付け、継続行の `\` をバッククォートへ変えるか1行で実行する。この明示PCドライバだけが本体プローブの `--user-pc` を付ける。独立密行列の上限は15,000 DOFであり、越える問題はこの参照器で実行しない。PC用の本体求解予算はこの例で1経路14,400秒であり、`--native-timeout-seconds` で利用者が変更できる。予算停止は不合格として途中場を残す。全履歴比較は配列をメモリに展開するため、密行列だけでなく保存場のRAMとディスクも確保する。レベル2は768 Tet10で、体積GP・接触GP・全履歴の保存が大きくなる。独立参照も全前史を新規に解く。これで一般バルク塑性の独立解になるわけではない。

続いて同じメッシュ・全体刻みで `--maximum 0.000025`、必要なら `0.0000125` の別出力を作る。空間だけ、全体刻みだけ、接触上限だけ、同時細分化の比較を混ぜず、全共通座標・18場の既存L2/絶対差、元の7端点ゲート、新指標を併記する。`measure.compare`、`measure.spatial`、`measure.richardson` は任意の保存経路に使える。Richardsonは同じ座標の3水準だけを使い、差の方向余弦と切替状態も確認する。負次数・方向不一致・履歴非単調性を外挿の合格根拠にしない。

独立4096/64・本体512などの過去の記録はその版のSHA・積分方式・接触分割・比較座標と対応付ける。今回の小さい系列へ過去の約5.08%差やユーザーPCの128/256診断L2 0.01883を合格値として注入しない。旧参照値は変更しない。

## Linux/Windows全件チェック

```sh
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release
cargo test --locked --release
python tools/run_ci_python_tests.py --group python-core-0 --output work/NEW-check/core0
python tools/run_ci_python_tests.py --group python-core-1 --output work/NEW-check/core1
python tools/run_ci_python_tests.py --group python-mean-evidence --output work/NEW-check/mean
python tools/run_ci_python_tests.py --group python-spatial-evidence --output work/NEW-check/spatial
```

現時点の全テストIDが重複なく実行されたことと、fail/error/skipped/ignored/filtered-out 0を確認する。保存証拠の全件監査には既存mean約5.49 GB・spatial約0.92 GB等が必要。sparse checkoutの場合は必要なZIPとpartを取得し、取得完了前に試験しない。旧証拠の再監査は新規大規模求解ではない。

追加系列は最大128全体増分・96 Tet10に限定した。全件試験の既存小ケースや保存証拠の再監査とは区別する。

## この検証で完了にしない残件

湿潤JOINTの点ごとの空間・全履歴塑性すべり収束、一般バルク塑性と構造を含む独立3D全場、同じ物理経路の完了までの別プロセス再開、低プレロードNC全撤去は未達のまま。NCチューニングは行っていない。実CAD、実断面SRM、Native GUI/GPU、別PC/GPU、Griffithsや既存MCC大規模系列の再求解・完成配布はユーザーPCの別作業であり、今回の証拠では認定しない。既存の入力・出典・許容値と各専用手順を使う。
