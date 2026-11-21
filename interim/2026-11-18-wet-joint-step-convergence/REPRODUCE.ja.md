# C-01の再現とWindows PCへの引継ぎ

保存先ラベル2026-11-18、実施日2026-11-11（日本時間）。起点 `b069ec84da432f677c33fab5c984b6e680835bc4`、PR #13の同じブランチ。PRはマージしない。出力先は毎回未存在のパスを使う。旧証拠・参照値・失敗を上書きしない。

## 準備

PowerShell、Python 3.13、Rust 1.99.1を使う。以下の `python` はロック依存を入れた仮想環境を指す。

```powershell
git fetch origin codex/plastic-spatial-convergence
git switch codex/plastic-spatial-convergence
git pull --ff-only origin codex/plastic-spatial-convergence
git merge-base --is-ancestor b069ec84da432f677c33fab5c984b6e680835bc4 HEAD
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
py -3.13 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements-reference.lock.txt
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:CARGO_BUILD_JOBS='2'
$env:RUST_TEST_THREADS='2'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
cargo build --locked --release
cargo build --locked --release --example wet_joint_step_probe
```

各コマンドの終了コードが0でなければ次へ進まず、ログと途中場を保持する。Linux/WSLでは `.venv/bin/python`、実行ファイルの `.exe` を外し、上の環境変数を `export` する。このクラウドの既存Rustは `/workspace/.cargo/bin` と `/workspace/.rustup` にあった。一般PCでそのパスを流用する必要はない。

sparse checkoutの場合は必要な証拠を先に取得する。全件検査用の既存CIと同じ集合は次の通り。新しい証拠ZIP partも含める。取得完了前に検査を始めない。

```powershell
git sparse-checkout set --no-cone '/*' '!/interim/2026-11-10-nc-postfold-branches/evidence.tar.gz.part*'
```

この除外は既存CIの未使用保存系列だけであり、試験skip・整合検査削除ではない。NC系列の新規計算は今回行わない。

## 小ケースを入力から全履歴新規求解

```powershell
python interim/2026-11-18-wet-joint-step-convergence/study.py --mode init --solver target/release/examples/wet_joint_step_probe.exe --output work/NEW-c01-solve
python interim/2026-11-18-wet-joint-step-convergence/study.py --mode fixed --output work/NEW-c01-solve
python interim/2026-11-18-wet-joint-step-convergence/study.py --mode cap --output work/NEW-c01-solve
python interim/2026-11-18-wet-joint-step-convergence/study.py --mode native --output work/NEW-c01-solve
python interim/2026-11-18-wet-joint-step-convergence/study.py --mode open-fine --steps 64 --output work/NEW-c01-solve
python interim/2026-11-18-wet-joint-step-convergence/study.py --mode open-coarse --steps 256 --output work/NEW-c01-solve
python interim/2026-11-18-wet-joint-step-convergence/study.py --mode restart --steps 256 --output work/NEW-c01-solve
python interim/2026-11-18-wet-joint-step-convergence/measure_steps.py work/NEW-c01-solve work/NEW-c01-measure
```

16/32/64/128/256の各7工程を解く。`fixed` は独立接触内8固定、`cap` は独立接触上限5e−5 m、`native` は同じ上限の本体。全体刻みの差と局所接触分割の差を混同しない。介入はopenだけ256、またはopenだけ64にして**全体系を再求解**する。時刻整列ではない。独立系に本体解・本体受理係数を注入しない。

観測用Rust exampleは既存の `solve_observed` / `Checkpoint::seal` / `resume_observed` を呼ぶ。同一プロセスの観測再開に加え、N=256の各工程途中7個のcheckpointを実ファイルから**別プロセス**で復旧する。入力・ソルバーSHAとchecksumを既存APIで検証し、完了までの全後続確定状態と全Runを直接経路と照合する。時間 `elapsed_seconds` 以外はJSONオブジェクトの完全一致を要求する。さらに同じ実行版の `--whole` による中断なし求解と比較する。

出力は `gN/{fixed8,independent,native}/manifest.json` と可逆NPZ、独立の `*-result.json`、全場監査 `*-audit.json`、`native-raw/accepted.jsonl`・`progress.jsonl`・`run.json`・`checkpoint-NNNN.json`、`restart-checkpoint-NNNN/`、`restart-verification.json`。計測先は `measurements.json`、`switch-events.json`、`all-entity-differences.json`、`chord-defect.json` の4件。

クラウド実測は本体N=16/32/64/128/256が約10/16/29/62/112秒、独立8固定が約2/3/6/13/25秒、独立上限固定が約2/3/5/9/18秒（保存・監査を含み、並行作業の影響あり）。全7別プロセス再開・中断なし対照・計測を含む一式は通常10〜30分を見込む。CPU・ディスク・ウイルス対策に依存する目安であり、Windows実測時間や成功保証ではない。RAMは小ケースで8〜16 GB、出力用ディスクは5〜10 GBを目安に確保する。`--timeout 1800` は本体1プロセスの予算であり、停止を成功にしない。

合否は二層に分ける。

- ソフトウェア・整合確認：全7工程、入力の `i/N` 受理列、全GP・全節点・18場が欠けず有限値、元の力/線形基準、MC弾性枝、摩擦円・水圧・仕事・累積すべり・独立再組立、全7別プロセス再開と中断なし対照の完全一致。欠落・改変・未収束・時間切れは失敗として保存する。
- C-01の数値受入：全体刻みを細かくした全接触履歴・全節点・全材料GPの差の単調性・同一座標の観測次数・方向と独立系を確認する。今回の系列は再増加があり**未達**。7端点ゲート、接触内だけの収束、整合監査、再開一致、負次数や切替整列を受入に転用しない。既存にない全履歴許容値を測定値に合わせて新設しない。

## 封印証拠の再計測（新規求解ではない）

```powershell
python interim/2026-11-18-wet-joint-step-convergence/evidence_steps.py --output work/NEW-c01-replay
```

全part・ZIP・全エントリーSHAを照合してから展開し、4つの測定JSONを全保存場から再計算してオブジェクト完全一致を要求する。`software_replay_passed=true` と `numerical_acceptance=false` を区別する。実checkpoint求解は前節の `--mode restart` であり、保存JSONの再読込で代用しない。実行版はアーカイブに入れないため、新規求解では再ビルドした版のSHAを新規入力に結び付ける。封印された旧checkpointへ別バイナリのSHAを付け替えない。

## PCだけ：長い全体刻み系列

次は今回のクラウドで実行していない。既存PCドライバを使い、固定12要素の全体256/512/1024とする。

```powershell
cargo build --locked --release --example wet_joint_history_probe
python interim/2026-11-15-wet-joint-convergence/user_pc.py --solver target/release/examples/wet_joint_history_probe.exe --output work/NEW-PC-c01-g256-512-1024 --levels 0 --steps 256 512 1024 --maximum 0.00005 --native-timeout-seconds 14400
python interim/2026-11-18-wet-joint-step-convergence/pc_summary.py --input work/NEW-PC-c01-g256-512-1024 --output work/NEW-PC-c01-long-measure --steps 256 512 1024 --h 0 --maximum 0.00005
```

同じ入力から両実装の全前史・全18場を新規求解し、隣接全体刻みの共通受理点、全点の絶対差・L2、同一座標3水準の観測次数・方向余弦を保存する。`pc_summary.py` は外挿値を出さず、`passed=false` を維持する。出力は `h0-gN-m5e-05-cut/` の全場・本体Runと `NEW-PC-c01-long-measure/measurements.json`。さらに同じ系列を上限2.5e−5 mで別ディレクトリに実行し、局所上限依存を別に測る。必要なら2048へ延長する。小ケースでも観測checkpointごとの封印費用は増えるため、一系列30分〜2時間・RAM16 GB・ディスク10〜30 GBを見込む（未実測の資源計画）。本体予算は1経路4時間。予算終了・独立系停止・欠落は失敗であり、最後の成功点から新しい受入値を作らない。

## PCだけ：C-02の空間・時空間データ取得

C-02は今回着手せず、未達のまま。この既存コマンドは12/96/768要素を固定全体64/128/256で交差させるための生データ取得手順である。

```powershell
python interim/2026-11-15-wet-joint-convergence/user_pc.py --solver target/release/examples/wet_joint_history_probe.exe --output work/NEW-PC-c02-h012-g64-128-256 --levels 0 1 2 --steps 64 128 256 --maximum 0.00005 --native-timeout-seconds 14400
```

出力は `mesh-h0.json`〜`mesh-h2.json`、各 `hH-gN-m5e-05-cut/` の全受理場、本体/独立系の全履歴、最上位 `measurements.json`。同じ固定Nでhだけを比較し、同じ固定hでNだけを比較する。同時細分化は別欄とする。768要素は独立の密行列と全材料GP履歴の保存・展開が大きい。RAM64 GB以上、空きディスク100 GB、数時間〜1日を計画する（未実測、PC性能依存）。本体の1経路4時間予算を超えるなら、利用者の資源計画で明示変更し、停止を合格にしない。

**この既存ドライバの空間比較は親面平均と元の共通節点だけであり、C-02の点ごとの合格コマンドではない。** 生データ取得後も、同一物理点での全接触履歴と全材料場の評価・対応付け、空間3水準以上の収束、独立系との照合が未実装・未検証ならC-02は未達。2水準の平均差・Nとhの同時変化だけで合格にしない。新しい空間合格値は作らない。一般塑性の独立3D全履歴（C-03）、変動水位/Richards（C-04）、NC-01、P6、GUI/GPU・実CADは別タスクのまま。

## Linux/Windows全件検査

```powershell
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release
cargo test --locked --release
python tools/run_ci_python_tests.py --group python-core-0 --output work/NEW-c01-check/core0
python tools/run_ci_python_tests.py --group python-core-1 --output work/NEW-c01-check/core1
python tools/run_ci_python_tests.py --group python-mean-evidence --output work/NEW-c01-check/mean
python tools/run_ci_python_tests.py --group python-spatial-evidence --output work/NEW-c01-check/spatial
```

全Python発見IDと4群の実行IDの和集合が完全一致し、重複・欠落・fail/error/skipped/expected-failure/unexpected-successが0、Rustのfailed/ignored/filtered-outが0であることを数える。保存証拠の全件監査には旧mean約5.49 GB・spatial約0.92 GB等の取得・展開も必要。ディスクを確保し、取得未完了をskipして通さない。通常の全件検査は数十分〜1時間以上を想定するが、保存証拠の検査と大規模問題の新規求解は区別する。今回push後のWindows CIの成否は別途新先端のrunで確認する。
