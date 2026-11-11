# Python湿潤JOINT参照器の再現（Windows / Linux）

保存先の2026-11-18は依頼のラベル。実施日は2026-11-11（日本時間）。
起点は `5f564621c7396c2a605edbfe5451cd800cf318d3`、作業ブランチは
`codex/py-wetjoint-reference`。既存ファイルを編集せず、出力には未存在のフォルダーを使う。
この参照器は12 Tet10・54節点の既存湿潤小柱を入力から解く。
96要素までの軽量対照と、768要素以上のPC計算を区別する。

## 環境と全Python試験

リポジトリ直下のPowerShellで実行する。Python 3.13、ロック済みNumPy/SciPyを使う。
依存を追加・更新する必要はない。以下のコマンドはGit checkout後に実行する。

```powershell
git fetch origin codex/py-wetjoint-reference
git switch codex/py-wetjoint-reference
git merge-base --is-ancestor 5f564621c7396c2a605edbfe5451cd800cf318d3 HEAD
py -3.13 -m venv .venv-reference
.\.venv-reference\Scripts\Activate.ps1
python -m pip install -r requirements-reference.lock.txt
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONDONTWRITEBYTECODE='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
python -m unittest discover -s tools/py_reference/wet_joint/tests -v
```

新規試験は数秒〜30秒を想定。既存全Python試験にはRust CLI・診断exampleと保存証拠が必要。
既存証拠partの取得が完了してから実行する。sparse checkout利用時は各証拠フォルダーも含める。
既存のCI分類器をそのまま使い、4群すべての終了コードを確認する。

```powershell
rustup toolchain install 1.99.1 --profile minimal
cargo build --locked --release
cargo build --locked --release --examples
python tools/run_ci_python_tests.py --group python-core-0 --output work/NEW-python-core0
python tools/run_ci_python_tests.py --group python-core-1 --output work/NEW-python-core1
python tools/run_ci_python_tests.py --group python-mean-evidence --output work/NEW-python-mean
python tools/run_ci_python_tests.py --group python-spatial-evidence --output work/NEW-python-spatial
```

各群の `python-test-inventory.json` の `tests_run == assigned_total`、終了0、
failures/errors/skipped/expected_failures/unexpected_successes がすべて0であることを確認する。
全 `assigned` IDの和集合が `discovered_test_ids` と一致し、重複0であることも確認する。
既存全件は保存証拠約7 GBの読込・展開を含み、通常10〜60分程度の予算を確保する
（PC性能・ストレージ・証拠キャッシュで変動。Windowsで今回未実測）。
新しいpackage初期化ファイルにより、既存分類器も新規試験を自動発見する。CI設定は変更しない。
全件の内訳（既存＋新規）を報告し、新規の単独実行を全件数へ二重計上しない。

Linuxでは `.venv/bin/python` を使い、環境変数は `export OPENBLAS_NUM_THREADS=1` 等で設定する。
試験skip、タイムアウト、キャンセル、未完了は成功に含めない。

## 軽量の新規求解と本体比較

Pythonだけで計算できる。Rustの出力・既存参照ソルバは参照解の生成に使わない。

```powershell
python -m tools.py_reference.wet_joint.study --output work/NEW-py-wet-reference --steps 16 32 64 128 --levels 0 1
```

Linuxで約1〜2分。一般的なWindows PCでは2〜10分を計画値とする。
12→96 Tet10、各工程16/32/64/128増分、接触内8固定、別途64増分で接触内16/32、
32増分で残差1e-11の対照を解く。全7工程で履歴を最初から生成する。
主出力は `measurements.json`、`diagnostics.json`、各経路の `history.json` / `fields.npz`、
入力JSONとSHA指紋、失敗時にも残る `*-partial.jsonl`。

本体の各受理点を比較する場合は、既存の観測exampleをビルドして指定する。
観測exampleは既存checkpoint APIの呼出しであり、本体の修正ではない。

```powershell
cargo build --locked --release --example wet_joint_history_probe
python -m tools.py_reference.wet_joint.study --output work/NEW-py-wet-native --steps 16 32 64 128 --levels 0 1 --native-solver target/release/examples/wet_joint_history_probe.exe --native-steps 16 32 64 --native-timeout-seconds 900
```

Linuxでは `.exe` を除く。本体6経路と独立経路の計算、受理点JSONの保存に、
Windowsでは5〜30分、空きディスク5 GB程度を見込む（計画値）。
900秒は本体1経路あたりの予算であり、超過した経路を完了や合格にはしない。

本体は接触増分上限5e-6 m、Pythonは接触内8固定なので、局所の分割列は同一ではない。
両者の差は実装差と局所積分差を含む。本体の力1e-8・真の線形1e-10、参照の力1e-10は維持する。
本体が違う受理格子へカットバックした場合、同じ全体刻みの比較として通さず停止する。

## 比較スクリプトを単独で使う

正規化出力を使う。SHA改変、非有限値、工程・受理点欠落、境界係数不一致はエラー。
各対のL2と、3水準の全共通座標に固定した差ベクトルからの次数・方向余弦を出す。

```powershell
python -m tools.py_reference.wet_joint.compare --axis load --paths work/NEW-py-wet-native/h0-g16-l8 work/NEW-py-wet-native/h0-g32-l8 work/NEW-py-wet-native/h0-g64-l8 work/NEW-py-wet-native/h0-g128-l8 --output work/NEW-load-comparison.json
python -m tools.py_reference.wet_joint.compare --axis space --paths work/NEW-py-wet-native/h0-g32-l8 work/NEW-py-wet-native/h1-g32-l8 --parents work/NEW-py-wet-native/parents.json --output work/NEW-space-comparison.json
python -m tools.py_reference.wet_joint.compare --axis pair --paths work/NEW-py-wet-native/h0-g32-native work/NEW-py-wet-native/h0-g32-l8 --output work/NEW-native-pair.json
```

GPが異なる空間比較は元の2親面ごとの面積積分平均と、土塊所属を保った共通節点変位。
点ごとのすべりの空間収束の代用にはしない。2空間水準では次数は空欄・理由付きになる。
ゼロ分母で差が非ゼロの相対誤差も空欄・理由付きにし、小さい分母で合格を作らない。

既存観測exampleの別出力を変換する場合：

```powershell
python -m tools.py_reference.wet_joint.native --folder work/OTHER-native-raw --mesh work/NEW-py-wet-native/mesh-h0.json --ground work/NEW-py-wet-native/ground.json --case work/NEW-py-wet-native/stages.json --steps 32 --output work/NEW-other-normalized
```

`verification.json`・`run.json`・`layout.json`・`accepted.jsonl` が必要。
入力物性・荷重・拘束・元の残差基準、全受理列、接触GPの位置/法線/面積を確認する。
別プロセスの実checkpointファイル再開を認定する機能ではない。

## 大きい空間3水準・長い系列はWindows PCだけ

以下はクラウド未実行。少なくとも32 GB RAM、空きディスク50〜100 GBを推奨する計画値。
密行列・受理場保存を使うため、実測を見ながら予算を確保する。多重起動しない。
最大15,000 DOF・768 Tet10という参照器の上限を越える問題は対象外。
全体256以上または空間level2には `--user-pc` が必須である。

```powershell
python -m tools.py_reference.wet_joint.study --output work/NEW-PC-py-wet-h012 --steps 64 128 256 --levels 0 1 2 --native-solver target/release/examples/wet_joint_history_probe.exe --native-steps 64 128 256 --native-timeout-seconds 14400 --user-pc
python -m tools.py_reference.wet_joint.compare --axis space --paths work/NEW-PC-py-wet-h012/h0-g256-l8 work/NEW-PC-py-wet-h012/h1-g256-l8 work/NEW-PC-py-wet-h012/h2-g256-l8 --parents work/NEW-PC-py-wet-h012/parents.json --output work/NEW-PC-space-256.json
python -m tools.py_reference.wet_joint.compare --axis load --paths work/NEW-PC-py-wet-h012/h2-g64-l8 work/NEW-PC-py-wet-h012/h2-g128-l8 work/NEW-PC-py-wet-h012/h2-g256-l8 --output work/NEW-PC-load-h2.json
```

想定総時間2〜12時間（8〜16コア、単一BLASスレッド、SSDを想定した予算目安。未実測）。
14,400秒は本体1経路あたりの上限。開始前に軽量版を実行し、RAM・ディスク・実測速度を確認する。
Python密行列側には強制時間制限を設けていない。キャンセル時は `failure.json` と確定点列を保存する。
失敗した出力を消さず、新しいフォルダーで再試行する。

## 合否の読み方と次の診断

- ソフトウェア試験：指定件数の全件実行、終了0、失敗/error/skip 0。保存場・入力SHA一致が必要。
- 求解の完了：7工程の全受理点、参照力1e-10（対照のみ1e-11）、本体元の力1e-8・真の線形1e-10、弾性域・接触則の適用範囲を満たすこと。
- 数値精度：`passed=false` / `numerical_acceptance=false` を維持する。塑性すべり全履歴・点ごとの空間収束には既存の受入基準が未整備で、新しい許容値を計測値に合わせて作らない。負次数や非単調差は漸近収束の根拠にしない。
- 既存の端点ゲート（変位0.005、traction0.01）は全履歴のすべりゲートではない。今回はその値を緩和・上書きしない。

PCの3空間水準を固定荷重刻みごとに比較し、次に同じ空間水準の3荷重刻みを比較する。
開口直前の平衡経路を全体で細分化する対照を優先し、接触jumpだけの内分と区別する。
時変水位・Richards・非ゼロ流れの時間積分は別問題として独立水理参照を用意してから調べる。
NC-01、一般バルク塑性・部材統合、C-01/C-02の全受入、P6の完了には読み替えない。

## 保存証拠だけを再監査する（新規求解なし）

```powershell
python -m tools.py_reference.wet_joint.archive --index interim/2026-11-18-py-wetjoint-reference/archive-sha256.json --output work/NEW-py-wet-unpacked
python -m tools.py_reference.wet_joint.compare --axis load --paths work/NEW-py-wet-unpacked/data/study/h0-g16-l8 work/NEW-py-wet-unpacked/data/study/h0-g32-l8 work/NEW-py-wet-unpacked/data/study/h0-g64-l8 work/NEW-py-wet-unpacked/data/study/h0-g128-l8 --historical interim/2026-11-15-wet-joint-convergence/measurements.json --output work/NEW-py-wet-remeasured.json
```

全part・ZIP全体・各エントリーのSHA/サイズを照合してから、新規フォルダーへ展開する。
`verification.json` の `software_integrity_passed=true` は保存バイトの照合であり、数値精度合格ではない。
`data/study/measurements.json` の `load_series.h0` と再計測の対応欄をJSONオブジェクトで比較する。
`historical_comparison` は旧JSONの荷重増分の意味・`passed=false` を確認し、
旧L2と新L2の差をそのまま表示する。旧参照値を書き換えたり、新しい受入許容値を作ったりしない。
空間・本体対参照の再計測は、上記の比較CLIのパスを `data/study/` 配下へ置き換える。
保存場の再計測は、入力からの新規求解の代用ではない。
