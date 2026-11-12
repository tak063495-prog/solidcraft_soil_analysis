# B-01 再現手順（Windows PC / Linux）

依頼フォルダー日付はラベル、実施日は2026-11-12（日本時間）。PRはマージしない。
リポジトリルートで実行し、出力は毎回未存在名にする。途中失敗も保持する。
時間切れ、キャンセル、未収束、出力欠落、未実行を成功に数えない。

## 環境・全件試験

```powershell
git fetch origin codex/plastic-spatial-convergence
git worktree add --detach D:\Work\b01-repro origin/codex/plastic-spatial-convergence
Set-Location D:\Work\b01-repro
git rev-parse HEAD
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
python -m venv .venv
.venv\Scripts\python -m pip install -r requirements-reference.lock.txt
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:CARGO_BUILD_JOBS='2'
$env:RUST_TEST_THREADS='2'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
$B='interim/2026-12-03-tri6-curved-b01'
$P='interim/2026-12-01-tri6-offset-composite-b02'
.venv\Scripts\python "$P/validate.py" --part rust --output runs/NEW-b01-rust
.venv\Scripts\python "$P/validate.py" --part python --output runs/NEW-b01-python
```

20〜120分、RAM 4〜16 GB（未実測目安）。完全checkoutが前提で、旧分割証拠は
元Git blobのまま使う。Rustはfmt、全target clippy -D warnings、locked release build、
全target試験、doc試験。Pythonはtests/toolsの全IDを発見して重複/欠落なく4群に分ける。
`validation.json`と各群ログの件数/終了コードを確認する。
旧NC mean/spatial各4件がsetUpClass厳密検査で止まった場合は未実行8・準備エラー2を
そのまま記録する。skip、許容差化、参照/封印書換え、試験除外で通さない。
過去の未実行・失敗記録を保持し、全件成功しなければvalidation.passed=false。

## 実CLIの受理/拒否（小型）

```powershell
.venv\Scripts\python "$B/probe.py" runs/NEW-b01-probes --cli target/release/soil-fem-3d.exe
```

数分目安。`scope.json`と各入力・PID/終了コード・診断・Run有無を保存する。
面内/面外中間節点曲がり、真円筒上へ中間節点を戻した曲面パッチ、曲面+正負e、
平面/法線/director偽装、範囲外e、共有IDの法線/e不一致は根拠どおり拒否すること。
`all_expected_contracts_passed`を確認する。通常の拒否だけでなく型検査の拒否を区別する。
旧mindlin/cst/既定CSTのeなし折れ目、平面offsetは完成Runを生成する回帰対照。
曲面を直線中点へスナップしたり未知フィールドを削除して拒否試験を成功へ変えない。

## 小〜中規模の新規求解・全場比較・別PID再開

```powershell
.venv\Scripts\python "$P/build_observer.py" work/NEW-b01-observer-build
.venv\Scripts\python "$B/study_b01.py" runs/NEW-b01-fresh --binary target/release/b02-observer.exe
.venv\Scripts\python "$B/curved_comparison.py" runs/NEW-b01-fresh runs/NEW-b01-curve
.venv\Scripts\python "$B/rigid_modes.py" runs/NEW-b01-fresh runs/NEW-b01-rigid.json
.venv\Scripts\python "$B/checkpoint_mutations.py" runs/NEW-b01-fresh runs/NEW-b01-mutations
```

計5経路：端部地盤付着の円筒多面体2/4/8周分割、古典屋根2/4周・軸分割。
最大32枚Tri6、81構造節点、地盤45節点/12 Tet10。各工程4刻み。
10〜60分、RAM 4〜8 GB、空き3 GB以上（未実測目安）。
`native/accepted.jsonl`、`independent.jsonl`は双方とも入力からの**新規求解**。
本体の解/行列/設置参照/受理経路を独立解生成へ渡さない。
同じ実行版の観測なしRunとも照合する。保存JSON再監査を新規求解に数えない。

各`comparison.json`の`all_point_component_comparisons`で全受理ID・全成分を確認する。
土全GP、全構造節点の並進/回転/設置参照、全活動板6回復点N/Mp/Ms/Qとエネルギー、
支持反力/モーメント、荷重仕事を比較する。未付着節点へ存在しない土補間を割り当てない。
付着節点は土との並進一致、各要素内点/辺点は付着面補間運動学を検査する。
既存ゲートはSPECどおりで、反力/仕事/曲面/一般細分化の新しい許容値は設定しない。
既存場ゲートが通っても、包括受入null、B-01全体passed=falseを維持する。

`restarts.json`の全行で別PID、後続全状態一致、Run全体一致が必要。
除外するのはelapsed_secondsだけで、符号付きゼロ、全試行/監査/設置参照/入力SHAを含む。
WindowsではそのPC自身の実行版からcheckpointを新規発行する。
Linux checkpointを別SHAのWindows実行版へ付け替えない。
`mutations.json`では外側hashを再計算しても受理経路改変が拒否されること。
checksum未更新と未知法線の型拒否は別区分として保持する。

`curved-comparison.json`はPR #20をバイト保持したQ4 director参照の**新規求解**と比較する。
古典屋根の同一物理頂点全並進差、符号付き自由辺中央uz、文献値-.3024との差を記録する。
Q4とTri6の形状/director/drilling/離散化差を含み、同方式の実装一致とは解釈しない。
文献値は厳密解ではなく、曲面変位許容値・最低次数も未定義。
円筒・半球の本体系列を実施したとは主張しない。

## より細かい系列（PCのみ・クラウド未実行）

まず古典屋根8分割と、端部付着の16/32分割を別出力へ実行する。

```powershell
.venv\Scripts\python "$B/study_b01.py" runs/NEW-b01-finer --binary target/release/b02-observer.exe --attached-levels 16 32 --roof-levels 8
.venv\Scripts\python "$B/curved_comparison.py" runs/NEW-b01-finer runs/NEW-b01-finer-curve
```

30〜180分、RAM 8〜16 GB、空き5 GB以上（未実測）。屋根8は128 Tri6/289節点。
`--steps 8`、`--steps 16`も別出力で実行し、共通progressの全場・瞬間仕事を比較する。
差の再増加、負次数、既存ゲート超過は元のまま保存する。
粗い2水準だけで漸近収束を宣言しない。

屋根16は既存の構造数/凝縮自由度予算を超え得るため、自動的に上限を広げない。
本体の既存入力から準備して容量拒否も記録する：

```powershell
.venv\Scripts\python "$B/b01_cases.py" runs/NEW-b01-roof16 --kind roof --n 16 --steps 4
$L='runs/NEW-b01-roof16'
& target/release/soil-fem-3d.exe solve-stages "$L/mesh.json" "$L/ground.json" "$L/stages.json" "$L/run.json"
```

RAM 16〜32 GB、数時間の可能性（未実測）。拒否/未完了を受入にしない。
大規模の独立全場/復旧がなければ本体だけの結果をB-01合格へ転用しない。
draft Q4だけの細かい3形状は、同梱した元プログラムで別に実行できる：

```powershell
.venv\Scripts\python "$B/draft-reference/reference.py" --levels 32 64 128 --allow-large --out runs/NEW-b01-draft-large
```

5〜60分、RAM 8〜32 GB（未実測）。CSV、results.json、符号付き差・非単調性・残差を保存。
これは本体曲面Tri6の受入でも、元draftの結果の置換でもない。

## 表示・GUI/GPU・実CAD（PC手順のみ）

```powershell
New-Item -ItemType Directory runs/NEW-b01-view | Out-Null
.venv\Scripts\python -c "import sys,json;from pathlib import Path;sys.path.insert(0,'tools');import plate_io;plate_io.export(json.loads(Path('runs/NEW-b01-fresh/attached-n4-s4/native/run.json').read_text()),Path('runs/NEW-b01-view'))"
$SolidCraftCli='D:\Apps\SolidCraft\solidcraft.exe'
.venv\Scripts\python tools/prepare_plate_wall.py runs/NEW-b01-CAD-baseline --cli $SolidCraftCli --mesh-sizes-m 4 3 --example plate-wall
```

表示10〜30分、CAD準備5〜60分（未実測）。VTU/CSV/PVDの付着面/中立面、法線、
e符号、各工程の設置参照・回転を確認し、法線平均や点mergeで多面体の折れ目を消さない。
CADコマンドは既存平板の基準だけで、曲面Tri6/offset Native準備は未対応のまま。
実行版はユーザーPCの実配置を指定する。GPUを選択する実装済みソルバーオプションはない。
実CAD曲面/director/荷重仕事/全場/復旧の独立検証、GUI/GPUの計測は今回未実行。

## 分割証拠の照合（新規求解ではない）

```powershell
.venv\Scripts\python "$B/package_evidence.py" --verify "$B/evidence-manifest.json" --scratch work/NEW-b01-verify
```

全part、連結ZIP、全memberのsize/SHAが一致すること。`verification.json`に件数を保存する。
source snapshotは求解開始時と後処理追加後を区別する。旧snapshot/hashは書き換えない。
Linuxは`.venv/bin/python`、拡張子なしの実行版、exportで読み替える。
push後のWindows CIは当該SHAの結果だけを確認し、旧成功やLinux結果を転用しない。
