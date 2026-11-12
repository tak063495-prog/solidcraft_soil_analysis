# B-03 再現手順（Windows PC / Linux）

保存先日付は依頼ラベル、実施日は2026-11-12（日本時間）。PR #13および参照PRをマージしない。
リポジトリルートで実行し、出力先は毎回未存在名とする。終了コードとログを保存する。
時間切れ・中断・未収束・結果欠落を成功にしない。

## 環境と全件試験

```powershell
git fetch origin codex/plastic-spatial-convergence
git worktree add --detach D:\Work\b03-repro origin/codex/plastic-spatial-convergence
Set-Location D:\Work\b03-repro
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
$B='interim/2026-12-02-tri6-offset-kink-b03'
$P='interim/2026-12-01-tri6-offset-composite-b02'
.venv\Scripts\python "$P/validate.py" --part rust --output runs/NEW-b03-rust
.venv\Scripts\python "$P/validate.py" --part python --output runs/NEW-b03-python
.venv\Scripts\python "$B/report_validation.py" runs/NEW-b03-rust runs/NEW-b03-python runs/NEW-b03-validation.json
```

Rust検査は `cargo fmt --all --check`、`cargo clippy --locked --all-targets -- -D warnings`、
locked release build、`cargo test --locked --release --all-targets`、doc試験をすべて実行する。
Pythonはtestsとtoolsの全発見集合を4群へ重複/欠落なく分け、実行IDも照合する。
PC想定20〜120分（未実測）。完全checkoutが前提。疎checkoutの場合、旧ZIPと分割tarを
元Git blobのまま復元する。省略したまま試験を除外しない。

正常な実行件数は本フォルダーの `validation.json` と照合する。
旧NC mean/spatial各4件のsetUpClass厳密検査が失敗した場合は、準備エラーと未実行8件を
そのまま残す。Linuxの比較を許容差化せず、元のOS分岐・参照・封印を変えない。
過去の未実行記録は保持する。全件成功しなければvalidation.passedもfalseである。

## 実CLIの受理/拒否

```powershell
.venv\Scripts\python "$B/probe.py" runs/NEW-b03-probe --cli target/release/soil-fem-3d.exe
.venv\Scripts\python "$B/supplemental.py" runs/NEW-b03-supplemental --cli target/release/soil-fem-3d.exe
.venv\Scripts\python "$B/geometry_audit.py" runs/NEW-b03-probe runs/NEW-b03-probe-geometry.json
```

数分程度の小型入力。`scope.json` と各入力・CLIログ・PID/終了コードを保存する。
非ゼロ腕の15/45/90/135/179度実折れ目、共有IDのe段差/法線のみ反転/eのみ反転、
±.051 mの範囲外eが現行の診断で拒否されること。
共平面の法線/e同時反転、e=0の90度折れ目、別IDで回転を分けた面、
異なる高さの独立平板は実CLIが完走する入力である。それを溶接折れ目の成功にしない。

土メッシュも変えた曲がり2入力はCAD体積整合検査が先に拒否するため、`probe.py`の
当初の期待診断との不一致falseを保持する。`supplemental.py`は未付着の板節点だけを
曲げ、板の直線/平面ガードによる拒否を独立に確認する。補助実行の成功で当初JSONを
書き換えない。旧mindlin/cst/既定CSTのe=0折れ目は回帰対照であり、新方式の一般受入ではない。

PR #25の294条件を新たに点幾何として評価する。`accepted=null`、未定義の幾何許容値を保持する。
点幾何、候補拘束の核、CLIの受理、独立力学求解の4つを混同しない。

## 小〜中規模の新規求解と別PID再開

```powershell
.venv\Scripts\python "$P/build_observer.py" work/NEW-b03-observer-build
.venv\Scripts\python "$B/study_b03.py" runs/NEW-b03-formal --binary target/release/b02-observer.exe
.venv\Scripts\python "$B/convergence_b03.py" runs/NEW-b03-formal runs/NEW-b03-convergence.json
.venv\Scripts\python "$B/orientation_audit.py" runs/NEW-b03-orientation-audit.json
.venv\Scripts\python "$B/rigid_modes.py" runs/NEW-b03-formal runs/NEW-b03-rigid-modes.json
.venv\Scripts\python "$B/checkpoint_mutations.py" runs/NEW-b03-formal runs/NEW-b03-mutations
.venv\Scripts\python "$B/verify_contract.py" runs/NEW-b03-formal runs/NEW-b03-probe/scope.json runs/NEW-b03-supplemental/supplemental.json runs/NEW-b03-mutations/mutations.json runs/NEW-b03-contract.json
.venv\Scripts\python "$B/collect_b03.py" runs/NEW-b03-formal runs/NEW-b03-probe/scope.json runs/NEW-b03-supplemental/supplemental.json runs/NEW-b03-mutations/mutations.json runs/NEW-b03-summary
```

21経路。平面対照1と、5種類×（12 Tet10・刻み4、96 Tet10・刻み4、12 Tet10・刻み8/16）。
45/225土節点、12/96 Tet10、乾燥MC弾性枝、厚さ.1 m。最大2面の平板群である。
各経路は観測付き本体、**入力だけの独立Python新規求解**、同じ実行版の観測なし本体。
`native/accepted.jsonl` と `independent.jsonl` は新たな全受理点の解で、保存JSON再監査ではない。
B-02の公開API観測器と独立実装を変更せず使用する。
PC想定30〜120分、RAM 4〜16 GB、空き5 GB以上（未実測）。

`comparison.json` は土全GP、全土/構造節点、全活動板の6回復点・全エネルギー・設置参照・
活動・反力・荷重仕事を照合する。既存ゲートはSPECのまま：再組立1e-10、非線形1e-8、
真の線形1e-10、土変位/応力相対1e-7、構造並進/回転絶対1e-9、
板生値絶対1e-8/相対1e-7、付着拘束1e-10 m。
反力合計・仕事・一般折れ目/細分化の許容値を新設せず、包括受入nullを保持する。

`convergence_b03.py` は共通土節点・土GP位置・板物理点、全体基底のN/M/Q、
粗い二次仮想変位の転置による反力制限を用いる。
法線/e同時反転ではNは不変、Mp/Ms/Qは向き付き断面量なので、元の粗い面法線へ
sign(n_f・n_c)を掛けて比較する。各点の生法線/e・面ID・変換係数を保存する。
`orientation-audit.json`の独立座標変換診断は求解でも新しい許容値でもない。
空間は2水準の差だけを報告し、次数を作らない。荷重は4→8→16の共通20点で比較する。
`work-audit.json`は瞬間F・dと、各履歴自身の設置参照からの受理端点台形積分を分ける。
施工保持荷重を含む離散観測であり、厳密な連続経路仕事の受入ではない。

`restarts.json`の全行で別PID、`run_exact_except_elapsed=true`、`all_later_states_exact=true`を確認する。
除外はelapsed_secondsだけ。符号付きゼロ、全監査・全試行・入力/実行版SHAを含む。
生成元と同じ固定コピー実行版を使い、WindowsではPC自身のcheckpointを新たに発行する。
LinuxのcheckpointをWindows版の別実行版へ直接付け替えない。
改変30件（内容再hash24、checksum未更新3、型拒否3）は全てRun未生成で拒否が必要。

## 追加細分化：PCだけ（クラウドで未実行）

768 Tet10を含む3水準を、選んだ幾何だけで新規求解する：

```powershell
.venv\Scripts\python "$B/study_b03.py" runs/NEW-b03-h2 --binary target/release/b02-observer.exe --kinds kink-zero --level-max 2
.venv\Scripts\python "$B/convergence_b03.py" runs/NEW-b03-h2 runs/NEW-b03-h2-convergence.json
```

12→96→768 Tet10、最大1377土節点、土自由度4131で既存独立直接系の6500上限内。
RAM 8〜24 GB、空き10 GB以上、1〜6時間（未実測）。受理点を間引かず全場比較する。
`--kinds step-split` / `kink-split` / `normal-sign-pair` / `height-step`で別々の出力へ実行できる。
差が増加、既存ゲート超過、未完了なら未達として保存し、参照・許容値を書き換えない。
B-02の768 Tet10再組立超過とせん断/仕事の負次数は過去記録のまま保持する。

6144 Tet10は独立直接参照の上限外であり、本体だけの追加診断に限る：

```powershell
.venv\Scripts\python "$B/prepare_pc.py" runs/NEW-b03-large --kind kink-zero --level 3 --steps 16
$L='runs/NEW-b03-large'
& target/release/soil-fem-3d.exe solve-stages "$L/mesh.json" "$L/ground.json" "$L/stages.json" "$L/run.json"
& target/release/soil-fem-3d.exe checkpoint-stages "$L/mesh.json" "$L/ground.json" "$L/stages.json" 40 "$L/checkpoint-40.json"
& target/release/soil-fem-3d.exe resume-stages "$L/mesh.json" "$L/ground.json" "$L/stages.json" "$L/checkpoint-40.json" "$L/resumed.json"
.venv\Scripts\python -c "import sys;sys.path.insert(0,'$P');from cases import read;from study import exact;assert exact(read('$L/run.json'),read('$L/resumed.json'))"
```

RAM 16〜32 GB、空き20 GB以上、各1〜8時間（未実測）。既存容量拒否が出たらそのまま残す。
`scope.json`に寸法・参照上限を出す。容量・独立上限を広げたり、拘束を加えて通したりしない。
独立全場比較がない大規模結果をB-03の受入へ転用しない。

## 表示・実CAD・GUI/GPU（PC手順のみ）

完成Runを既存VTU/CSV出力へ渡す：

```powershell
New-Item -ItemType Directory runs/NEW-b03-view | Out-Null
.venv\Scripts\python -c "import sys,json;from pathlib import Path;sys.path.insert(0,'tools');import plate_io;plate_io.export(json.loads(Path('runs/NEW-b03-formal/kink-split-h0-s4/native/run.json').read_text()),Path('runs/NEW-b03-view'))"
```

`plates.pvd`、`plate-samples.pvd`、工程ごとのVTU/CSV、`plate-fields.json`を外部Viewerで確認する
（10〜30分目安）。付着面triangle6位置/並進、中立面reference位置/total並進、生値位置、
符号付きe、局所基底、設置基準、撤去後の空セルを確認する。
別IDの同位置点は回転が別であり、表示時のmerge/weld/法線平均で消してはならない。
これはNative GUI/GPUの受入ではない。Native準備は新方式を拒否したままであり、GPUを
選ぶ実装済みソルバーオプションもない。存在しないオプションや成功結果を作らない。

実CADの**オフセットなし基準**は既存手順で別ジョブとして用意する：

```powershell
$SolidCraftCli='D:\Apps\SolidCraft\solidcraft.exe'
.venv\Scripts\python tools/prepare_plate_wall.py runs/NEW-b03-CAD-baseline --cli $SolidCraftCli --mesh-sizes-m 4 3 --example plate-wall
```

PCの実配置に合わせて実行ファイルを指定する。準備5〜60分、計算はCAD/格子依存。
元CAD・捕捉・メッシュ・向き・境界/材料/荷重のSHAを保存する。
現行コンパクト/Native経路から非ゼロオフセット折れ目への変換と表示は未対応であり、
この準備コマンドを新方式の成功に数えない。
実CADの付着面/中立面/結合条件、剛体6モード・余剰零モード、全場・仕事・設置参照・
別PID再開・GUI/GPU全工程の独立検証が揃うまで未達とする。

Linuxでは`.venv/bin/python`と拡張子なしの実行版、`export NAME=value`に読み替える。
新先端Windows CIは、そのコミットの実結果だけを別途確認する。


## 封印証拠の再読込（新規求解ではない）

```powershell
.venv\Scripts\python "$B/package_evidence.py" --verify "$B/evidence-manifest.json" --scratch work/NEW-b03-evidence-verify
```

数分、空きZIPの2倍程度。全part・ZIP全体・全memberのサイズ/SHAが一致し、
`scratch/verification.json`に検証数が出ること。保存JSON/ZIP再監査を新規求解に数えない。
`source-snapshot-manifest.json`は実行前後の743ファイルと、その後の補助コード/文書を区別する。
旧証拠のsnapshot/hashを書き換えない。新規求解時の本体/入力だけの独立の証拠は
`fresh/*/native/accepted.jsonl`と`fresh/*/independent.jsonl`に保持する。
