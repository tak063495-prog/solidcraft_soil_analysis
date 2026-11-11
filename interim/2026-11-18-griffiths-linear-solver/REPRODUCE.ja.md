# Griffiths 線形ソルバ診断の再現

保存先ラベルは依頼どおり 2026-11-18。実施日は 2026-11-11（日本時間）。起点は `5f564621c7396c2a605edbfe5451cd800cf318d3`。本体の前処理は実装しない。元の力許容値 `1e-8`、線形許容値 `1e-10`、参照値、失敗記録を変更しない。

## Linux：追加試験と既存 Python 全件

リポジトリのルートで実行する。既存全件には本体 CLI を呼ぶ試験と過去の証拠 archive の完全再計測がある。後者は数 GB の展開を伴うが、大規模な新しい FEM 求解ではない。空き容量 30 GB 以上を目安とする。環境設定・依存取得を除き、release build は数分、追加試験は数秒、既存全件はおおむね 10～30 分（PC・ディスク依存）。

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-reference.lock.txt
export PATH="$PWD/.venv/bin:$PATH"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
export PYTHONUTF8=1 PYTHONIOENCODING=utf-8
export CARGO_BUILD_JOBS=2
cargo build --locked --release
mkdir -p runs/griffiths-linear-validation
python -m tools.py_reference.griffiths_solver.run_tests --report runs/griffiths-linear-validation/new-tests.json
python tools/run_ci_python_tests.py --group python-core-0 --output runs/griffiths-linear-validation/core0
python tools/run_ci_python_tests.py --group python-core-1 --output runs/griffiths-linear-validation/core1
python tools/run_ci_python_tests.py --group python-mean-evidence --output runs/griffiths-linear-validation/mean
python tools/run_ci_python_tests.py --group python-spatial-evidence --output runs/griffiths-linear-validation/spatial
```

各コマンドの終了コードを個別に確認する。各 `python-test-inventory.json` の実行 ID を合併し、重複・欠落なく発見全件と一致すること、fail/error/skip/expected-failure/unexpected-success がすべて 0 であることが合格条件。追加試験は既存 CI の探索外の新規フォルダーに置いているため、**追加 runner を必ず別に実行**する。CI 設定は変更しない。

## 小～中規模の行列比較（Linux / Windows 共通）

```bash
python -m tools.py_reference.griffiths_solver.compare --output runs/griffiths-linear-xy1-seed
python -m tools.py_reference.griffiths_solver.compare --system runs/griffiths-linear-xy1-seed/system --output runs/griffiths-linear-xy1-repeat
python -m tools.py_reference.griffiths_solver.refinement_probe --output runs/griffiths-small-elastic-h.json
```

所要時間はそれぞれ通常 10～60 秒程度、メモリーは数百 MB 以下を見込む。単位ブロックの追加 h 対照は 24/96/384 要素の弾性行列だけで、通常数秒、Griffiths h 系列ではない。新しい出力先を使う。既存出力は拒否する。既定は既存の固定ジオメトリから作る xy1-ny1、1,200 Tet10、全 y 変位固定、乾燥・零応力種の**弾性初期接線**。E=100000 kPa、ν=0.3、γ=20 kN/m³。初期種では強度低減 F=1 と元参照 F=1.35 の弾性接線は同じだが、後続塑性経路は異なる。

出力は `comparison.json` と `system/{matrix.npz,vectors.npz,metadata.json}`。元行列の真の線形残差と力正規化した線形残差、内部反復数、Krylov 系の残差履歴（`row_scaling` で生の系/行スケールを区別）、分解 nnz、時間、条件数推定、固有値標本を保存する。直接法の反復数 0 は GMRES 未使用を表す。全要求手法が元基準を満たす場合だけ exit 0。1 手法でも失敗すると行を保持して exit 1。ピボット失敗をシフトや別ソルバで隠さない。

保存済み任意の元接線 A と Newton RHS b も、`save_system` の形式で比較できる。`vectors.npz` は `rhs`（縮約自由度順）と `free_dofs`（元節点の `3*node+axis`）を含む。`metadata.json` には元の全外力自由度ノルム `external_force_norm` と、入力 SHA・進捗・反復・求解元 SHA を記録する。外力ノルムを省略した一般入力は `||b||` を使うので、Newton 残差 RHS では必ず指定する。`matrix_sha256` と `vectors_sha256` は実バイトの SHA256。読み込み時に照合する。

元の GP アルゴリズム接線が別途ある場合は `--mesh MESH --ground GROUND --case CASE --gp-tangents GP.npz` で同じ積分・拘束に再組立できる。NPZ の `tangent` は `[要素, GP, 6, 6]`、順序は元の element/local GP 順、成分は xx,yy,zz,xy,yz,zx、ひずみは engineering shear。非対称成分を保持する。この再組立の RHS は全重力であり、失敗した Newton の RHS そのものではない。**実際の失敗系を比較するには A と b を一緒に凍結する**。Run の確定応力だけから失敗反復のアルゴリズム接線を作ったことにはしない。

## Windows PC：xy2 本体の再現（クラウドでは実行しない）

PowerShell、Rust 1.99.1、既存 lock の Python 依存を使う。既存の失敗入力・出力は保持し、専用の別 checkout と新しい `runs` フォルダーで実行する。**失敗を再現する主経路は、失敗時と同じ mesh/ground/case 3 ファイルをそのまま渡すこと**。下記 `C:\soil-inputs\griffiths-xy2\...` を実際の保存場所へ置き換える。材料低減係数・ny・刻み・Newton35・GMRES3000・restart60 等を元入力/既定値と照合する。変更を混ぜない。

```powershell
git fetch origin codex/py-griffiths-linear-solver
git switch --detach origin/codex/py-griffiths-linear-solver
py -3.13 -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements-reference.lock.txt
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
cargo build --locked --release
if ($LASTEXITCODE -ne 0) { throw 'build failed' }
$mesh='C:\soil-inputs\griffiths-xy2\mesh.json'
$ground='C:\soil-inputs\griffiths-xy2\ground.json'
$case='C:\soil-inputs\griffiths-xy2\case.json'
$out='runs\griffiths-xy2-original-replay-01'
New-Item -ItemType Directory $out -ErrorAction Stop
Get-FileHash $mesh,$ground,$case,target\release\soil-fem-3d.exe | Format-List | Out-File "$out\sha256.txt"
$sw=[Diagnostics.Stopwatch]::StartNew()
& target\release\soil-fem-3d.exe solve-stages $mesh $ground $case "$out\result.json" 2>&1 | Tee-Object "$out\solve.log"
$solverExit=$LASTEXITCODE
$sw.Stop()
@{exit_code=$solverExit; elapsed_seconds=$sw.Elapsed.TotalSeconds} | ConvertTo-Json | Out-File "$out\process.json"
```

想定：4,800 要素・129,600 材料 GP、メモリー数 GB、数十分～数時間。これは未計測の予算目安で、完了を保証しない。元失敗の **0.44 付近の GMRES3000、geostatic 0.4717 の Newton35、相対力残差 1.109e-7 > 1e-8** を探し、受理/棄却、元状態保持、全 exit code を保存する。失敗再現は原因診断の観測であり、数値合格ではない。手動停止・時間切れ・キャンセルは未完了として記録する。

元入力が手元にない場合の別実験だけは、単一レベル準備器を使える。これは失敗入力との同一性を認定しない。

```powershell
.venv\Scripts\python.exe -m tools.py_reference.griffiths_solver.prepare_pc --xy 2 --output runs\griffiths-xy2-prepared-F1
# 元の公開 2D 固定低減材なら、別名で明示的に指定する：
.venv\Scripts\python.exe -m tools.py_reference.griffiths_solver.prepare_pc --xy 2 --reference-2d tests\data\griffiths_2d_reference.json --output runs\griffiths-xy2-prepared-F135
# 上の solve-stages の3入力を選んだ準備先の mesh.json / ground.json / case.json にする。
```

どちらも xy2-ny1 のみを作り、xy4 や長い h 系列を作らない。xy2 行列の Python 再組立は今回のガード（1,200 要素まで）で拒否する。PC で凍結した xy2 CSR 系は `--system` で比較可能。直接法の fill/メモリーをまず xy1 で確認し、xy2 の直接法は PC の余力に応じた診断に限る。

## 次タスクの Rust 実験案（本 PR では未実装）

1. 独立実験ブランチで、`src/stages.rs` の縮約後・`Tangent::solve` の直前と `src/linalg.rs` のスケーリング前に観測点を設ける。xy1 の零種・0.40・0.44 付近、失敗反復直前の A/b、全自由度対応、外力ノルム、確定/試行 u、材料旧状態、固定 partition、Newton 番号、進捗、入力・ソース・実行版 SHA を保存する。確定点の再計算を失敗反復の保存と混同しない。
2. まず同じ凍結系で現行 ILU(0) と Python ILU(0) の自然順序・構造ゼロ・行最大スケーリングを照合する。丸めや Arnoldi 再直交化が違うため反復数の完全一致は必須としない。元 A/b の残差は必ず照合する。
3. 同じ GMRES 右前処理インターフェイスへ ILUT（drop=1e-4、fill 目安10）と COLAMD 等の順列を実験導入する。順列は解と残差を元順へ戻す。小さいピボットを無断で正則化しない。分解失敗、nnz/時間/メモリー、再構築頻度を報告する。現行 ILU(0) を対照として保持する。
4. 同一物理入力の xy1 で ILUT と直接法の Newton 補正・線探索・受理/棄却を比較し、その後 PC の xy2 で進捗 0.44～0.4717 と 1.0 完了までを確認する。非対称・不定になる塑性接線を CG/SPD 前提の AMG にそのまま渡さない。AMG は別途、剛体 near-nullspace（並進/回転）、2～3 成分結合、非対称 Krylov 外側、粗格子補正を設計してから比較する。
5. 合格は線形全 solve の元残差 `||A du-b||/||b|| <= 1e-10` と、受理全増分の本体 `||r_free||/max(||F||,1) <= 1e-8`、geostatic/本体完了、履歴・棄却保持・再開照合をすべて満たすこと。線形診断の `force_normalized_linear_residual` は非線形力残差とは別量。元失敗の単なる前進や直接法の成功だけでは全条件の受入にしない。NC-01・P6全体・元 Griffiths 2D 変位誤差の未達は別に残す。
