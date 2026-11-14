# C-04 の再現（Windows PC / Linux）

判定は **未達**。小型規定水位ケースの独立3D照合・Richards単独照合・実再開の完了と、C-04/A-04の受入を分ける。保存先ラベル2026-11-22、実施2026-11-12（日本時間）。PR #13 の同じブランチを使い、マージしない。

## 準備

Windows 11 / PowerShell、Python 3.13、Rust 1.99.1、CPU 4〜8コア、RAM 16 GB以上。全件試験で使う過去証拠も取得する。リポジトリ＋一時展開に空き80 GB程度を計画する。今回の新規小型データだけは数百MB程度である。

```powershell
git switch codex/plastic-spatial-convergence
git pull --ff-only
# 既存CIと同じ、試験に使わない旧postfold生アーカイブだけを省く。
git sparse-checkout set --no-cone '/*' '!/interim/2026-11-10-nc-postfold-branches/evidence.tar.gz.part*'
py -3.13 -m venv .venv-c04
.\.venv-c04\Scripts\Activate.ps1
python -m pip install -r requirements-reference.lock.txt
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
$env:CARGO_BUILD_JOBS='2'
$env:RUST_TEST_THREADS='2'
cargo build --locked --release
if ($LASTEXITCODE -ne 0) { throw 'CLI build failed' }
cargo build --locked --release --example wet_joint_transient_probe
if ($LASTEXITCODE -ne 0) { throw 'observer build failed' }
$Study='interim/2026-11-22-wet-joint-transient-water'
$Run='work/NEW-C04-small'
```

全件試験はmean/spatial等の旧生証拠も必要。取得未完了で試験を始めたり、NCの整合試験を省いたりしない。旧 `passed=false`、失敗、参照値、hash/snapshotを保持する。初回ビルド10〜30分程度、過去証拠の取得は回線に応じて追加時間を見込む。

## 小型行列を入力から新規求解

```powershell
python "$Study/run_all.py" --output $Run --probe target/release/examples/wet_joint_transient_probe.exe --solver target/release/soil-fem-3d.exe --mechanical-steps 32 --water-steps 8 16 32 --hydraulic-steps 4 8 16 --timeout 900
if ($LASTEXITCODE -ne 0) { throw 'C-04 matrix incomplete; retain every partial file' }
python "$Study/saturated_control.py" --output $Run --steps 4 8 16 --timeout 900
if ($LASTEXITCODE -ne 0) { throw 'saturated PDE control incomplete' }
python "$Study/diagnose.py" --output $Run --report "$Run/cause-diagnostics.json"
if ($LASTEXITCODE -ne 0) { throw 'cause diagnostics incomplete' }
```

出力先は未存在の名前にする。失敗後はフォルダーを消して再利用せず、別名へ保存する。Linuxでは実行版の `.exe` を外す。

```bash
cargo build --locked --release --example wet_joint_transient_probe
cargo build --locked --release
python interim/2026-11-22-wet-joint-transient-water/run_all.py \
  --output work/NEW-C04-linux --mechanical-steps 32 \
  --water-steps 8 16 32 --hydraulic-steps 4 8 16 --timeout 900
python interim/2026-11-22-wet-joint-transient-water/saturated_control.py \
  --output work/NEW-C04-linux --steps 4 8 16
python interim/2026-11-22-wet-joint-transient-water/diagnose.py \
  --output work/NEW-C04-linux --report work/NEW-C04-linux/cause-diagnostics.json
```

今回のLinux実測は主行列約424秒、飽和対照約41秒（ビルド・全件試験を除く）。Windowsでは10〜30分程度を初期予算にする。1子プロセス900秒の上限を超えたら未完了として保持する。未完了や時間切れを成功へ変更しない。

水位ケースの力学境界/初期平衡は32増分固定、水位工程だけ8/16/32。13工程の全前史、H=.75→1.25→.25→.9→.4→1.1→.75 m、CONTINUOUS、MC弾性バルクを使う。RichardsのVG物理時計は別の0〜7秒、Δt=.25/.125/.0625秒である。規定水位にはPDE質量収支がなく、RichardsとJOINTを勝手に接続しない。

独立水理の局所演算は `numpy.longdouble` を用いる。Linuxの厳密再現には同じ固定依存版のWSL/Linuxを使う。Windowsネイティブで残差/丸め停止が生じた場合も失敗を保持し、参照の `<.01` 基準を緩めない。Windowsでの新しい全系列は今回未確認である。

## 出力と確認基準

|出力|内容・確認|
|---|---|
|`identity.json`、`mesh.json`、`ground.json`、各`case.json`|起点SHA、実行版SHA、明示した入力。新規物理/許容値の変更なし。|
|`native-w*/raw/accepted.jsonl`・`constitutive.jsonl`|全受理状態と本体の実内部力/適用荷重・実分割。272/320/416状態、切戻し0、入力progress列そのもの。|
|`independent-w*/fields.jsonl`・`result.json`|入力だけから新規求解した全33場。保存JSON・本体分割数を解の入力にしない。|
|両側`fields.jsonl`|全節点u/反力/水荷重、全材料GPの有効/全応力・ひずみ・間隙比・水圧/ゼロ塑性履歴、全接触GPのjump/traction/塑性すべり/累積/仕事/開閉/圧力。|
|`hydraulic-h*/native.json`・`independent.json`|全受理物理時刻の全セル/面水頭・流量・貯留・収支。各28/56/112受理点、完了・棄却0。節点の接続面平均は診断投影で、主未知量と区別する。|
|`measurements.json`|全場・全工程・全点寄与、最大絶対差・L2、規定水位とPDE刻み感度を別欄に保存。|
|`saturated-controls/measurements.json`|既存Richardsの飽和極限と入力のみの級数解。代数一致を連続PDE精度へ読み替えない。|
|`scope-refusals/verification.json`|ONE_WAY＋CONTINUOUS/SOIL/PRESCRIBEDを実CLIで拒否、結果ファイルなし。未対応・双方向・世代を明記。|
|`mechanical-restarts/verification.json`|13実checkpoint・別PID、後続2704全状態とRun全体が時間以外完全一致。|
|`hydraulic-restarts/verification.json`|.5/2.5/5.5秒の実checkpoint・別PID、後続104/72/24全受理フレームとRun全体が時間以外完全一致。|
|`whole/raw/run.json`|観測なし直接求解と観測経路のRun全体が時間以外一致。|
|`cause-diagnostics.json`|局所診断、stick/slip不一致、共通28物理時刻での連続PDE差。新規求解へ数えない。|
|`completion.json`・`matrix-failure.json`・`failure-*.json`・ログ|完了・失敗・タイムアウトを区別し、削除しない。主行列完了でも`C04='未達'`, `passed=false`。|

元の本体力1e-8・真の線形1e-10・独立力1e-10、Richards相対流量1e-8/絶対1e-16 m³/s/水頭1e-8 mを維持する。参照RT0は同じ無次元残差式で `<.01` を要求する。再開は許容差ゼロの全状態・全Run一致（`elapsed_seconds`だけ除外）。欠けた場、カテゴリ差、非単調差、未完了を許容値やskipで隠さない。

全場の数値受入許容値は未定義であり、小さい差だけでC-04を合格にしない。接触仕事/カテゴリの差、連続PDEとの差、Richards→JOINT写像/物理時刻/開閉/世代/連成質量が残れば **C-04/A-04は未達**。C-01/C-02/C-03・NC-01・P6にも結果を転用しない。

## 保存点からの再開だけを再実行

```powershell
python "$Study/study.py" --mode restart --output $Run --source native-w32 --folder mechanical-restarts-SECOND --timeout 900
if ($LASTEXITCODE -ne 0) { throw 'mechanical real-file resume differs' }
python "$Study/study.py" --mode hydraulic-restart --output $Run --source hydraulic-h16 --folder hydraulic-restarts-SECOND --timeout 900
if ($LASTEXITCODE -ne 0) { throw 'hydraulic real-file resume differs' }
```

そのRunを作った同一SHAの実行版を使う。Linux保存点を別のWindows実行版へ付け替えず、Windowsでは上の新規行列から作った同一Windows実行版/保存点を使う。

## PC：大きめのGP数と長い系列

クラウドの正式系列は12 Tet4まで。既存12 Tet10は54節点・324材料GP・32接触GPで、新規の大規模メッシュ・CAD入力へ拡張した仕様ではない。力学刻み32を固定し、水位更新だけを長くする例：

```powershell
python "$Study/run_all.py" --output work/NEW-PC-C04-tet10 --probe target/release/examples/wet_joint_transient_probe.exe --solver target/release/soil-fem-3d.exe --element tet10 --mechanical-steps 32 --water-steps 32 64 128 --hydraulic-steps 16 32 64 --timeout 14400
if ($LASTEXITCODE -ne 0) { throw 'Tet10 matrix incomplete; keep failures' }
python "$Study/saturated_control.py" --output work/NEW-PC-C04-tet10 --steps 16 32 64 --timeout 14400
```

RAM 32 GB、追加空き20 GB以上を計画する。Tet10新系列は未実測で、1経路0.5〜4時間、全行列は半日〜数日を初期予算にする（完了保証ではない）。水理側はここでも別の12 Tet4であり、Tet10水理の独立参照を実装済みとはしない。本体/独立が切戻しや未収束なら、その時点で固定列の検証は失敗として残す。力学刻みまで変える別実験はC-01へ混入させない。

小型Tet4で長い水位/水理時間系列を測る場合：

```powershell
python "$Study/run_all.py" --output work/NEW-PC-C04-long --probe target/release/examples/wet_joint_transient_probe.exe --solver target/release/soil-fem-3d.exe --element tet4 --mechanical-steps 32 --water-steps 64 128 256 --hydraulic-steps 32 64 128 --timeout 7200
if ($LASTEXITCODE -ne 0) { throw 'long matrix incomplete' }
python "$Study/saturated_control.py" --output work/NEW-PC-C04-long --steps 32 64 128 --timeout 7200
```

1〜4時間を初期予算、RAM 16〜32 GB・追加空き10 GB以上を計画する（未実測）。これは同じ0〜7秒を細分化する系列で、長時間の別物理を作っていない。大規模3D連成・一般地質・実CAD・GUI/GPUには既存未対応の接続/独立参照が必要で、本ドライバから実行可能な合格コマンドを捏造しない。PCでも現行CLIの拒否をまず保持し、必要な入力・独立参照・資源と受入基準を別途定義する。

## 全件ソフトウェア検査

```powershell
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release
cargo test --locked --release
python tools/run_ci_python_tests.py --group python-core-0 --output work/NEW-C04-check/core0
python tools/run_ci_python_tests.py --group python-core-1 --output work/NEW-C04-check/core1
python tools/run_ci_python_tests.py --group python-mean-evidence --output work/NEW-C04-check/mean
python tools/run_ci_python_tests.py --group python-spatial-evidence --output work/NEW-C04-check/spatial
```

各コマンドが終了コード0であることを確認し、4群の発見ID/実行IDの和集合が重複・欠落なく一致すること、failure/error/ignored/skipped/filtered-out/expected-failure/unexpected-successが0であることを確認する。過去のNC証拠群も実行する。Linux成功・旧版CI成功を今回のWindows CI成功へ転用しない。

## 保存生証拠の再監査（新規求解とは別）

```powershell
python "$Study/archive.py" --mode replay --archive $Study --output work/NEW-C04-sealed-replay
if ($LASTEXITCODE -ne 0) { throw 'archive SHA or isolated measurement differs' }
```

全part・ZIP・全エントリーのサイズ/CRC/SHAを検査して新規展開し、封印ソースを隔離Pythonで実行する。`measurements.json`と原因診断がJSONオブジェクトとして完全一致することが基準。必要ディスクは数GB、数分程度を計画する。Windowsの浮動小数点差もそのまま失敗として保持し、厳密再現は同版WSL/Linuxで行う。これは新規求解・実checkpoint復旧を実行するコマンドではない。


今回Linuxの全件検証は初回Rustテストビルド・旧証拠の配置不足による初回群と再検証を含む監視プロセス全体で約49.3分。Windowsで4群を直列に実行する場合は60〜120分を初期予算とし、CPU/ディスク/初回ビルド差を見込む。全件実測値はvalidation.jsonを参照。タイムアウトを設けて成功扱いにはしない。
