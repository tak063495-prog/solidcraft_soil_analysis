# A-04原因診断の再現とPC追加手順

PR #13の同じブランチ、起点1c921b1df0a3da6784476a81c68e6353e39b3678の子孫を使う。
マージしない。全出力は未存在のパスを指定し、旧入力・値・失敗・封印を上書きしない。
A-04は未達、passed=false。新しい診断は明示実行だけ。通常CLIへ自動接続しない。

## Windows準備（PC未実行）

GPD WIN Max 2、使えるメモリ約48GB。追加計算は既存12 Tet4、乾燥24 Tet10、静水圧12 Tet10だけ。
新たな768/6144要素、見込み40GB超、1経路1日超のPC手順は作らない。
封印/ソースのバイト照合用の新規checkoutはCIと同じcore.autocrlf=falseで作る。
文書末尾の試験自体はCRLFでも一致し、封印SHAは元バイトのまま検査する。
この上限を超える確認が必要なら「規模上限のため確認不可・未達」。許容値は変えない。
全旧証拠checkout約28.4GBに加え、ビルド/展開/出力用のディスクを確保する。
これはディスク量であり、常駐メモリ量ではない。

```powershell
git status --short
git fetch origin codex/plastic-spatial-convergence
git switch codex/plastic-spatial-convergence
git merge-base --is-ancestor 1c921b1df0a3da6784476a81c68e6353e39b3678 HEAD
if ($LASTEXITCODE -ne 0) { throw '指定起点の子孫ではありません' }
$Py = '.venv/Scripts/python.exe'
if (-not (Test-Path $Py)) { py -3.13 -m venv .venv }
& $Py -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw '依存配置失敗' }
$env:CARGO_BUILD_JOBS='2'
$env:RUST_TEST_THREADS='2'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
if ($LASTEXITCODE -ne 0) { throw 'toolchain配置失敗' }
cargo build --locked --release --all-targets
if ($LASTEXITCODE -ne 0) { throw 'ビルド失敗' }
$Diag = 'interim/2026-12-23-a04-richards-wet-joint-cause-diagnosis'
```

ビルド予算10〜60分、概算メモリ2〜8GB（Windows未実測）。必要な証拠partがsparse checkoutで
欠落していないか確認する。欠落なら固定Git blob/SHAで復元し、skip/空ファイルで回避しない。

## 小型診断の全行列

```powershell
$env:A04_PROCESS_BUDGET_SECONDS='3600'
$env:A04CAUSE_PROCESS_BUDGET_SECONDS='7200'
try {
  & $Py "$Diag/a04cause_driver.py" --output runs/NEW-a04-cause-PC
  if ($LASTEXITCODE -ne 0) { throw '診断未完了。failure/process/log/途中状態を保持' }
  & $Py "$Diag/a04cause_report.py" --run runs/NEW-a04-cause-PC --output runs/NEW-a04-cause-PC-report.json
  if ($LASTEXITCODE -ne 0) { throw '集計未完了。求解成功と区別' }
} finally {
  Remove-Item Env:A04_PROCESS_BUDGET_SECONDS -ErrorAction SilentlyContinue
  Remove-Item Env:A04CAUSE_PROCESS_BUDGET_SECONDS -ErrorAction SilentlyContinue
}
```

計画枠10〜60分、同時1経路、概算ピーク1〜6GB、子プロセス上限1時間、親工程上限2時間。
Windows未実測の予算であり、完走保証ではない。Richards＋JOINTの拒否先へ進める計算ではない。
12 Tet4のVG/飽和でN4/8/16・4秒・2秒DEATH、別の乾燥/静水圧JOINT、同入力/キー逆順、
入力のみの独立Jacobi3/5/7積分感度、旧16拒否、旧境界不整合漏れを分離する。
既存内部スクリプトの復旧・再hash改変検査も行い、新規完走経路へ加算しない。

|出力|見る値|
|---|---|
|`source-binding.json`, `source-snapshot/`|求解前のソースSHAとコピー。途中で差替えない|
|`soil/bin/`, `soil/binaries.json`、各JOINT`bin/`|固定実行版SHA。同一バイナリで生成/復旧/内容改変を確認|
|`driver-result.json`, `*.process.json`, `*.log`, `failure.json`|実コマンド、PID、終了/時間切れ、完了した工程と準備失敗|
|`soil/VG-n*`, `soil/saturated-n*`|生入力、全Run、全受理場、旧定義49場の比較、実checkpointと後続/改変|
|`*-repeat`, `*-reverse-keys`|時間以外のRun厳密一致。入力の物理値は同じ、丸め差を許容しない|
|`*-quadrature5/7`|独立だけの積分変更。共通セル/面水頭・流量・水量・積分力・節点場の感度。異なるGPを付け替えない|
|`dry/`, `hydrostatic/`|別々の独立3D平衡と面荷重・世代/開閉。水量/流量/物理時計のnullを保持|
|`refusals/`|16件すべて終了1・結果なし・時間切れなし。拒否を対象機能の受入へ数えない|
|`boundary-gap/`|本体完走と独立の初期検査拒否。既知拒否漏れを未達として保持|
|`NEW-a04-cause-PC-report.json`|固定順L2/絶対差/生分母、物理dt、積分感度、凍結・符号/単位/力会計、interface各null、passed=false|

元の相対1e-7・節点力1e-10 kN・幾何/拘束1e-10 m、本体水理scaled残差1、独立.01、
力/真の線形残差と支持・拘束を保持する。旧SPECの全場包括基準未定義も残す。
同じ16物理時刻で比較し、2秒の4つの解放progressを物理秒へ読み替えない。
土収支は水量＋除去水−初期水−正味流入。界面の台帳へ転用しない。

WindowsのNumPy longdoubleがLinuxの拡張精度と異なり、独立側の厳しい残差.01を
満たせない場合は、生値/反復/途中出力を失敗記録として保持する。参照値・許容値を変えて通さない。

## 任意の長系列：同じ12 Tet4だけ（PC未実行）

```powershell
$env:A04_PROCESS_BUDGET_SECONDS='21600'
$env:A04CAUSE_PROCESS_BUDGET_SECONDS='21600'
try {
  & $Py "$Diag/a04cause_driver.py" --steps 32 64 128 --output runs/NEW-a04-cause-long-PC
  if ($LASTEXITCODE -ne 0) { throw '長系列未完了。途中/失敗を保存' }
  & $Py "$Diag/a04cause_report.py" --run runs/NEW-a04-cause-long-PC --output runs/NEW-a04-cause-long-PC-report.json
  if ($LASTEXITCODE -ne 0) { throw '長系列集計失敗' }
} finally {
  Remove-Item Env:A04_PROCESS_BUDGET_SECONDS -ErrorAction SilentlyContinue
  Remove-Item Env:A04CAUSE_PROCESS_BUDGET_SECONDS -ErrorAction SilentlyContinue
}
```

物理4秒をdt=1/32,1/64,1/128秒へ細分化する。空間・保持則・境界値・力学解放4等分はそのまま。
このオプションでは土6経路だけを両実装で解き、乾燥/静水圧・キー順・積分対照・16拒否は再反復しない。
最細のcheckpoint復旧/改変は元studyのまま含む。同じ128物理時刻だけを比較する。
計画枠は全体30分〜6時間、1経路5〜90分、メモリ概算1〜8GB（小型からの参考外挿、Windows未実測）。
時間切れなら未完了・未達を記録し、1経路1日超やメモリ40GB超へ延長しない。
追加後も有限dt/固定空間の対照であり、Richardson補正や次数推定で厳密参照・合格を作らない。

C-04接触小刻み/全体経路の交差対照、土の空間h収束は別の診断である。
今回のコマンドがそれらを計算したとはしない。未定義の連成仕様はPCへ回しても解消しない。
一般3D空間/界面参照の確認に上限超が必要なら「規模上限のため確認不可・未達」。

## 保存証拠・全件検証

```powershell
& $Py "$Diag/a04cause_evidence.py" --verify $Diag --extract runs/NEW-a04-cause-unsealed
if ($LASTEXITCODE -ne 0) { throw 'part/ZIP/entry SHA不一致' }
& $Py "$Diag/a04cause_report.py" --run runs/NEW-a04-cause-unsealed/formal --output runs/NEW-a04-cause-remeasured.json
if ($LASTEXITCODE -ne 0) { throw '保存場再集計失敗' }
& $Py interim/2026-12-06-p6-01-status-reconciliation/audit_status.py
if ($LASTEXITCODE -ne 0) { throw '状態台帳integrity不合格' }
cargo fmt --all --check
Write-Output "fmt exit=$LASTEXITCODE"
cargo clippy --locked --all-targets -- -D warnings
Write-Output "clippy exit=$LASTEXITCODE"
cargo test --locked --release
Write-Output "Rust main exit=$LASTEXITCODE"
cargo test --locked --release --examples
Write-Output "Rust examples exit=$LASTEXITCODE"
foreach ($group in @('python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence','python-c01-switch','python-c01-path','python-c01-cause')) {
  & $Py tools/run_ci_python_evidence_shards.py --group $group --output "runs/NEW-a04-validation/$group"
  Write-Output "$group exit=$LASTEXITCODE"
}
```

全件はRust10〜90分、Python全群30〜120分、メモリ概算2〜8GBの計画枠（Windows未実測）。
発見/実行/成功/失敗/error/未実行/skip/ignore/filterを各inventoryと生ログで照合する。
最初の失敗や準備ログを保持し、再実行の成功を同じIDへ二重加算しない。
文書はuniversal newlineで読み、実CRLF試験も行う。ノルム順序の差を許容値で吸収しない。
旧NCの過去失敗や8件未実行を書き換えず、今回の実測と区別する。
保存証拠の照合/再集計は新規求解・物理受入ではない。

Linuxは`source /workspace/soil-environment/activate.sh`後に同じPythonコマンドを使う。
push後のWindows CIは未確定。A-04、A-01〜A-03、B各項、C-01〜C-04、NC-01、P6は未達。
