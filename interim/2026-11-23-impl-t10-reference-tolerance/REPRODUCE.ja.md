# T10 の小型比較の再現

以下は Windows PC 用の手順。今回の PC 実行は0件。リポジトリのルートで PowerShell を開く。使用する Git コミットを記録し、通常の Git checkout で同じコミットの原票を用意する。Python 3.13、Rust は rust-toolchain.toml の版を使う。

```powershell
$ErrorActionPreference = 'Stop'
$PSNativeCommandUseErrorActionPreference = $true
git rev-parse HEAD
py -3.13 -m venv .venv-t10
$py = '.\.venv-t10\Scripts\python.exe'
& $py -m pip install -r requirements-reference.lock.txt
$env:OPENBLAS_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$env:CARGO_BUILD_JOBS = '2'
$env:PYTHONUTF8 = '1'
cargo build --locked --release --bin soil-fem-3d --example joint_birth_points --example t10_layered_cantilever
$study = 'interim/2026-11-23-impl-t10-reference-tolerance'
$out = Join-Path 'runs' ('t10-' + (Get-Date -Format 'yyyyMMdd-HHmmss'))
New-Item -ItemType Directory -Path $out | Out-Null
& $py "$study/t10_compare.py" "$out/comparisons.json"
& $py -m unittest discover -s tests -p test_t10_reference_tolerance.py -v
```

|項目|PC で扱う規模と比較内容|
|---|---|
|B-01|8・32・128 Tri6 の保存測定を再集計。規則・歪みごとに荷重点変位を比較し、欠測の診断を保持する。|
|B-02|10 m × 5 m、各層100 Tri6、計200要素・231共有節点。固定辺を除く共通自由度で線形系を1回解く。A3 の DZ と A1 の各層の QX、基準面まわりの MXX を読む。|
|B-03|1 m の二次辺、各翼2 Tri6。剛・ヒンジ・ばねを順に解き、同じ入力の独立乗数解と比較する。|
|C-04|乾燥の1点構成則を計算。湿潤は12 Tet4・6接触点、13工程・416受理点の本体と独立解の原票を再読込みする。|
|A-03|12 Tet10 の切断・配置と12 Tet4 の過渡施工、固定領域の土柱の保存場を再読込みし、物理秒ごとの収支を照合する。|

各経路を逐次実行する。最大の新規構造計算は B-02 の200要素で、全剛性行列1枚は約15.4 MB。プロセス全体のピーク RAM をタスクマネージャーで確認し、40 GB または1経路1日の上限に達する前に中止する。768／6144要素を使う検証はこの手順の対象外。

原票は T8 の `observations.zip` と `observations.sha256`、T2 の `evidence.zip` と `acceptance.json`、C-04 旧検証・原因調査の分割 ZIP と索引を使う。比較器は読み込む原票のハッシュを確認する。結果 JSON の `passed: null` は個別幅の保留に対応する。B-01 の停止行は `diagnostic`、A-03 の定常収支は `steady_domain_balance` を見る。ゼロ参照のパッチは物理単位の絶対差で判定する。
