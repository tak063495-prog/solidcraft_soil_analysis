# A-04 の Windows 再実行

対象は PR #13 の本記録を含むコミット。GPD WIN Max 2（利用可能メモリ約 48 GB）で、旧 failed ID を保存したまま新規出力を作る。以下はクラウドから PC に対して実行していない。Git の固定証拠を検査できる checkout を使う。初回 checkout 前に core.autocrlf=false を指定する。既存 checkout の改行を一括変換しない。

## 準備と軽量回帰

リポジトリ直下の PowerShell で実行する。既存の固定 requirements と Rust toolchain を使う。追加の数値ライブラリは不要。すべての旧ソースが必要になるため、少なくとも git ls-files に含まれる .py、.rs、.ps1、.toml、SPEC.ja.md、Cargo.lock は実体を置く。

```powershell
$ErrorActionPreference = 'Stop'
$Base = '6b147977ef80f1ac41a9c9cfb14cd1eb63c263b2'
git merge-base --is-ancestor $Base HEAD
if ($LASTEXITCODE -ne 0) { throw 'PR #13 の指定先端を含む履歴が必要' }
$D = 'interim/2027-01-12-a04-fixed-order-portable-reference'
$Old = 'interim/2026-12-23-a04-richards-wet-joint-cause-diagnosis'
$Py = '.venv/Scripts/python.exe'
if (-not (Test-Path $Py)) {
  py -3.13 -m venv .venv
  if ($LASTEXITCODE -ne 0) { throw 'venv 作成失敗' }
}
& $Py -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw '固定依存の導入失敗' }
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:RUST_TEST_THREADS='2'
$env:CARGO_BUILD_JOBS='2'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
cargo build --locked --release --bin soil-fem-3d --example a04_checkpoint_rehash --example joint_a02_history_probe --example wet_birth_a03_probe --example joint_ext_checkpoint_rehash
if ($LASTEXITCODE -ne 0) { throw '本体・既存プローブのビルド失敗' }
& $Py -m unittest discover -s tests -p test_a04ref_portable.py
if ($LASTEXITCODE -ne 0) { throw 'binary64 強制の軽量回帰失敗' }
```

## 新規 ID で元の matrix を再実行

水理は従来の 12 Tet4、既存 matrix の乾燥・静水圧対照も元の小型入力を使う。刻みは 4/8/16、直交積分は元の 3/5/7。旧 driver のプロセス予算（親 1800 秒、study の子は既定 900 秒等）を引き継ぐ。既存環境の診断用予算上書きがあれば別の初期 PowerShell セッションで行う。1 経路 1 日を超える手順、40 GB 超の経路、一般空間収束のための大規模入力を計画しない。

```powershell
$Id = 'a04-matrix-portable-' + [DateTime]::UtcNow.ToString('yyyyMMddTHHmmssfffffffZ')
$Run = Join-Path 'runs' $Id
if (Test-Path $Run) { throw '出力 ID が衝突' }
& $Py "$D/a04ref_matrix.py" --force-binary64 --output $Run
if ($LASTEXITCODE -ne 0) { throw "matrix 失敗。$Run のログ・途中状態を保持" }
& $Py "$Old/a04cause_report.py" --run $Run --output "$Run/a04ref-original-report.json"
if ($LASTEXITCODE -ne 0) { throw '元ゲートの報告未完了' }
Get-Content -Raw "$Run/driver-result.json" | ConvertFrom-Json | Select-Object execution_complete, passed, a04_complete
```

soil/* と VG/saturated の reverse-keys、quadrature5/7 各フォルダの a04ref-line-search.json に全試行の残差、半減回数、反復番号がある。boundary-gap の独立側は元の想定拒否で、その理由は元 report が検査する。soil/measurements.json の cases ごとに existing_bounded_gates_passed を確認する。execution_complete、元の部分ゲート、A-04 の全体判定は別の欄を読む。失敗時は failure.json、各 *.log、*.process.json と最終反復を含む出力全体を新規 ID のまま保存する。

## クラウドの実測を再現する入口

リポジトリ直下。すべて未存在の出力先を指定する。Linux long double の対照だけは仮数が拡張された環境を使う。乾燥経路は 12-31 に使った既存の 24 Tet10、接触 32 点をそのまま使う。水理 12 Tet4 のメッシュへ置き換えると当該第 5 反復を比較できない。

```bash
D=interim/2027-01-12-a04-fixed-order-portable-reference
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
python "$D/a04ref_dry_build.py" --output work/NEW-a04ref-build
python "$D/a04ref_dry_run.py" --probe work/NEW-a04ref-build/target/release/a04ref_dry --output work/NEW-a04ref-dry
python "$D/a04ref_dry_report.py" --raw work/NEW-a04ref-dry --output work/NEW-a04ref-dry.json
python "$D/a04ref_host.py" --host work/NEW-a04ref-host
python "$D/a04ref_matrix.py" --force-binary64 --output runs/NEW-a04ref-matrix
python "$D/a04ref_host.py" --host work/NEW-a04ref-host --compare runs/NEW-a04ref-matrix --output work/NEW-a04ref-comparison.json
python interim/2026-12-23-a04-richards-wet-joint-cause-diagnosis/a04cause_report.py --run runs/NEW-a04ref-matrix --output work/NEW-a04ref-original-report.json
```

## 封印済み生証拠の読取り

本フォルダの archive-index.json にある全 part を配置する。Linux 全件用に以前から必要だった 1,076 パスは 01-09 の linux-results.zip 内 archive-paths.tsv に固定されている。sparse checkout で必要な実体が欠けた場合は、固定コミットの tree の blob ID と実体の git hash-object を照合する。partial clone の promisor だけに存在する場合は、必要なディスク量を確認して git fetch --refetch --no-filter を使うか、固定コミットの該当パスだけを取得して同じ blob ID・バイト数を検証する。欠落時は検査を中止する。

```powershell
$Verify = 'interim/2026-12-31-a04-richards-wet-joint-cloud-followup/a04follow_evidence.py'
& $Py $Verify --verify $D
if ($LASTEXITCODE -ne 0) { throw '封印照合失敗' }
```

展開する場合は --extract に未存在のディレクトリを指定する。formal 以下が実際の portable matrix、retained/dry-raw が E/F/Q ごとの全反復、retained/a04-observations/dry-iterations.json が各節点の固定順・fsum・厳密和・同一状態独立再組立である。retained/a04-host が Linux 既存参照の全場。封印の再照合に全経路の再計算は不要である。
