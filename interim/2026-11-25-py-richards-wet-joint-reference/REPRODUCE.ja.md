# 再現：独立Richards湿潤JOINT参照

A-04/C-04は未達。保存先ラベル2026-11-25、実施2026-11-12（日本時間）。本参照の数学的問題、出典、未確認、判定の範囲は [README.ja.md](README.ja.md)。全コマンドはリポジトリの根から実行し、出力は未存在名を使う。失敗したフォルダーを消して再利用しない。

## Windows PCの準備

Windows 11、PowerShell、Python 3.13。小型参照はCPU2コア/RAM4 GB/空き1 GB程度。既存全件・旧封印証拠まで行うならRAM16 GB以上/空き80〜100 GBを計画する。初回取得時間は回線依存。以下は既存checkoutを想定する。

```powershell
git fetch origin codex/py-richards-wet-joint-reference
git switch codex/py-richards-wet-joint-reference
git pull --ff-only
# 未コミット変更のない状態で実施し、実際に検証する先端を記録。
git status --short
git rev-parse HEAD
git rev-parse HEAD^
# このPRは指定79d839e3...を直接の親にする。以後の他作業を混ぜない。
py -3.13 -m venv .venv-rj
.\.venv-rj\Scripts\Activate.ps1
python -m pip install -r requirements-reference.lock.txt
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
$Study='interim/2026-11-25-py-richards-wet-joint-reference'
$Run='work/NEW-PC-richards-joint-small'
```

NumPy/SciPyの版は固定requirementsと実環境を記録する。状態水頭・貯留演算にlongdoubleを使い、仮数ビット数も記録する。Windowsネイティブは通常52ビット、今回Linuxは63ビット。小流量台帳や旧証拠の厳密再現は同じ依存版のWSL/Linuxでも確認する。Windowsの丸め差/台帳拒否は失敗のまま保持し、許容値を緩めない。新参照の求解にRustビルドは不要。旧証拠まで含む既存全件はRustと既存CLIを必要とする。

## 小型参照・診断・新規試験

```powershell
python "$Study/inspect_inputs.py" --output work/NEW-rj-input-inventory.json
if ($LASTEXITCODE -ne 0) { throw 'cited input/source differs' }
python "$Study/study.py" --output $Run
if ($LASTEXITCODE -ne 0) { throw 'reference incomplete: retain failure and partial results' }
python "$Study/validate.py" --output work/NEW-rj-validation
if ($LASTEXITCODE -ne 0) { throw 'new software checks failed' }
```

Linux（環境内では既存venvを使用）:

```bash
source /workspace/soil-environment/activate.sh
Study=interim/2026-11-25-py-richards-wet-joint-reference
python "$Study/inspect_inputs.py" --output work/NEW-linux-rj-inputs.json
python "$Study/study.py" --output work/NEW-linux-rj-small
python "$Study/validate.py" --output work/NEW-linux-rj-validation
```

小型行列は46経路、1D最大16セル、2D最大128セル、各0〜1秒（固有モード対照は.02秒）。本クラウドLinuxの最終行列は約157秒、Windowsの初期予算は2〜10分、RAM4 GB/出力数十MB程度（Windows未実測）。新規44試験はLinux約4秒（同時実行の影響あり）、Windows10〜60秒程度を初期予算とする。時間切れ/例外/キャンセル/部分出力を成功に数えない。

単独条件を実行する例:

```powershell
python "$Study/reference.py" --nx 8 --ny 3 --law vg --mode TWO_WAY --policy CONTINUOUS --dt .025 --output work/NEW-PC-rj-2d-two-way
python "$Study/reference.py" --nx 8 --ny 1 --law vg --mode ONE_WAY --policy PRESCRIBED --dt .025 --output work/NEW-PC-rj-1d-prescribed
python "$Study/reference.py" --nx 8 --ny 1 --law vg --mode TWO_WAY --dry --dt .025 --output work/NEW-PC-rj-dry-mechanical-control
```

`--dry`は力学の水圧を0にする対照で、水理の無水化ではない。SOIL過渡・接触が起動する双方向PRESCRIBED・世代変更は未確認で、対応済みのオプションとして追加しない。

## 出力と判定

|出力|内容と確認|
|---|---|
|`summary.json`|46経路の完了、各ケースの収支/流量/状態、時間/格子の生差と観測次数、±.0125秒の開閉時刻感度、ONE_WAY/双方向の差。`accepted=null`, A04/C04未達。|
|`vg-*.json`, `linear-*.json`, `dry-*.json`|1D/2D、CONTINUOUS/PRESCRIBED、乾燥力学対照。実物理時刻の全セル/接合水頭・全貯留・流量・接触履歴。|
|`time-*.json`, `space-*.json`, `space2d-*.json`|時間/格子を分離。格子値の制限は明示した体積平均だけ。異なる接合点を黙示補間しない。|
|`eigen-*.json`|閉鎖飽和の連続PDE固有モード対照。固定時間誤差を残した生のRMS/次数。|
|`hydraulic_frames`|各時刻の全主水頭、実貯留、外面外向き流量、水源流入、累積正味流入、収支欠損、生残差。初期主水頭も保存。|
|各`rows`|両土面トレース/接合中心H、a→joint/joint→b流量、薄層貯留率/水源流量、水理開閉と接触状態、世代0、gap/有効接触圧/物理圧/せん断/塑性すべり/散逸。|
|`contact_switch_brackets`|解に依存する接触開閉の物理時刻区間。正確なイベント根としない。|
|`failure.json`|数値未完了/例外を保持。完了した他の出力があっても行列成功にしない。|
|`validation.json`・全ID/ログ・前後SHA|新規44件の実実行件数と状態、sourceの前後一致。skip/xfail/失敗を成功に加算しない。|

ソフトウェア確認の合否は全件実行・失敗/error/skip/expected-failure/unexpected-success=0、ID重複欠落0、ソース前後一致。数値の物理受入はこれと別。既存C-04の相対流量1e-8・絶対1e-16 m³/s・水頭1e-8 mは元の適用範囲で維持する。比較器の符号付き流量台帳には、既存保存式と同じ絶対＋相対×絶対値和を使う（src/richards.rs）。全連成の受入許容値は未定義なので**accepted=null**。新しい観測値から許容値を作らず、負/低次数や差を保存する。A-04/C-04を完了にしない。

## 将来の本体CSV/JSON

現行本体のRichards＋JOINTは拒否され、本体側adapter/exporterは未確認。以下は将来、同じ検証入力/全時刻/全IDを出力できた際の比較コマンドである。本体schemaや実装を本PRでは追加しない。

```powershell
# 新参照を交換契約へ包むだけ。本体exporterではない。
python "$Study/pack.py" "$Run/vg-1d-ONE_WAY-CONTINUOUS.json" work/NEW-rj-reference.json
python "$Study/pack.py" "$Run/vg-1d-ONE_WAY-CONTINUOUS.json" work/NEW-rj-reference.csv
# CSVには work/NEW-rj-reference.metadata.json が自動生成される。

# 自己比較の契約検査。正常でも未定義の数値受入なのでexit 2が正しい。
python "$Study/compare.py" work/NEW-rj-reference.json work/NEW-rj-reference.csv work/NEW-rj-contract-report.json
if ($LASTEXITCODE -ne 2) { throw 'exchange contract check failed' }

# 将来の生本体出力を外部adapterで明示変換できた場合だけ実行。
# このパスはユーザーが用意するファイルで、本PRが生成済みとはしない。
python "$Study/compare.py" work/NEW-rj-reference.json work/FUTURE-body-joint.json work/NEW-rj-body-comparison.json
if ($LASTEXITCODE -eq 1) { throw 'body contract/ID/time/sign mismatch; do not interpolate' }
if ($LASTEXITCODE -ne 2) { throw 'unexpected comparator exit' }
```

`producer=future_body_export`を実本体の外部変換にだけ指定する。人工例を本体証拠としない。CSVは同名`.metadata.json`必須。単位・符号・時計・連成方式・全ID・世代・補間なしを明示し、入力SHAとともに照合する。比較差が0でも採否はnull。出力の時刻数が異なる場合、無断補間や最寄り時刻合わせで通さない。

## PC専用：長い時間・格子系列

クラウドでは実行していない `--pc` は、時間刻み.00625/.003125/.0015625秒、1D格子32/64/128（固定dt=.00078125秒）、2D格子(16,8)/(32,16)/(64,32)（固定dt=.0015625秒）を使う。物理時間は同じ1秒で、別の長期物理現象を作らない。16 GB RAM、空き5 GB以上、初期予算2〜8時間を見込む（未実測、完了保証なし）。

```powershell
$PcRun='work/NEW-PC-rj-fine-series'
# 8時間を上限にし、終了コードとは別にタイムアウトを記録する。
$PcLog='work/NEW-PC-rj-fine-series.stdout.log'
$PcErr='work/NEW-PC-rj-fine-series.stderr.log'
$Proc=Start-Process -FilePath '.\.venv-rj\Scripts\python.exe' -ArgumentList @("$Study/study.py",'--pc','--output',$PcRun) -PassThru -RedirectStandardOutput $PcLog -RedirectStandardError $PcErr
if (-not $Proc.WaitForExit(28800000)) {
  Stop-Process -Id $Proc.Id
  '{"complete":false,"reason":"PC timeout; not accepted"}' | Set-Content -Encoding utf8 work/NEW-PC-rj-fine-series-timeout.json
  throw 'PC reference timed out; preserve all partial data'
}
$Proc.Refresh()
if ($Proc.ExitCode -ne 0) { throw 'PC reference failed; preserve logs' }
```

出力と合否の考え方は小型と同じ。次数低下、開閉ブラケット、累積水量、モード差を報告し、最細値への差が小さいだけで合格としない。停止はこの自分が起動したPIDだけを対象にし、並行の杭/NC等のプロセスへ触れない。

## 大きな3D・本体C-04はPCで扱う

この新規参照は1D/2Dまでで、独立3Dのバルク/接合世代連成は**未実装・未確認**。同じ交換契約だけで大きな3D GP群へ補間して比較したとはしない。大規模3Dの入力・写像・独立参照・同定/採否基準が必要であり、未対応の本体Richards-JOINTを成功させるコマンドは作らない。

既存C-04に用意された3D Tet10の長い規定水位系列と**別の水理単独系列**をPCで確認する場合は、元の手順のコマンドをそのまま使える。これは本参照との同一問題比較でも、Richards-JOINT双方向の3D受入でもない。今回クラウドで再実行しない。

```powershell
rustup toolchain install 1.99.1 --profile minimal
cargo build --locked --release
if ($LASTEXITCODE -ne 0) { throw 'CLI build failed' }
cargo build --locked --release --example wet_joint_transient_probe
if ($LASTEXITCODE -ne 0) { throw 'observer build failed' }
$C04='interim/2026-11-22-wet-joint-transient-water'
python "$C04/run_all.py" --output work/NEW-PC-C04-tet10-long --probe target/release/examples/wet_joint_transient_probe.exe --solver target/release/soil-fem-3d.exe --element tet10 --mechanical-steps 32 --water-steps 32 64 128 --hydraulic-steps 16 32 64 --timeout 14400
if ($LASTEXITCODE -ne 0) { throw 'C04 incomplete; retain all failures' }
python "$C04/saturated_control.py" --output work/NEW-PC-C04-tet10-long --steps 16 32 64 --timeout 14400
if ($LASTEXITCODE -ne 0) { throw 'C04 PDE control incomplete' }
```

元REPRODUCEの初期予算：RAM32 GB、空き20 GB以上、1経路.5〜4時間、全行列半日〜数日（未実測）。同例は54節点・12 Tet10までであり、一般の大規模3Dメッシュ生成器を備えたとはしない。さらに大きい実CAD/3Dは既存の入力準備と独立参照の設計から必要で、A-04/C-04/P6未達を維持する。

## Linux/Windowsの既存Python全件

必要な旧証拠を元Git blobとして取得する。全件は元CIの4群であり、この新規チェックを既存収集へ混ぜない。疎checkoutを使う場合も、試験が読むarchive/source snapshot/measurementsは省略しない（未配置エラーをskipで通さない）。

```powershell
# 旧証拠をすべて配置する。回線/ディスク依存で数十分以上。
git sparse-checkout disable
cargo build --locked --release
if ($LASTEXITCODE -ne 0) { throw 'build failed' }
python "$Study/validate.py" --existing --output work/NEW-PC-rj-all-validation
if ($LASTEXITCODE -ne 0) { throw 'some discovered tests failed or did not execute; inspect every log' }
```

または4群を直接実行してそれぞれの終了コード/ログを保持する：

```powershell
python tools/run_ci_python_tests.py --group python-core-0 --output work/NEW-rj-core0
python tools/run_ci_python_tests.py --group python-core-1 --output work/NEW-rj-core1
python tools/run_ci_python_tests.py --group python-mean-evidence --output work/NEW-rj-mean
python tools/run_ci_python_tests.py --group python-spatial-evidence --output work/NEW-rj-spatial
```

Linuxも同じPythonコマンド。Windowsは60〜120分以上を初期予算とし、封印展開/初回ビルド/ディスク性能による追加を見込む。旧NCの2群は保存証拠の再計測で、新規NC求解やNCの数値受入ではない。setUpClassエラーなら各群の4メソッドは未実行で、成功件数へ加算しない。失敗したモジュールを再実行する場合も、元ログを残し、成功ID和集合で件数を集計する。全件成功でなければ明記し、Windows CIの未完了・キャンセル・時間切れを成功としない。


## 保存された46経路の生値を展開

圧縮によってPRの巨大なJSON差分を避け、集計/交換例だけを直接表示する。アーカイブは全46JSON（各時刻の全主水頭・貯留・接合点履歴）を損失なく保持する。再求解とは別のバイト検証：

```powershell
$Manifest=Get-Content "$Study/results/evidence-manifest.json" -Raw | ConvertFrom-Json
$Digest=(Get-FileHash "$Study/results/evidence.zip" -Algorithm SHA256).Hash.ToLower()
if ($Digest -ne $Manifest.archive_sha256) { throw 'archive SHA differs' }
python -m zipfile -t "$Study/results/evidence.zip"
if ($LASTEXITCODE -ne 0) { throw 'archive CRC differs' }
$Stored='work/NEW-rj-stored-evidence'
if (Test-Path $Stored) { throw 'choose a new extraction directory' }
python -m zipfile -e "$Study/results/evidence.zip" $Stored
if ($LASTEXITCODE -ne 0) { throw 'archive extraction failed' }
foreach ($Entry in $Manifest.entries) {
  $File=Join-Path $Stored $Entry.name
  if ((Get-Item $File).Length -ne $Entry.bytes) { throw "entry size differs: $File" }
  if ((Get-FileHash $File -Algorithm SHA256).Hash.ToLower() -ne $Entry.sha256) { throw "entry SHA differs: $File" }
}
```

Linuxは `python -m zipfile -t` / `-e` と `sha256sum` を使える。最終ソースSHAはmanifestと新規検証の前後指紋に保存する。開発中の最初の台帳拒否と当時の全生値/ソースは `results/development.zip` に別名で残し、最終値に混ぜない。
