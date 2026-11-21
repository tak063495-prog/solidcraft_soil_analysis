# C-03 の再現手順（Windows PC / Linux）

判定は **未達**。以下の正常完了は固定した塑性MCケースの統合観測・監査の完了であり、C-01/C-02/C-04・NC-01・P6や一般材料の受入ではない。保存先ラベル2026-11-20、実施2026-11-12（日本時間）。PR #13 の同じブランチを使い、マージしない。

## 環境と全証拠の取得

Windows 11/PowerShell、Python 3.13、Rust 1.99.1。小型Tet4の実行はCPU 4〜8コア・RAM 16 GBを目安とする。全件Pythonには既存NC保存証拠も必要で、リポジトリと一時展開を含め空き80 GB程度を計画する。新規C-03だけのデータは数百MB以下。GUI/GPU・実CADはこの数値ドライバの対象外で、実CAD入力・独立参照の定義がない状態で受入結果を作らない。

```powershell
git switch codex/plastic-spatial-convergence
git pull --ff-only
# 既存CIと同じ除外。未使用postfold保存系列以外を取得する。
git sparse-checkout set --no-cone '/*' '!/interim/2026-11-10-nc-postfold-branches/evidence.tar.gz.part*'
py -3.13 -m venv .venv-c03
.\.venv-c03\Scripts\Activate.ps1
python -m pip install -r requirements-reference.lock.txt
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
$env:CARGO_BUILD_JOBS='2'
$env:RUST_TEST_THREADS='2'
cargo build --locked --release --example wet_joint_plastic_probe
if ($LASTEXITCODE -ne 0) { throw 'probe build failed' }
$Study='interim/2026-11-20-wet-joint-plastic-bulk-integration'
```

旧証拠の取得前に試験を開始しない。NCの試験・整合検査を除外せず、元の `passed=false` やhash/snapshotを変更しない。過去のRust保存記録テストが読むFAILED.jsonも必要。全件検査の初回ビルドには10〜30分程度、旧証拠の取得は回線に応じてさらに時間を見込む。

## 小型塑性ケースを入力から新規求解

```powershell
$Run='work/NEW-C03-tet4-g32'
python "$Study/run_all.py" --output $Run --solver target/release/examples/wet_joint_plastic_probe.exe --steps 32 --element tet4 --timeout 900
if ($LASTEXITCODE -ne 0) { throw 'C-03 matrix incomplete; keep all failure files' }
```

新規フォルダーのみ受け付ける。既存フォルダーを消して再実行せず、異なる名前へ保存する。Linuxでは `.exe` を外す。通常のcheckoutでのLinuxコマンドは次の通り。

```bash
cargo build --locked --release --example wet_joint_plastic_probe
python interim/2026-11-20-wet-joint-plastic-bulk-integration/run_all.py \
  --output work/NEW-C03-linux-g32 \
  --solver target/release/examples/wet_joint_plastic_probe \
  --steps 32 --element tet4 --timeout 900
```

初期化は既存12 Tet4・16節点・2層入力から作り、両MCのc=30 kPa、open係数=.12へ明示変更する。それ以外のE/ν/φ/ψ・初期応力/間隙比・水位・JOINT則・拘束方向は既存入力に従う。元open=1.0の失敗を上書きした再現ではない。全体は各7工程1/32、材料内上限1e-4/5e-5/2.5e-5、接触内上限1e-4/5e-5/2.5e-5を別々に変える。加えて分割1対照を解く。

本体6・入力のみの独立系6経路、基準本体の各7工程途中の実ファイル再開、観測なし対照、全場監査を順に実行する。今回Linuxでは合計約271秒（ビルドを除く）。Windows PCでは5〜20分程度を初期予算とし、1経路900秒を超えたら失敗として残す。時間切れ・未完了を成功にしない。

主出力：

- `mesh.json`・`ground.json`・`stages.json`・`identity.json`：入力と実行版のSHA。
- `source-sha256.json`：求解時の数値/試験ソース。不変を確認する。
- `*-native/raw/accepted.jsonl`：全受理点の本体確定状態。`constitutive.jsonl`には材料/接触の実分割数、`checkpoint-*.json`には各工程途中の実保存ファイル。
- `*-independent/fields.jsonl`・`result.json`：入力のみから独立新規求解した全受理場と分割数。本体の保存状態や受理列を解の入力にしない。
- 各経路の`fields.jsonl`：全節点変位/反力、水荷重、全材料GPの応力・工学ひずみ・塑性ひずみ・相当塑性ひずみ・塑性仕事・間隙比・水圧、全接触GPのtraction/jump/塑性すべり/累積すべり/仕事/散逸/開閉・カテゴリ。
- `measurements.json`：全23場の全受理点比較、工程別/全成分の最大絶対差・二乗和/L2、材料/接触を分けた感度、同一分割での局所診断。
- `restarts/verification.json`：7地点の実PID・checkpoint SHA・全後続状態数・Run完全一致。
- `completion.json`：計12新規経路と監査の完了。`C03_numeric_acceptance=false`のまま。
- `failure-*.json`・`matrix-failure.json`・各log：失敗を保持。部分出力を新規完全解として扱わない。

`Model`は小型の固定静水圧CONTINUOUS・MC・全活動・境界係数サイクルに限定する。MCC/IDH・過渡水理・構造・BIRTH等を入力しても自動で別モデルへ置換しない。

本体のGP記述子・診断ひずみ増分・固定静水圧の観測量の由来を追加監査する場合：

```powershell
python "$Study/audit_native_descriptors.py" --output $Run --report work/NEW-C03-native-observable-provenance.json
if ($LASTEXITCODE -ne 0) { throw 'native observable provenance mismatch' }
```

これは保存された新規求解出力の監査であり、独立3Dを再度求解したものではない。

## 再開だけをもう一度実行

```powershell
python "$Study/study.py" --mode restart --output $Run --source baseline-native --folder restarts-SECOND --steps 32 --timeout 900
if ($LASTEXITCODE -ne 0) { throw 'real-file restart mismatch' }
```

基準の7個のcheckpointを別PIDで復旧し、最終工程まで解く。後続状態数は208/176/144/112/80/48/16。全後続状態と全Runが`elapsed_seconds`以外完全一致することが合格条件。許容差を設けて通さない。基準と同じ実行版SHAが必要で、再ビルドした別バイナリへcheckpointを付け替えない。

## PC：大きい材料GP数と長い系列

クラウドでは正式行列を12 Tet4・全体32に限定した。12 Tet10は54節点・324材料GP・32接触GPとなり、独立塑性積分と数値接線の負荷が増える。まずこれを別出力へ解く。

```powershell
python "$Study/run_all.py" --output work/NEW-PC-C03-tet10-g32 --solver target/release/examples/wet_joint_plastic_probe.exe --steps 32 --element tet10 --timeout 14400
if ($LASTEXITCODE -ne 0) { throw 'Tet10 not accepted; retain partial data' }
```

RAM 32 GB、空き20 GB（旧証拠を除く）を計画する。Tet10の新しい全系列は未実測で、1経路0.5〜4時間・全行列は半日〜数日を予算として確保する。これは実測値/完了保証ではない。元の大きいopen経路を試したTet10の切戻し記録は残す。収束しなければその入力・受理列・失敗を保存し、新しい増分や許容値へ暗黙変更しない。

同じ小型Tet4で長い全体系列を測る場合は次をPCで実行する。**各出力内では全体刻みを固定して内分割だけを比較する。** C-01の受入やC-02の空間次数に読み替えない。

```powershell
foreach ($n in 64,128,256) {
  python "$Study/run_all.py" --output "work/NEW-PC-C03-tet4-g$n" --solver target/release/examples/wet_joint_plastic_probe.exe --steps $n --element tet4 --timeout 7200
  if ($LASTEXITCODE -ne 0) { throw "g$n incomplete; keep evidence" }
}
```

各行列20分〜2時間程度を初期予算とし、checkpoint/Runの履歴長による増加も見込む（未実測）。同じ条件でTet10の64/128へ進む場合は`--element tet10`と別の出力名、上の14400秒予算を使う。メモリ/時間予算を先に決める。96/768以上のdense独立系・実CAD・GUI/GPUは本ドライバの小型検証範囲を超え、入力/独立実装・資源設計が別途必要で未検証。利用できない大規模独立求解を成功した手順として記載しない。

## 合否の区別

統合観測の完了には、全7工程、全受理列`i/n`（切戻し0）、全場欠落0/有限値、非ゼロ塑性GP・相当塑性ひずみ/塑性仕事、MC/接触法則、時間不変の入力水圧/水荷重、既存の本体力1e-8・線形1e-10・独立力1e-10、実checkpointの完全一致を要求する。これらは`run_all.py`と`measure.py`が検査する。元の7端点変位0.005・traction0.01のゲートは`conclusion.json`に別欄であり、全場の新しい許容値ではない。

**C-03の数値判定は未達**。塑性すべりの接触内分割感度は本体3.88396e-5→1.11421e-4、独立4.51397e-5→7.90003e-5で再増加する。非単調・大きい仕事差・未定義の全場受入基準が残れば完了扱いにしない。数値散逸は物理摩擦散逸と分ける。許容値、旧参照、物理、既定動作、NCの設定を変えて通さない。

## 保存証拠の再計測（新規求解とは別）

```powershell
python "$Study/archive.py" --mode replay --archive $Study --output work/NEW-C03-sealed-replay
if ($LASTEXITCODE -ne 0) { throw 'sealed evidence replay differs' }
```

全4part・ZIP・全1,363エントリーをSHA照合し、新規展開した封印ソースを隔離プロセスで実行して測定JSONのオブジェクト完全一致を確認する。必要ディスク約1 GB、数十秒〜数分を計画する。WindowsのBLAS/丸め差による不一致も失敗として残し、閾値を足さない。Linuxの厳密再現は同じ固定依存版のWSL/Linuxを使う。保存済みデータの再監査を、入力からの新規求解や実checkpoint再開へ数えない。

## 全件ソフトウェア検査

```powershell
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo test --locked --release
cargo build --locked --release
python tools/run_ci_python_tests.py --group python-core-0 --output work/NEW-C03-check/core0
python tools/run_ci_python_tests.py --group python-core-1 --output work/NEW-C03-check/core1
python tools/run_ci_python_tests.py --group python-mean-evidence --output work/NEW-C03-check/mean
python tools/run_ci_python_tests.py --group python-spatial-evidence --output work/NEW-C03-check/spatial
```

4群の発見ID/実行IDの和集合が重複・欠落なく一致し、fail/error/ignored/skipped/filtered-out/expected-failure/unexpected-successが0であることを確認する。旧NC証拠群も全件実行し、CIの時間切れ・キャンセル・未実行を成功にしない。新先端Windows CIはLinux成功や旧CIから推定しない。
