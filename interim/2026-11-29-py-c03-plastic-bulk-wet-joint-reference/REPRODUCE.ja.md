# 再現手順：小型独立参照とWindows PCでの後続検証

C-03は未達。正常終了は計算・形式検査の完了を意味し、物理的な数値受入とは別。全履歴の許容値未定義の項目は `accepted=null` を維持する。実施2026-11-12（日本時間）、フォルダーの日付はラベル。

## Windows：独立参照だけを実行

リポジトリルートからPowerShellで実行。Python3.13、NumPy/SciPyの既存固定依存を使う。新規の参照・診断・比較器はRustや既存Pythonソルバを起動しない。出力先は未使用名を選び、失敗出力を消さない。

```powershell
git fetch origin codex/py-c03-plastic-bulk-wet-joint-reference
git switch codex/py-c03-plastic-bulk-wet-joint-reference
git log -1 --format='%H'
py -3.13 -m venv work/c03-venv
& work/c03-venv/Scripts/Activate.ps1
python -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw 'dependency installation failed' }
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
$Study='interim/2026-11-29-py-c03-plastic-bulk-wet-joint-reference'
$Run='work/NEW-py-c03-small'
python "$Study/checks.py" -v
if ($LASTEXITCODE -ne 0) { throw 'independent self-check failed' }
python "$Study/evidence.py" --out "$Run/source-evidence.json"
if ($LASTEXITCODE -ne 0) { throw 'original evidence missing' }
python "$Study/diagnose.py" --out $Run --levels 8 16 32 64
if ($LASTEXITCODE -ne 0) { throw 'diagnostics incomplete' }
python "$Study/export_reference.py" --history "$Run/column-history.json" --out "$Run/consumer-example"
if ($LASTEXITCODE -ne 0) { throw 'reference export failed' }
python "$Study/compare.py" --manifest "$Run/consumer-example/manifest.json" "$Run/consumer-example/reference.json" "$Run/consumer-example/reference.csv" --out "$Run/json-csv-difference.json"
if ($LASTEXITCODE -ne 0) { throw 'JSON/CSV strict comparison failed' }
```

原資料のSHA検査は指定Git blobのLFバイト列を要求する。WindowsでCRLF変換されたcheckoutは失敗として扱う。既存ファイルを書き換えて通さず、必要なら `core.autocrlf=false` の別cloneを使う。

小型計算のLinux実測は約5秒、17自己試験は約1秒。Windowsで10秒〜1分を初期目安とする（依存取得を除く、未実測）。CPU2コア・RAM2 GB程度。`diagnostics.json` は次数・接線・仕事・力・収支差、`stress-point-history.json` は最細の全規定ひずみと両方式の応力/塑性ひずみ履歴、`column-history.json` は6条件×24状態の全履歴、`source-evidence.json` は原資料の出典pointer・SHA、`consumer-example/` は明示manifestとCSV/JSON。独立小型出力は1 MB未満を目安とする。

自己試験は17件すべて成功、失敗/error/skip/expected failureは0が形式上の成功条件。アルゴリズムの丸め・反復停止検査値はchecks.pyで明示し、C-03の許容値へ転用しない。原資料の許容値はREADMEの行別表を参照する。応力差・次数・塑性仕事・接線差・全体エネルギーの新しい物理ゲートは未定義。0差の自己比較も `accepted=null`。

## 将来の本体CSV/JSONを厳密比較

新規ファイル `export_reference.py` と `compare.py` が入力仕様の完全な例。manifestは参照側が事前に宣言した全工程、各工程の増分1..N、増加する物理時刻、全bulk/contact IDを持つ。例示模型では無次元cycleを**24秒**に対応させることを明記している。元3Dの工程時間を24秒と推定しない。

本体出力はこのPRでは追加していない。次は同じ物理入力・同じ離散問題の参照と本体出力が揃った後のコマンドであり、現在実施済みの比較とはしない。

```powershell
# manifest・reference・nativeはユーザーが確認済みの同一問題のもの。
python "$Study/compare.py" --manifest work/READY/manifest.json work/READY/reference.json work/READY/native.csv --out work/NEW-native-diff.json
if ($LASTEXITCODE -ne 0) { throw 'metadata/ID/increment/time mismatch; do not interpolate' }
```

JSONは `schema, metadata, records`、CSVは `kind,stage,increment,time,id,pore` と `stress_0..5,plastic_strain_0..5,traction_0..2,plastic_slip_0..2`。CSVには隣接する `native.csv.meta.json` の `schema,metadata` が必須。kind別の不要列は空欄、必要列は有限数。JSONの重複キー、CSVの重複ヘッダーも拒否する。

metadataの単位・符号・成分順は `compare.META` の全キーを省略せず指定する。応力はkPa、変位/すべりm、時間s、テンソルせん断、xx/yy/zz/xy/yz/zx、x/y/z。`stress_kind` はeffectiveまたはtotal、`biot_alpha` は明示数値。pore_modelはfixed_hydrostatic_nonnegativeで負のpを拒否する。接触tractionはeffective_skeleton、jumpはu_b−u_aに共役なa側力/面積。元の圧縮正応力・塑性ひずみや工学せん断は、供給側で出典を残した明示変換が必要。比較器は黙示変換しない。全応力を有効応力とラベルする誤申告は入力作成側で検査する必要がある。

未知のID・欠損増分・時刻相違があれば終了コード非0。全/有効応力やpの符号相違も拒否する。補間・最近傍GP対応・端点だけへの縮約はない。差は場ごとの最大絶対値、絶対L2、参照非ゼロの場合の相対L2（参照0ならnull）。元R:48の1e-12分母床を新しい受入へ持ち込まない。許容値が未定義なので比較の正常終了も合格認定ではない。

## PC：局所参照の長い分割系列

クラウドでは8..64分割と24増分の小型問題のみ実行した。Windows CPU4コア/RAM8 GBで次を別出力へ実行できる。1〜5分を初期予算とする（未実測）。

```powershell
python "$Study/diagnose.py" --out work/NEW-PC-py-c03-fine --levels 64 128 256 512
if ($LASTEXITCODE -ne 0) { throw 'fine reference did not complete' }
```

すべて正の2倍系列が必要。次数を出すために誤差を0へ丸めたり、負の仕事・降伏driftを削除しない。最後の分割でも仕事・全体収支差が残る可能性があり、未定義ゲートを補って成功扱いにしない。

## PC：既存12 Tet4/Tet10の長い全履歴と実再開

ここからは**既存追補ドライバを使う別の本体＋既存独立3D作業**。本追加のNumPy応力点・柱参照とは分ける。今回クラウドでは実行していない。入力条件・実分割数・GP位置・工程・全増分を保持したまま将来比較に必要な出力が揃うか、まず確認する。出力を揃えるための本体改変は本PRの対象外。

Windows11、Python3.13、Rust1.99.1、RAM32 GB、空き80 GB程度。既存の全履歴/再開出力と旧証拠を含めた予算。過去の未使用postfold保存partのcheckout除外は既存CIと[元REPRODUCE](../2026-11-20-wet-joint-plastic-bulk-integration/REPRODUCE.ja.md)と同じ設定を使えるが、必要な証拠や試験は除外しない。

```powershell
rustup toolchain install 1.99.1 --profile minimal
$env:CARGO_BUILD_JOBS='2'
$env:RUST_TEST_THREADS='2'
cargo build --locked --release --example wet_joint_plastic_probe
if ($LASTEXITCODE -ne 0) { throw 'native probe build failed' }
$Old='interim/2026-11-20-wet-joint-plastic-bulk-integration'
foreach ($n in 64,128,256) {
  python "$Old/run_all.py" --output "work/NEW-PC-C03-tet4-g$n" --solver target/release/examples/wet_joint_plastic_probe.exe --steps $n --element tet4 --timeout 7200
  if ($LASTEXITCODE -ne 0) { throw "g$n incomplete; retain every partial/failed record" }
}
python "$Old/run_all.py" --output work/NEW-PC-C03-tet10-g32 --solver target/release/examples/wet_joint_plastic_probe.exe --steps 32 --element tet10 --timeout 14400
if ($LASTEXITCODE -ne 0) { throw 'Tet10 incomplete; not accepted' }
```

元手順による未実測予算はTet4各行列20分〜2時間、12 Tet10（54節点、324材料GP、32接触GP）は1経路0.5〜4時間・全行列半日〜数日。ビルド10〜30分と証拠取得時間を別に見込む。96/768以上や任意大規模3Dの実行可能性は未確認で、用意のない系列をこのコマンドで扱えるとはしない。

主出力は `*-native/raw/accepted.jsonl`、材料/接触実分割数、`*-independent/fields.jsonl`、各入力/SHA、checkpoint、`restarts/verification.json`、measurements、completion、失敗/部分ログ。原資料での限定確認条件は全7工程固定列i/n・棄却0・全場欠落0/有限・非ゼロ塑性履歴、本体力1e-8/真の線形1e-10/独立力1e-10、実再開の時間以外完全一致。変位.005/traction.01は7端点だけ。元open=1.0の失敗を.12の成功で置換しない。MCCへの自動置換は行わない。

再開だけの具体例（上で完成したg64を使い、同一実行版SHAを要求）：

```powershell
python "$Old/study.py" --mode restart --output work/NEW-PC-C03-tet4-g64 --source baseline-native --folder restarts-SECOND --steps 64 --timeout 7200
if ($LASTEXITCODE -ne 0) { throw 'separate-process checkpoint resume failed' }
```

保存JSONの再監査を新規求解/実再開と数えない。プロセス終了・最後の工程までの後続全状態・同一バイナリ・時間以外のRun完全一致を確認する。時間切れ・キャンセル・未収束・固定列からの切戻しは未達として保持する。いずれもC-03全体の完了とはしない。

## Linux：今回の新規試験と既存Python全件

既存収集設定と4群の割当ては変更しない。環境の既存有効化（PATHを含む）後、ルートで以下を実行した。既存試験がCLIを必要とするため、全件前にCLIをrelease buildする。これは独立参照から本体を呼ぶこととは別の既存全件検証。

```bash
source /workspace/soil-environment/activate.sh
python interim/2026-11-29-py-c03-plastic-bulk-wet-joint-reference/checks.py -v
cargo build --locked --release
python tools/run_ci_python_tests.py --group python-core-0 --output work/c03-validation/python-core-0
python tools/run_ci_python_tests.py --group python-core-1 --output work/c03-validation/python-core-1
python tools/run_ci_python_tests.py --group python-mean-evidence --output work/c03-validation/python-mean-evidence
python tools/run_ci_python_tests.py --group python-spatial-evidence --output work/c03-validation/python-spatial-evidence
```

Linuxで4群を別プロセスとして実行する場合も全群の終了を確認し、各群の終了コード・発見ID・assigned_total・tests_run・errors/failures/skipを保存する。旧証拠2群のsetUpClassエラーは成功に数えず、対象8件未実行も別記する。既知エラーをskipへ置換しない。既存のテスト全件成功要件はエラーがあれば未達のまま。Windows新PR先端CIは、実際の終了結果を確認するまでは未確定。

今回の確定件数は `results/validation.json`、生ログと各群inventoryは `results/validation/`。既存ファイルの保持検査は `results/preservation.json`。PRはdraftのままとし、マージしない。

Linux最終件数は、新規17件成功、既存573件発見・通常565件成功・旧NC証拠2群のsetUpClassエラー2件・8件未実行。通常群は300件/265件、証拠群は各4件割当て/実行0件。全群skip0で、既存全件成功要件は未達。CLI build終了0、既存4群終了コードは0/0/1/1。新規試験約1秒、既存通常群約407/526秒、mean/spatial群約742/192秒。40分の作業目安は旧証拠取得を含め超過した。保存証拠の再監査を新しい3D求解へ加算していない。

## draft PR作成の継続用コマンド

push後、APIの接続制限がある環境ではPR作成が別途必要になる。本文は `PR_BODY.ja.md` に具体化してある。APIへ到達でき、対象リポジトリに認証済みのPCからのコマンドは次の通り。既存PRが作成済みなら重複作成しない。

```powershell
gh pr list --head codex/py-c03-plastic-bulk-wet-joint-reference --state all
gh pr create --draft --base codex/plastic-spatial-convergence --head codex/py-c03-plastic-bulk-wet-joint-reference --title 'C-03: 独立Python MC塑性バルク・湿潤摩擦接合参照（未達維持）' --body-file "$Study/PR_BODY.ja.md"
if ($LASTEXITCODE -ne 0) { throw 'draft PR creation failed; no merge' }
```
