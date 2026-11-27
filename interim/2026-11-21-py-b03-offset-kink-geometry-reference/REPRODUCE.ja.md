# B-03独立幾何参照の再現とユーザーPCでの未確認事項

依頼の保存先ラベルは2026-11-21、実施日は2026-11-12（日本時間）。B-03・NC-01・P6は未達のまま。以下の終了コード2は診断完了・受入未決定、終了1は拒否/エラーである。時間切れ・キャンセル・途中出力を成功にしない。

## Windowsの独立参照（本体を呼ばない）

既存作業と分離したチェックアウトのリポジトリルートで実行する。Python 3.13、固定依存を使用。仮想環境と出力を追跡済みファイルの外へ置く。保存先が既存の場合は新しい名前に替える。

```powershell
git fetch origin codex/py-b03-offset-kink-geometry-reference
git worktree add ..\soil-b03-review origin/codex/py-b03-offset-kink-geometry-reference
Set-Location ..\soil-b03-review
$D = 'interim/2026-11-21-py-b03-offset-kink-geometry-reference'
$V = Join-Path $env:TEMP ('soil-b03-venv-' + (Get-Date -Format yyyyMMdd-HHmmss))
py -3.13 -m venv $V
$Py = Join-Path $V 'Scripts/python.exe'
& $Py -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw 'dependency installation failed' }
$env:OPENBLAS_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
$Out = 'runs/b03-pc-' + (Get-Date -Format yyyyMMdd-HHmmss)
New-Item -ItemType Directory $Out | Out-Null
& $Py "$D/run_reference.py" "$Out/reference" --repository-root .
if ($LASTEXITCODE -ne 2) { throw 'independent diagnostic did not finish normally' }
& $Py "$D/run_checks.py" "$Out/checks"
if ($LASTEXITCODE -ne 0) { throw 'new checks failed' }
& $Py "$D/geometry_reference.py" compare "$Out/reference/synthetic-reference.json" "$Out/reference/synthetic-reference.csv" "$Out/reference/csv-roundtrip.json"
if ($LASTEXITCODE -ne 2) { throw 'synthetic roundtrip rejected' }
Get-Content "$Out/reference/summary.json"
Get-Content "$Out/checks/new-checks.json"
git diff --exit-code f377e8eda41777c84c3bb849af556bf61683b4cb -- . ":(exclude)$D"
if ($LASTEXITCODE -ne 0) { throw 'existing tracked files changed' }
```

参照294条件は通常1〜10秒、新規42チェックは通常1〜10秒を想定（Windows未実測、依存導入時間を除く）。出力はseries、summary、legacy-refusals、人工CSV/JSON、人工往復診断、新規試験のログ・ID・件数。生成先の小計は数MB程度。重い力学・殻・水理求解は含まない。

判定：全294条件と範囲外ケースを欠落なく保存。全候補の `accepted=null`、`b03_complete=false`。原入力/manifestのSHA不一致は拒否。新規試験は42件、failure/error/skip/expected failure/unexpected success=0を要求する。rankの丸め規約を物理受入へ転用しない。数値照合誤差や候補の幾何誤差から新しい許容値を作らない。`synthetic_roundtrip_consistent=true` は人工交換契約だけで、本体の合格ではない。

## Linuxで同じ参照を再現

既存の固定依存が入った `.venv/bin/python` を使う場合：

```bash
D=interim/2026-11-21-py-b03-offset-kink-geometry-reference
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
.venv/bin/python "$D/run_reference.py" runs/b03-reproduce/reference --repository-root .
# 直前の終了コードは2でなければならない。
.venv/bin/python "$D/run_checks.py" runs/b03-reproduce/checks
# 終了0、全件・skip等0をJSONで確認。
.venv/bin/python "$D/audit_existing.py" . runs/b03-reproduce/existing-files.json
# 終了0、existing_tree_preserved_and_available_bytes_identical=true（疎チェックアウトの省略一覧は再検証した内容に数えない）。
```

## 既存Python全件の回帰（独立参照とは別）

現行CIの4群入口を変更せず使う。既存試験には小さい本体CLI求解が含まれる。本独立参照の実装は本体を呼ばない。過去アーカイブを省略したチェックアウトでは全件検証できない。元の追跡ファイルを取得し、欠損をskipや参照値変更で回避しない。疎チェックアウトで不足する場合はPC上で `git sparse-checkout disable` により元の追跡ファイルを展開する。

```powershell
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
if ($LASTEXITCODE -ne 0) { throw 'toolchain setup failed' }
cargo build --locked --release
if ($LASTEXITCODE -ne 0) { throw 'CLI build failed' }
cargo build --locked --release --examples
if ($LASTEXITCODE -ne 0) { throw 'examples build failed' }
foreach ($Group in @('python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence')) {
    & $Py tools/run_ci_python_tests.py --group $Group --output "$Out/existing/$Group"
    if ($LASTEXITCODE -ne 0) { throw "existing test group failed: $Group" }
}
```

Linuxでは `.venv/bin/python` と同じ引数を使用する。元の全534件（tests532＋tools2）、割当265/261/4/4の全IDがちょうど1回ずつ実行され、各群終了0、failure/error/skip/expected failure/unexpected success=0であることを `python-test-inventory.json` から照合する。新規42件を足した576件であり、集中実行や途中失敗を重複加算しない。既存証拠2群は保存されたNC証拠の監査で、新たなNC全撤去求解/合格ではない。

Windows未実測の目安：build 5〜30分、全Python 20〜90分、指定起点の追跡データは約21.7 GB。Gitオブジェクト・作業コピー・封印証拠の展開を含め60 GB以上の空き容量を目安とする（Windows未実測）。PC性能とディスクに依存する。長い回帰はユーザーPCで実行し、途中停止なら未完了と記録する。クラウドで行った今回の実測・件数はREADME末尾とresults/validationを参照。Rust全試験・Windows CIの成功をこの手順から推定しない。

## 将来の本体出力との比較（機能実装後の手順）

現行本体に折れ目/段差の受理や専用CSV/JSON出力があるとは主張しない。実装前に拒否を外してこの手順だけを通さない。信頼できる元入力から、[CONTRACT.ja.md](CONTRACT.ja.md) の全ID・単位・座標・節点順法線・符号付きe・t・中立点/付着点の参照を作って固定する。

```powershell
# パスはユーザーが保存した未加工の元入力参照と将来の本体出力。
$Trusted = Read-Host '固定済み元入力参照JSONのフルパス'
$Body = Read-Host '未加工の本体出力CSVまたはJSONのフルパス'
Get-FileHash -Algorithm SHA256 $Trusted,$Body | Format-List | Out-File "$Out/input-output-sha256.txt"
& $Py "$D/geometry_reference.py" compare $Trusted $Body "$Out/body-comparison.json"
$CompareExit = $LASTEXITCODE
Get-Content "$Out/body-comparison.json"
if ($CompareExit -eq 1) { Write-Host '拒否：ID・単位・基底・符号・e範囲・座標差を確認。補正して再提出しない。' }
if ($CompareExit -eq 2) { Write-Host '交換契約の整合のみ。許容値未定義、accepted=null。' }
if ($CompareExit -notin @(1,2)) { throw 'unexpected comparison exit' }
```

付着点だけ一致する法線/eの同時反転、平均化、範囲内へのクリップも元入力との差として拒否する。丸め差も今は受入許容値がないため拒否する。許容値の上書き引数はない。バイナリSHA、元の入力・全出力、ログ、使用基底と向きの記録を保持し、出力から参照を逆作成しない。

## コンパクト面宣言・Native準備・表示・GUI/GPU・実CAD（PC専用、今回未確認）

この節は実装/表示の合格ではなくユーザーPCで確認する手順。現行の新方式拒否を保持する。ソースの `prepare_native_plates.py` は新方式/非ゼロeを明示拒否し、既存の平面・直接CLI入力の限定受入とは別である。

1. 実CADとNativeのgeometry-only inventoryを別フォルダーへコピーし、元ファイルと面ID・共有節点ID・CAD面の向き/単位・段差位置・折れ目角・厚さ・符号付きeを記録する。実入力がない場合は「未確認」のまま。外見が合うように面の平坦化、節点複製、平均法線、符号反転、eクリップをしない。
2. 次のコマンドで**実際に保存した** `prepared-model.json`、mesh/stages/native-captureとSHAを持つ入力を検査する。コンパクト宣言は `plate_surfaces`。新方式を含む拒否例の原入力を保ち、準備先は存在しない別フォルダーにする。

```powershell
$Source = Read-Host '保存したgeometry-only Native inventoryフォルダー'
Get-ChildItem $Source -File | Get-FileHash -Algorithm SHA256 | Export-Csv "$Out/native-source-sha256.csv" -NoTypeInformation
& $Py tools/prepare_native_plates.py $Source "$Out/native-prepared" *> "$Out/native-prepare.log"
$NativeExit = $LASTEXITCODE
Get-Content "$Out/native-prepare.log"
Write-Host "Native準備の終了コード=$NativeExit"
# 現行の新方式/非ゼロeは非ゼロ終了と明示拒否を期待。出力欠落を表示成功にしない。
Get-CimInstance Win32_VideoController | Select-Object Name,DriverVersion,AdapterRAM | ConvertTo-Json | Out-File "$Out/gpu.json"
```

3. 現行は拒否を確認するだけで、候補の表示受入は未確認。将来Native側の意味論を別途実装・検証した時点で、実際のアプリを起動し（`$NativeExe=Read-Host '使用中Nativeアプリ実行ファイル'; Start-Process -FilePath $NativeExe`）、UIの通常の開く操作で保存モデルを選ぶ。架空のCLI引数でロードしない。
4. 付着面と中立面を別表示し、全共有ID・各面法線・±e・t・局所軸を数値出力と照合する。0/90°近傍/反平行近傍、正負段差、e=0/±t/2/±0.6t、節点順反転、単位違い、ID欠損を別ジョブにする。範囲外と向き不一致は拒否のまま。スクリーンショット、未加工CSV/JSON、モデル/バイナリSHA、GPU/driverを保存する。表示誤差許容値は未定義で `accepted=null`。
5. 設置/表示更新/Undo/再読込で面の平均化、古いIDの残留、符号逆転、表示だけの補正が起きないかを追う。実CAD別モデル・別GPUはそれぞれ未確認から始める。元の限定平面CLI合格をNativeの合格へ転用しない。

想定時間（すべてWindows未実測）：拒否/準備確認1〜10分、実GUI手動照合30〜120分/モデル、実CAD整備・複数GPU確認は数時間以上。今回この節をクラウドで実行しない。大型求解の計算時間・収束基準はモデル・仕様が未定義のため未定義であり、架空の成功基準や見込み時間を作らない。PCで時間切れ・キャンセルになれば未確認/未完了とし、CI成功へ数えない。
