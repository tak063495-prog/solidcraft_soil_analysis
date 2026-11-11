# Windowsでの再現とPC専用の追加検証

PowerShell、Python 3.13、Rust 1.99.1を使用する。新しい短い作業パスを使い、
他のPR #13作業ツリーのブランチを切り替えない。以下の時間はPC向けの計画値で、
Windowsでの実測ではない。全コマンドは終了コードを確認し、失敗・キャンセル・時間切れの出力を残す。

## 小型全試験

```powershell
git clone --branch codex/acceptance-tests-todo --single-branch https://github.com/tak063495-prog/solidcraft_soil_analysis.git C:\soil-acceptance
if ($LASTEXITCODE -ne 0) { throw 'clone failed' }
Set-Location C:\soil-acceptance
git merge-base --is-ancestor 5f564621c7396c2a605edbfe5451cd800cf318d3 HEAD
if ($LASTEXITCODE -ne 0) { throw 'wrong base commit' }
git rev-parse HEAD
py -3.13 -m venv .venv-acceptance
if ($LASTEXITCODE -ne 0) { throw 'venv failed' }
$py = (Resolve-Path .\.venv-acceptance\Scripts\python.exe).Path
& $py -m pip install -r requirements-reference.lock.txt pytest==9.0.2
if ($LASTEXITCODE -ne 0) { throw 'dependency install failed' }
$env:OPENBLAS_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:CARGO_BUILD_JOBS='2'
$env:RUST_TEST_THREADS='2'
rustup run 1.99.1 cargo build --locked --release --bins --examples
if ($LASTEXITCODE -ne 0) { throw 'build failed' }
$out = 'runs/acceptance-todo-' + (Get-Date -Format 'yyyyMMdd-HHmmss')
if (Test-Path $out) { throw 'output must be new' }
New-Item -ItemType Directory $out | Out-Null
git rev-parse HEAD | Out-File "$out/head.txt"
& $py -m pip freeze | Out-File "$out/python-packages.txt"
Get-FileHash target/release/soil-fem-3d.exe | Format-List | Out-File "$out/solver-sha256.txt"
& $py -m pytest tests/acceptance_todo -ra --junitxml="$out/acceptance.xml" --basetemp="$out/pytest-data" 2>&1 | Tee-Object "$out/acceptance.log"
if ($LASTEXITCODE -ne 0) { throw 'acceptance tests failed, including strict XPASS' }
& $py -m unittest discover -s tests -p 'test_*.py' -v 2>&1 | Tee-Object "$out/existing-tests.log"
if ($LASTEXITCODE -ne 0) { throw 'existing tests failed' }
& $py -m unittest discover -s tools -p 'test_*.py' -v 2>&1 | Tee-Object "$out/existing-tools.log"
if ($LASTEXITCODE -ne 0) { throw 'existing tools tests failed' }
```

想定：依存取得/ビルド5〜30分、新規小型試験1〜5分、既存Python全件5〜60分
（CPU・既存ビルドキャッシュに依存）。ソース/過去の証拠のcloneは大容量なので通信時間を別に見積もる。
想定件数はREADMEと検証JSONを参照し、Windows固有条件の差はIDとともに記録する。
新規試験の本来のskipは0。JUnitではxfailも `skipped type="pytest.xfail"` として保存されるため、
XMLのskipped合計を通常skipの件数と混同しない。
`--runxfail`で全項目を強制合格させたり、xfail数をpassへ加算したりしない。

各CLI小型試験の入力・stdout/stderr・終了コード・checkpoint・再開/改変データは
`pytest-data`以下へ残る。各CLIの120秒制限を超えると通常のFAILになり、xfailへ変換されない。
既存試験の収集設定もCIの件数条件も変更しない。試験用ディレクトリーへ `__init__.py` を追加しない。

## 現行の限定受入の再現（新しい組合せの成功ではない）

新しい出力パスで実行する。ビルド後、想定各1〜15分。

```powershell
& $py tools/verify_joint_ext_combos.py --solver target/release/soil-fem-3d.exe --rehash target/release/examples/joint_ext_checkpoint_rehash.exe --output "$out/joint-ext"
if ($LASTEXITCODE -ne 0) { throw 'limited JOINT study failed' }
& $py tools/plate6_offset_study.py "$out/offset-plus" --binary target/release/soil-fem-3d.exe --offset 0.05
if ($LASTEXITCODE -ne 0) { throw 'positive offset study failed' }
& $py tools/plate6_offset_study.py "$out/offset-minus" --binary target/release/soil-fem-3d.exe --offset -0.05
if ($LASTEXITCODE -ne 0) { throw 'negative offset study failed' }
```

出力：入力、全Run、別CLI再開、独立測定、改変拒否ログ、VTU等。
合否：(a)/(b)の元の全節点1e-10・非線形1e-8・線形1e-10・独立弾性相対1e-7・
構造1e-9・板生値絶対1e-8/相対1e-7・拘束1e-10 m・仕事1e-10を維持する。
JOINTドライバの `requested_all_combinations_passed=false`、板の `p6_complete=false` を保持する。
陽性対照が通っても湿潤先端/曲面/Richards湿潤JOINT/NC/P6は完了にならない。

## 大規模メッシュ・長系列（このクラウド作業では未実行）

以下は、将来の実装と固定済み入力/独立参照がそろった場合のPC側コマンドである。
**現時点ではB-01の古典殻h系列用入力と許容差、A-04の双方向/開閉/世代移送の仕様がないため、
実行可能な受入ベンチマーク一式があるとは扱わない。** 次の存在検査で不足を明示的に止める。
ガードを外す、範囲外eを小さくする、許容値を新設/緩和する作業はこの手順に含めない。

各研究のPC入力を `C:\soil-acceptance-inputs\<study>\h0,h1,h2` に固定し、各水準に
`mesh.json`、`ground.json`、`case.json`、元参照/出力量/許容差/単位を記した `reference-contract.json` を用意する。
製造解・h列の水準間の対応点も固定し、refファイルとsolverのSHA256を出力へ保存する。
研究名は `wet-tip`、`darcy-birth`、`richards-joint`、`scordelis-lo`、`pinched-cylinder`、
`hemisphere`、`offset-combination`。現在のリポジトリーにはこの新規入力セットを捏造して同梱しない。

```powershell
$study = 'scordelis-lo'  # 固定契約が準備された研究だけ選ぶ
$inputRoot = "C:\soil-acceptance-inputs\$study"
foreach ($h in @('h0','h1','h2')) {
    $src = Join-Path $inputRoot $h
    foreach ($f in @('mesh.json','ground.json','case.json','reference-contract.json')) {
        if (!(Test-Path (Join-Path $src $f))) { throw "UNVERIFIED: missing $src/$f" }
    }
}
$dest = Join-Path $out $study
New-Item -ItemType Directory $dest | Out-Null
foreach ($h in @('h0','h1','h2')) {
    $src = Join-Path $inputRoot $h
    $dst = Join-Path $dest $h
    New-Item -ItemType Directory $dst | Out-Null
    Copy-Item "$src/*.json" $dst
    Get-FileHash "$src/*.json" | Format-List | Out-File "$dst/input-sha256.txt"
    $timer = [Diagnostics.Stopwatch]::StartNew()
    & target/release/soil-fem-3d.exe solve-stages "$src/mesh.json" "$src/ground.json" "$src/case.json" "$dst/run.json" 2>&1 | Tee-Object "$dst/solve.log"
    $code = $LASTEXITCODE
    $timer.Stop()
    @{exit_code=$code;elapsed_seconds=$timer.Elapsed.TotalSeconds} | ConvertTo-Json | Out-File "$dst/exit.json"
    if ($code -ne 0) { throw "FAILED: $study/$h exit=$code" }
    $run = Get-Content "$dst/run.json" -Raw | ConvertFrom-Json
    if ($run.complete -ne $true) { throw "INCOMPLETE: $study/$h" }
}
```

`richards-joint`は `solve-stages`の代わりに、ONE_WAY契約なら `solve-one-way` を使う。
双方向コマンド/入力契約は未実装なので提示しない。
時間計画：小型を先に実測し、h0/h1/h2で概ね10分〜2時間／1〜8時間／4〜24時間以上を確保する。
これは保証値や完走記録ではない。時間切れ・キャンセル・未収束はFAILED/INCOMPLETEのまま保存する。
3水準が完走したことだけでは受入できない。元の固定参照・全場・仕事・水収支・h/Δt感度・
実checkpoint再開を照合する必要がある。古典殻の許容差は未定義なので、曲面の参照差が小さく見えても判定は未検証。
既知の元ゲートを持つ物理量にはREADMEの値をそのまま適用し、新しい収束合格値を作らない。

## 実CAD・Native表示

実CADの元データと固定メッシュ/先端辺/所有要素、Nativeの版・GPU/ドライバー・画面記録はPC上で保存する。
VTUとJSONの付着面/中立面、符号付きe、法線、JOINT世代、水頭/水圧、変位/回転を点と時刻で照合する。
想定30分〜数時間（入力整備は別）。比較対象がない、表示仕様が未定義、GUIが未起動なら未検証である。
データ往復だけをGUI/GPUの成功へ読み替えない。NC-01、P6、却下済みNC調整はこの手順の実行対象に追加しない。
