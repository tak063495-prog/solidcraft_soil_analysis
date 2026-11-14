# A-03限定診断の再現

保存先の日付は依頼ラベル。PR #13の指定起点
`3b1133d3be9905dc24f62c27b5b1eb456a83f2e2` の子孫で再現する。
PRや参照draft PR #16/#19/#21/#24/#31をマージしない。
既存出力を消さず、毎回未存在の出力フォルダーを指定する。
対象・符号・収支の未定義は [SPEC.ja.md](SPEC.ja.md) を参照。

## Windows PCの準備・全件検証

リポジトリルートで実行する。作業中の変更があれば別worktreeで再現する。

```powershell
git status --short
git fetch origin codex/plastic-spatial-convergence
git switch --detach origin/codex/plastic-spatial-convergence
git rev-parse HEAD
git merge-base --is-ancestor 3b1133d3be9905dc24f62c27b5b1eb456a83f2e2 HEAD
if ($LASTEXITCODE -ne 0) { throw '指定起点の子孫ではありません' }
$Py = '.venv/Scripts/python.exe'
if (-not (Test-Path $Py)) { py -3.13 -m venv .venv }
if (-not (Test-Path $Py)) { throw 'Python 3.13のvenvを作成できませんでした' }
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
& $Py -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw '固定依存の配置失敗' }
cargo fmt --all --check
if ($LASTEXITCODE -ne 0) { throw 'fmt失敗' }
cargo clippy --locked --all-targets -- -D warnings
if ($LASTEXITCODE -ne 0) { throw 'clippy失敗' }
cargo test --locked --release
Write-Output "Rust main exit_code=$LASTEXITCODE"
cargo test --locked --release --examples
Write-Output "Rust examples exit_code=$LASTEXITCODE"
cargo build --locked --release --examples
if ($LASTEXITCODE -ne 0) { throw '観測用examplesビルド失敗' }
cargo build --locked --release
if ($LASTEXITCODE -ne 0) { throw '標準CLIビルド失敗' }
foreach ($group in @('python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence')) {
  & $Py tools/run_ci_python_tests.py --group $group --output "runs/NEW-a03-validation/$group"
  Write-Output "$group exit_code=$LASTEXITCODE"
}
```

Rust初回ビルド・全件は10〜60分、Python通常群は数分〜20分、旧封印証拠2群は
追加で数分〜数十分を計画する（Windows未測定の予算目安）。元証拠の読込み・展開に
十数GB以上の空きディスクを用意する。各群のinventoryで全発見IDと割当の和集合、
重複0、未実行数を確認する。過去のLinuxでの旧NC setUpClass失敗をskip・許容差化・
参照値変更・Windows分岐の移植で通さない。欠落アーカイブなら元Git blobと一致する
バイトを配置し、初回失敗を保持して該当する通常群全体を再実行する。
今回の環境では旧NC両群も厳密比較のまま8件成功したが、過去の8件未実行は書き換えず、
実際の今回件数と区別する。環境差の原因は未特定であり、NCの数値受入へ転用しない。

## 小型2ケースを本体・独立の双方で新規求解

```powershell
$Study = 'interim/2026-11-26-wet-soil-birth-darcy-a03/study.py'
& $Py $Study --steps 4 --output runs/NEW-a03-PC-n4
if ($LASTEXITCODE -ne 0) { throw '求解/比較/復旧/拒否が未完了。途中出力を保持' }
& $Py interim/2026-11-26-wet-soil-birth-darcy-a03/water_cache_audit.py --run runs/NEW-a03-PC-n4
if ($LASTEXITCODE -ne 0) { throw 'N=4保存場の水圧cache/浮力監査失敗' }
& $Py interim/2026-11-26-wet-soil-birth-darcy-a03/pointwise_audit.py --run runs/NEW-a03-PC-n4
if ($LASTEXITCODE -ne 0) { throw '各受理点/準備点の既存ゲート検査失敗' }
```

既存形状、静水位1.5 m/2.0 mの2ケース。2.0 mだけ既存の外水面集合exteriorを明示する。各ケースは初期平衡→土BIRTH→JOINT BIRTH/B
→上載圧20 kPa、入力の各4等分、切戻しなし。32受理状態対と2準備状態を比較する。
本体観測あり・観測なし・標準CLI対照と、別Python PIDによる独立求解を分けて記録する。
5〜20分程度のPC計画枠（未実測）。実時間は各process.jsonとログで確認する。

|出力|内容|
|---|---|
|`measurements.json`|A-03全体のpassed=false、各対象の受理/拒否、全場差、復旧、改変|
|`source-fingerprints.json` / `source-snapshot.zip`|新規求解前のソース/SPECとSHA|
|`binaries.json` / `bin/`|固定実行版とSHA。処理中にビルドし直さない|
|各ケースの`native/run.json`, `accepted.jsonl`|全Run、全受理状態、途中実checkpoint|
|`native/constitutive.jsonl`, `preparations.json`|生産公開APIによる確定状態の力分解・面再生成観測。独立求解ではない|
|`independent-fields.jsonl`, `independent-prepared.jsonl`|入力のみから新規求解した独立3D全場|
|`native-fields.jsonl`, `native-prepared.jsonl`|本体の対応全場。ひずみ等の派生観測を含む|
|`comparison.json`, `prepared-comparison.json`|全成分の最大絶対差・相対L2・工程別差・カテゴリ不一致|
|`restarts/`|実ファイル別PID復旧、全後続状態と全Run|
|`tamper/`|外側hash再計算後の同一実行版による改変拒否|
|`refusals/`|CLIの終了コード/理由、生入力、Darcy完了prefix・MCC対照|

合格条件を分離する。全状態が完走し、入力のi/N格子、全ID/座標、全場の欠落なし、
カテゴリ一致、別PIDで全後続状態とRun全体がelapsed_secondsだけ除いて厳密一致、
狙ったガードで拒否され結果ファイルなしを検査する。
既存の独立弾性相対1e-7、節点力1e-10 kN、拘束1e-10 m、線形残差1e-10を緩めない。
全場の包括許容値は未定義のまま。`numerical_acceptance=false`、A-03 `passed=false`を
消したり、水量のnullを0へ変えて水収支合格にしたりしない。
時間切れ・未収束・キャンセル・旧CIを新CI成功へ転用しない。

## PC専用の長い系列

クラウドでは4等分の小型だけを実行する。長い系列は次のコマンドで実行する。
メッシュ・物理・材料・許可範囲を変えず、工程の格子だけを細かくする。

```powershell
$env:A03_PROCESS_BUDGET_SECONDS='21600'
try {
  foreach ($n in @(16,32,64)) {
    & $Py $Study --steps $n --output "runs/NEW-a03-PC-n$n"
    if ($LASTEXITCODE -ne 0) { throw "n=$n が未完了。生結果・ログを保持" }
  }
} finally { Remove-Item Env:A03_PROCESS_BUDGET_SECONDS -ErrorAction SilentlyContinue }
```

出力は上表と同じで、各Nの受理状態対は8N（128/256/512）。補助water_cache_auditは配布N=4用の監査であり、この長系列には呼ばない。
Study自身が全格子・全場を比較し、
共通進行度の差と工程端点を別に確認する。時間の目安はN16で20〜60分、N32で1〜3時間、
N64で3〜8時間（Windows未測定の計画枠）。毎受理点のseal/resume認証があるため
単純なN比例ではない。1プロセス6時間で停止した場合も失敗・未達のまま保持する。

独立密行列参照の350自由度予算を外して大きなメッシュを流す手順ではない。
一般3Dの空間収束・新設MCC・旧履歴・Darcy BIRTHには未定義仕様が残る。
現行拒否を解除してPCで「試す」ことはしない。実CAD/GUI/GPUや大規模の新たな受入は
今回の検証外であり、この固定メッシュの長系列をその代用としない。

## Linuxで同じ検証

```bash
source /workspace/soil-environment/activate.sh
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo test --locked --release
cargo test --locked --release --examples
cargo build --locked --release --examples
cargo build --locked --release
python interim/2026-11-26-wet-soil-birth-darcy-a03/study.py --steps 4 --output runs/NEW-a03-linux
for group in python-core-0 python-core-1 python-mean-evidence python-spatial-evidence; do
  python tools/run_ci_python_tests.py --group "$group" --output "runs/NEW-a03-tests/$group"
  echo "$group exit_code=$?"
done
```

## 保存証拠の照合

```powershell
& $Py interim/2026-11-26-wet-soil-birth-darcy-a03/package_evidence.py --verify interim/2026-11-26-wet-soil-birth-darcy-a03/evidence-manifest.json --scratch runs/NEW-a03-archive-check
if ($LASTEXITCODE -ne 0) { throw 'part/結合ZIP/全エントリーのSHA不一致' }
```

これは保存証拠の整合監査。上記の新規求解や実ファイル復旧の代わりにはならない。
