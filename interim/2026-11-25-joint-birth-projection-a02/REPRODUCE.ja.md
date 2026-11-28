# A-02の再現（Windows PC / Linux）

保存先の日付は依頼ラベル。実施日は2026-11-12（日本時間）。PR #13をマージしない。
既存出力を削除・上書きせず、毎回未存在の出力フォルダーを指定する。
以下の長い系列、GUI/GPU、実CAD、大規模計算はクラウドで実行していない。
範囲は [SPEC.ja.md](SPEC.ja.md)、実測値は [README.ja.md](README.ja.md) を参照する。

## Windowsの準備と全件検査

リポジトリのルートで実行する。作業中の変更を確認し、必要なら別worktreeを使う。
既存ブランチを保持したままdetached checkoutで再現できる。

```powershell
git status --short
git fetch origin codex/plastic-spatial-convergence
git switch --detach origin/codex/plastic-spatial-convergence
git rev-parse HEAD
git merge-base --is-ancestor 702a19a7ce189f1173caed727ee2a0059298a607 HEAD
if ($LASTEXITCODE -ne 0) { throw '指定起点の子孫ではありません' }
$Py = '.venv/Scripts/python.exe'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
& $Py -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw '固定依存の配置失敗' }
& $Py -c "import sys,numpy,scipy; print(sys.version);print(numpy.__version__,scipy.__version__)"
cargo fmt --all --check
if ($LASTEXITCODE -ne 0) { throw 'fmt失敗' }
cargo clippy --locked --all-targets -- -D warnings
if ($LASTEXITCODE -ne 0) { throw 'clippy失敗' }
cargo test --locked --release
if ($LASTEXITCODE -ne 0) { throw 'Rust試験失敗。ログを保持してください' }
cargo test --locked --release --examples
if ($LASTEXITCODE -ne 0) { throw 'example内のRust試験失敗' }
cargo build --locked --release --examples
if ($LASTEXITCODE -ne 0) { throw '観測用examplesビルド失敗' }
cargo build --locked --release
if ($LASTEXITCODE -ne 0) { throw '本番CLIビルド失敗' }
foreach ($group in @('python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence')) {
  & $Py tools/run_ci_python_tests.py --group $group --output "runs/NEW-a02-validation/$group"
  # 失敗しても他群を実行し、全群の終了コード・inventoryを保存する。
  Write-Output "$group exit_code=$LASTEXITCODE"
}
```

初回release全ターゲットのビルド・Rust全件はPC性能により10〜60分程度を見込む
（PC実測ではない）。Python通常群は数分〜20分、旧証拠2群は元データ約6.4GBを
読み、展開・厳密再計測のため別に数分〜数十分と十分な一時ディスク空間が必要。
全群の `python-test-inventory.json` で発見IDの和集合・重複0・未実行数を確認する。
Linuxの旧NC mean/spatialは既知の厳密差によるsetUpClass失敗を保持する。
Windowsの既存プラットフォーム分岐をLinuxへ移植しない。skip、基準緩和、参照上書きで通さない。

## 同じ小型3経路を新規求解

```powershell
$Study = 'interim/2026-11-25-joint-birth-projection-a02/study.py'
& $Py $Study --steps 4 --output runs/NEW-a02-PC-n4
if ($LASTEXITCODE -ne 0) { throw '求解・比較・復旧・拒否検査の異常。途中出力を保持' }
& $Py interim/2026-11-25-joint-birth-projection-a02/production_cli_controls.py --run runs/NEW-a02-PC-n4 --output runs/NEW-a02-PC-standard-cli
if ($LASTEXITCODE -ne 0) { throw '標準CLIの追加対照失敗' }
foreach ($script in @('geometry_audit.py','projection_diagnostics.py','roundtrip_controls.py')) {
  & $Py "interim/2026-11-25-joint-birth-projection-a02/$script" --run runs/NEW-a02-PC-n4
  if ($LASTEXITCODE -ne 0) { throw "後段監査失敗: $script" }
}
```

5〜20分程度を目安にする（クラウド測定と異なりPCでは未測定）。本体3経路と独立3経路を
新規求解し、観測なしの本体直接solveも別実行する。A/B閉鎖固着は各32点、Bの
開離・再接触経路は28点、合計92点の対を比較する。実checkpoint12地点から
別PIDで再開する。保存JSONを読むだけの監査を新規求解と数えない。

主な出力は次のとおり。

|ファイル|内容|
|---|---|
|`measurements.json`|全場差、独立平衡/仕事恒等式、別PID復旧、実入力拒否、改変拒否。A-02全体の`passed=false`を保持|
|`source-fingerprints.json` / `source-snapshot.zip`|新規求解前に固定した実行・試験ソースとSPEC|
|`binaries.json` / `bin/`|実行版SHAと固定コピー。実行中にビルドし直さない|
|各経路の`native/accepted.jsonl`, `native/run.json`|本体の全受理状態・全Run、途中checkpoint|
|`native/constitutive.jsonl`|本体公開APIによる確定状態の全節点力再組立。Newton反復診断とは異なる|
|`native-fields.jsonl`, `independent-fields.jsonl`|両側の全受理点・全50場。独立側は入力のみの新規3D求解|
|`comparison.json`, `independent-result.json`|全成分の生差、各工程差、独立平衡残差|
|`restarts/`, `whole/`|別PID復旧・全後続状態、観測なし直接Run|
|`refusals/`, `projection/`, `tamper/`|本番CLI拒否、公開API診断、外側hash再計算後も拒否したcheckpoint|

合否は「プロセス終了0」「本体/独立完走」「比較の数値判定」「A-02の完了」を分ける。
全受理格子が入力の `i/N`、欠落状態なし、世代/開閉カテゴリ一致、Runの時間以外一致、
所定の理由による拒否を確認する。旧独立弾性相対1e-7、節点力1e-10 kN、
拘束1e-10 m、接触仕事恒等式1e-10を緩めない。
全場の未定義な新許容値は作らず、非線形履歴の差は生値のまま評価する。
`comparison.json` の包括受入falseを削除したり、総合passedをtrueへ変更したりしない。
時間切れ、未収束、CIキャンセルは成功でない。

## PC専用の長い履歴系列

同じ81節点・24 Tet10・648材料GP・32接触GPの形状を使い、時間/荷重格子だけを細かくする。
今回の密行列独立参照の350自由度予算、先端の単一・乾燥・全所有要素活動のガードを
外して大きいメッシュを渡す手順ではない。

```powershell
$env:A02_PROCESS_BUDGET_SECONDS='21600' # 各子プロセス最大6時間。時間切れは失敗のまま
try {
  foreach ($n in @(16,32,64)) {
    & $Py $Study --steps $n --output "runs/NEW-a02-PC-n$n"
    if ($LASTEXITCODE -ne 0) { throw "n=$n が未完了。出力を保持して原因を確認" }
  }
} finally { Remove-Item Env:A02_PROCESS_BUDGET_SECONDS -ErrorAction SilentlyContinue }
```

16/32/64系列はそれぞれ368/736/1472状態対。observerが毎受理点でseal/resume認証を行うため
単純なN比例ではなく、20〜90分、1〜4時間、3〜12時間程度の計画枠を取る
（PC未測定の予算目安）。各プロセスの実測時間・終了コードを出力で確認する。
同じ進行度の全場、接触の実分割、外力仕事・塑性/数値散逸を比較し、
前系列との差が残れば収束未達とする。新しい収束許容値を観測後に作らない。
より細かい先端メッシュの既存局所h検証を行う場合は、別問題として次を実行する。

```powershell
& $Py tools/verify_joint_crack_tip_refinement.py --solver target/release/soil-fem-3d.exe --output runs/NEW-a02-PC-tip-local-h
if ($LASTEXITCODE -ne 0) { throw '先端局所h検証失敗' }
```

この既存h系列はBIRTH初期力投影の受入へ転用しない。実CAD/GUI/GPUや一般の大規模
3D施工には、この小型参照で対応していない入力・出力契約が必要であり、今回の受入外。
大規模計算をクラウドで代行していない。

## 配布した生証拠の整合確認

```powershell
& $Py interim/2026-11-25-joint-birth-projection-a02/package_evidence.py --verify interim/2026-11-25-joint-birth-projection-a02/evidence-manifest.json --scratch runs/NEW-a02-archive-check
if ($LASTEXITCODE -ne 0) { throw '主証拠のSHA不一致' }
& $Py interim/2026-11-25-joint-birth-projection-a02/package_evidence.py --verify interim/2026-11-25-joint-birth-projection-a02/production-cli-evidence/evidence-manifest.json --scratch runs/NEW-a02-standard-cli-archive-check
if ($LASTEXITCODE -ne 0) { throw '標準CLI対照のSHA不一致' }
```

全part・結合ZIP・全エントリーを照合する。これは保存証拠の整合確認であり、
上の本体/独立の新規求解とは別である。

## Linuxの今回と同じ実行

```bash
source /workspace/soil-environment/activate.sh
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo test --locked --release
cargo test --locked --release --examples
cargo build --locked --release --examples
cargo build --locked --release
python interim/2026-11-25-joint-birth-projection-a02/study.py --steps 4 --output runs/NEW-a02-linux
python tools/run_ci_python_tests.py --group python-core-0 --output runs/NEW-a02-tests/core0
python tools/run_ci_python_tests.py --group python-core-1 --output runs/NEW-a02-tests/core1
python tools/run_ci_python_tests.py --group python-mean-evidence --output runs/NEW-a02-tests/mean
python tools/run_ci_python_tests.py --group python-spatial-evidence --output runs/NEW-a02-tests/spatial
```

既存checkoutを部分取得している場合、各試験が読む元アーカイブ・元ソースを元Git blobと
一致するバイトで配置する。欠落による失敗を数値失敗へ読み替えず、初回失敗も保持し、
配置後は該当する通常群全体を再実行する。旧hash/snapshotは作り直さない。
