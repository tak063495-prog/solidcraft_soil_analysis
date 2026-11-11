# Windows PCでの再現・大きな系列・将来比較

A-01・NC-01・P6は未達のまま。以下は乾燥・固定世代の独立参照の手順である。水圧、Richards、先端BIRTH投影、板オフセットは対象外。既存ソース・試験・CI・許容値・参照値を編集しない。出力先は毎回新しいディレクトリを使い、失敗した出力を上書きしない。

## 分離した作業場所と依存

PowerShell、Git、Python 3.13を使う。既存リポジトリから別worktreeを作る。下記パスは利用者の配置に合わせる。`origin` の対象ブランチはpush済みのコミットを確認してから用いる。

```powershell
git fetch origin codex/py-multi-joint-structure-composite-reference
git worktree add ..\dry-multi-joint-reference origin/codex/py-multi-joint-structure-composite-reference
Set-Location ..\dry-multi-joint-reference
$refDir = 'interim/2026-11-21-py-multi-joint-structure-composite-reference'
py -3.13 -m venv .venv-dry-ref
.\.venv-dry-ref\Scripts\Activate.ps1
python -m pip install -r "$refDir/requirements.txt"
if ($LASTEXITCODE -ne 0) { throw 'dependency installation failed' }
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
New-Item -ItemType Directory -Force work | Out-Null
python -V
python -m pip freeze | Out-File -Encoding utf8 work/dry-environment.txt
```

新規参照の実行にはRust・gmsh・既存Pythonソルバは不要。既存全Python試験は本体と旧証拠を必要とするため、後述の別手順へ分ける。

## 新規試験と小型80系列

```powershell
python -m pytest "$refDir/checks_reference.py" -q -rA --junitxml=work/dry-new-tests.xml
if ($LASTEXITCODE -ne 0) { throw 'new reference tests failed' }
python "$refDir/study.py" --output work/dry-small
if ($LASTEXITCODE -ne 0) { throw 'small study incomplete' }
python "$refDir/report.py" --studies work/dry-small --output work/dry-small-summary
if ($LASTEXITCODE -ne 0) { throw 'raw history audit failed' }
python "$refDir/audit.py" --input work/dry-small/cross3d-steps32-gp2-identity.json --output work/dry-cross3d-audit.json
if ($LASTEXITCODE -ne 0) { throw 'independent audit failed' }
```

新規試験は26件（subtestは別試験として加算しない）、skip/xfailなし。小型studyは10ケース×8 Run＝80 Run、6,480確定増分、各工程8/16/32刻みとGP次数1/2/3、接触逆順・JOINT逆順・節点置換を実求解する。クラウド実測は参照80系列合計約10分（他の試験との並行実行を含む）；単独のWindows PCでの時間は未測定で、概算2～20分、監査1～5分、テキスト生出力約1.1 GBを見込む。環境により増える。正確な時間は各Runとstudyの `elapsed_seconds` を残す。

出力は全節点の変位・bulk/contact/structure内力・外力/反力/残差、全接触GPの位置/法線/面積/traction/力/jump/すべり/開閉/世代、全bulk要素のstrain/stress、全部材の軸力/端力/ゼロ端モーメント/弾性エネルギー、全仕事台帳である。`results/study.json` と `observations.json` が集計で、個別生出力を代替しない。

## 合否・診断の扱い

- ソフトウェア検査は終了コード0、全26件実行、fail/error/skip/xfail 0を要求する。`complete=false`、例外、時間切れ、キャンセル、部分Runは成功ではない。
- 再組立の全節点力絶対1e−10 kN、接触仕事絶対1e−10 kN mは、元 `2026-11-13-joint-ext-combos/acceptance.json:tolerances` とREADMEの引用である。範囲拡張の受入へ転用しない。Newtonの2e−11 kNは求解停止制御であって、新しい数値受入値ではない。
- 全節点平衡・応力/接触力の別散布・軸力・離散仕事恒等式を個別に読む。外力仕事は必ず生の `W−U−D` を確認する。離散恒等式の追加項を外力仕事へ足して合格にしない。今回のすべり履歴は仕事1e−10へ未達である。
- すべり・開閉・モーメント・断面力・順序・刻み/GP収束の新しい許容値は**未定義**、`accepted=null`。順序の丸め差、負の観測次数、非単調、GP非一致点を消さない。系列差の減少だけで合格にしない。
- 既存の限定受入・未達・失敗記録を保持する。独立モデルの完了を先端BIRTH、一般3D、実CAD、A-01、NC-01、P6の完了へ読み替えない。

## PC限定の長い時間／GP系列

以下は**クラウド未実行**。各工程32/64/128とGP次数1/2/3/4を全ケースで実行する。小型mesh2/mは固定し、まず時間差と接触積分差を分ける。途中から良い水準だけ抜き出さない。

```powershell
python "$refDir/study.py" --steps 32 64 128 --gp-orders 1 2 3 4 --output work/dry-long-mesh2
if ($LASTEXITCODE -ne 0) { throw 'long study incomplete; retain logs and partial data' }
python "$refDir/report.py" --studies work/dry-long-mesh2 --output work/dry-long-summary
if ($LASTEXITCODE -ne 0) { throw 'long study audit failed' }
```

Windows所要時間は未測定・概算20分～3時間、出力5～15 GB。最細の128/256/512などへ進む前に、32/64/128全履歴の平衡、開閉切替、仕事差、すべり差を記録する。新しい合格値や外挿参照を作らない。

## PC限定の大きな3D規則格子

bulkのhと時間/GPを同時に変えず、まず各工程64・GP次数2でmesh2/4/6/mを保存する。独立モデルは密行列なので、mesh6を超える前に実測RAMと時間を確認する。32 GB以上のRAMと十分な空きディスクを推奨する。数字は必要性の目安で、実測保証ではない。

```powershell
foreach ($caseName in @('parallel3d','cross3d','t3d','shared3d','anchor3d')) {
  foreach ($bulkMesh in @(2,4,6)) {
    python "$refDir/reference.py" --case $caseName --steps 64 --gp 2 --mesh $bulkMesh --output "work/$caseName-mesh$bulkMesh-step64-gp2.json"
    if ($LASTEXITCODE -ne 0) { throw '3D reference incomplete' }
    python "$refDir/audit.py" --input "work/$caseName-mesh$bulkMesh-step64-gp2.json" --output "work/$caseName-mesh$bulkMesh-audit.json"
    if ($LASTEXITCODE -ne 0) { throw '3D audit failed' }
  }
}
```

mesh4で概算10分～2時間、mesh6で1～12時間／ケース、全ケースは数時間～数日を見込む。**未実測の概算**で、時間切れを数値停止や合格にしない。`reference.py` は例外時に成功JSONを出さないのでstderrも保存する。meshを変えた場合、節点/GPは一致しない。現比較器は自動補間をせず拒否する。全場のh誤差ノルム・対応点・許容値は未定義であり、この手順の出力だけではh収束合格にならない。

## 既存拒否証拠の読取再現

```powershell
python "$refDir/scope_audit.py" --root . --archive interim/2026-11-13-joint-ext-combos/evidence.zip --output work/dry-scope-audit.json
if ($LASTEXITCODE -ne 0) { throw 'historical evidence changed or missing' }
```

元archiveのSHAと48範囲拒否エントリーを既存manifestに照合する。12の保存拒否ログは湿潤BIRTH/Richards。乾燥の複数JOINT／構造複合の具体的な旧拒否入力・ログは未確認を維持する。現行ガードのソース読取と、CLIを実行して得た拒否は同じ証拠ではない。ここでは本体を呼ばない。

## 将来の本体CSV/JSONとの比較

既存本体が対象複合を拒否している現段階では、対応する本体の完成Run/exporterは**未確認**。次は規約を変更せずに生exportを実装した別作業の出力を受け取る手順であり、今すでに本体比較が成功したという意味ではない。

`README`のschemaを満たす `production-dry.json` と、接触評価順を変えて独立に本体を実行した `production-dry-reordered.json` を用意する。同じ物理・離散化・負荷端点・節点/GPの物理IDを保持する。JOINT IDだけは次の明示対応を作る。実際のIDで全単射を記載し、未対応や最近傍対応を作らない。

```powershell
@'
{"ProductionJointA":"Jx=1","ProductionJointB":"Jx=2"}
'@ | Set-Content -Encoding utf8 work/joint-id-map.json
python "$refDir/compare.py" --reference work/dry-small/beam2d-steps32-gp2-identity.json --candidate work/production-dry.json --order-run work/production-dry-reordered.json --joint-map work/joint-id-map.json --output work/production-dry-comparison.json
$comparisonExit = $LASTEXITCODE
if ($comparisonExit -eq 2) { Write-Host '拒否または許容値未定義。JSONのreasonと全生差を確認。' }
elseif ($comparisonExit -ne 0) { throw 'comparison execution failed' }
```

UTF-8 BOM付きJSONも通常の文字コードとして読める。値の変換・単位推測は行わない。

CSVの場合は4列の生payloadと `.manifest.json`（全ID/規約/評価順/フレーム台帳/CSV SHA）を用意して次を実行する。exportを生成する本体CLIの具体的コマンドは未確認なので捏造しない。

```powershell
python "$refDir/compare.py" --reference work/dry-small/beam2d-steps32-gp2-identity.json --candidate work/production-dry.csv --manifest work/production-dry.manifest.json --order-run work/production-dry-reordered.csv --order-manifest work/production-dry-reordered.manifest.json --joint-map work/joint-id-map.json --output work/production-dry-csv-comparison.json
```

単位/符号、欠損/重複/未対応ID、GP/世代、平滑化/クリップ宣言、raw欠落、力台帳、順序別出力欠落を拒否する。数値差の許容値が未定義なので、丸め差も含めて全差を保存しexit 2とする。自己往復のexit 0も `accepted=null`。JSON/CSV最終値だけから秘密の補正全てを検出できるとは主張せず、入力・出力ハッシュと信頼できる生exportを保存する。

## 実CAD

実CADは今回のクラウドで読込・求解していない。現在の独立参照は上記の明示規則格子だけを受け付け、CADを最近傍節点・平面化・GP補間へ自動変換しない。実CAD→本参照へのadapter／具体的CLI、独立参照値、対応する一般3D許容値は**未実装・未確認**。

PCで別途用意するものは、原CADとSHA、単位・座標系、各弾性体のTri3/Tet4接続と境界条件、JOINT両側/物理ID/法線/GP/積分重み、共有節点か複製節点かの宣言、梁/アンカーのEAと取付点、全負荷端点、rawexportと順序別exportである。これらを保持するadapterと独立メッシュ検証を先に用意する。元本体が拒否する組合せは拒否ログを保存し、受理できたことにしない。データが揃っても本比較器は物理・離散化・ID不一致を拒否する。この未確認手順を実CADの合格や実行済みのコマンドとして報告しない。

## Linux／Windowsでの既存全Python

既存CIと同じ `tests/`＋`tools/` のunittest全発見を4群で実行する。全発見543件、割当272/263/4/4（この起点）。新しい26件はこの件数に含めない。アーカイブを省略したsparse checkoutで走らせない。既存lockの依存とRust 1.99.1を用意する。

```powershell
python -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw 'existing dependencies failed' }
cargo build --locked --release
if ($LASTEXITCODE -ne 0) { throw 'existing release build failed' }
foreach ($testGroup in @('python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence')) {
  python tools/run_ci_python_tests.py --group $testGroup --output "work/existing-$testGroup"
  if ($LASTEXITCODE -ne 0) { throw "existing group failed: $testGroup" }
}
```

Linuxでは同じpythonコマンドに `OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1` を指定する。この環境のようにRustを別配置した場合は `CARGO_HOME`、`RUSTUP_HOME`、PATHも既存インストールへ向ける。通常のWindowsでは利用者のインストールを用いる。

既存群は概算15～60分／通常群、封印証拠群は5～30分／群、旧アーカイブと展開で10～30 GB以上が必要。新しい大規模NCを求解する手順ではなく、既存の小型試験とSHAで封印した旧証拠の再監査である。各群の `python-test-inventory.json` で全ID、実行件数、エラー・失敗・skip等を確認する。setUpClassエラーで未実行となった試験を成功件数へ加算しない。Linuxのバイト一致検査で差が出た場合、`allow_platform_float_differences` を使って通したり旧JSONを書き換えたりしない。環境差と生差を保存し、不合格として報告する。

draft PR作成後のWindows CIは、当該先端のhead/merge SHAと各jobの結果・件数を別に確認する。旧CIや時間切れ・キャンセル・途中の成功を今回の成功に転用しない。
