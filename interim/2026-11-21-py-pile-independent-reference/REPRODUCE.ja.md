# Windows PC での再現・将来出力比較

リポジトリルートでPowerShell 7を使用する。既存のruns/interim成果物は上書きしない。依頼の起点は `7ae0a2a669e0e1df93f3dd9791bbeda056334c3e`。このPRのコミットをcheckoutして実行する。参照生成はNumPy/SciPyだけで完結し、Rust/既存Pythonソルバのビルドや呼出しは不要。既存Python全件試験は本体CLIを必要とする場合があるので後述の別工程で実行する。

元証拠をSHAで束縛するため、Gitの改行変換も差分として拒否する。Windowsで既存checkoutがCRLFへ変換されている場合は、元証拠を編集/再封印して合わせず、次のように新しい作業ツリーを作る。新規保存先が既に存在する場合は別名にする。

```powershell
git fetch origin codex/py-pile-independent-reference
if ($LASTEXITCODE -ne 0) { throw 'fetch failed' }
$PileCommit=(git rev-parse FETCH_HEAD)
git -c core.autocrlf=false worktree add --detach ../py-pile-independent-check $PileCommit
if ($LASTEXITCODE -ne 0) { throw 'worktree failed' }
Set-Location ../py-pile-independent-check
```

## 環境と小規模参照（数秒〜1分、初回導入を除く）

```powershell
$ErrorActionPreference = 'Stop'
$PSNativeCommandUseErrorActionPreference=$false
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
$PileDir='interim/2026-11-21-py-pile-independent-reference'
$PileRun='runs/py-pile-independent-' + (Get-Date -Format 'yyyyMMdd-HHmmss')
if (Test-Path $PileRun) { throw '既存出力を上書きしない' }
New-Item -ItemType Directory $PileRun | Out-Null
py -3.13 -m venv "$PileRun/venv"
if ($LASTEXITCODE -ne 0) { throw 'venv failed' }
$PilePython=(Resolve-Path "$PileRun/venv/Scripts/python.exe").Path
& $PilePython -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw 'dependency install failed' }
& $PilePython -m pip freeze | Out-File "$PileRun/packages.log" -Encoding utf8
& $PilePython --version 2>&1 | Out-File "$PileRun/python-version.log" -Encoding utf8
& $PilePython "$PileDir/pile_reference.py" --output "$PileRun/reference-results.json"
if ($LASTEXITCODE -ne 0) { throw '1D reference failed' }
& $PilePython "$PileDir/solid_pile_reference.py" --output "$PileRun/small-3d-results.json" 2>&1 |
    Tee-Object "$PileRun/small-3d.log"
if ($LASTEXITCODE -ne 0) { throw '3D reference failed' }
& $PilePython "$PileDir/diagnose_records.py" --small-3d "$PileRun/small-3d-results.json" --output "$PileRun/historical-diagnosis.json"
if ($LASTEXITCODE -ne 0) { throw 'historical diagnosis failed' }
& $PilePython -m unittest discover -s $PileDir -p checks_pile.py -v 2>&1 |
    Tee-Object "$PileRun/new-tests.log"
if ($LASTEXITCODE -ne 0) { throw 'new tests failed' }
```

初回依存導入は回線等により約2〜10分（未実測の目安）。参照に必要なのは既存ロックのNumPy/SciPy、上記は既存Python試験のため他のロック依存も導入する。Linux実測は新規44試験0.217秒、小規模3D系列約1.44秒。PC速度の保証ではない。新規の結果は `reference-results.json`、`small-3d-results.json`、`historical-diagnosis.json` とログ。既存出力名への再実行は拒否するので新しい `$PileRun` を用意する。

**判定**：新規44件が全件終了・失敗0・skip0。式の境界条件と釣合い、h次数、CSV/JSON拒否検査を確認する。1DはREADME記載の次数に近いこと、小規模3Dは体積/面積/荷重/反力と真残差を検討する。これらは検算法の確認。元の工学的しきい値は未定義で全出力 `accepted=null`、旧 `spatial_convergence_accepted=false`、P6未達を保持する。結果の単なる存在、exit code 0、小さい代数残差だけで空間収束を受入しない。

## 大きな独立3D系列（ユーザーPCだけ、クラウド未実行）

幾何/物理は変えず、基点区間の等分数を倍増する。m=8で9216 Hex8、m=16で73728 Hex8。三水準以上を指定し、同じ系列内で観測次数を出す。小規模の最大1152要素を越える指定には `--pc-series` が必須。

```powershell
# まずm=8まで。目安1〜10分、RAM 2〜8 GB以上（Windows未実測）。
& $PilePython "$PileDir/solid_pile_reference.py" --pc-series --subdivisions 1 2 4 8 `
    --output "$PileRun/pc-3d-through-m8.json" 2>&1 | Tee-Object "$PileRun/pc-3d-through-m8.log"
if ($LASTEXITCODE -ne 0) { throw 'm8系列未完了。成功扱いにしない' }

# m8実測でRAM/時間の余裕を確認してから別出力で実施。
# m16目安20分〜4時間、RAM 16〜64 GB以上（疎LU fill-in次第、未実測）。
& $PilePython "$PileDir/solid_pile_reference.py" --pc-series --subdivisions 4 8 16 `
    --output "$PileRun/pc-3d-through-m16.json" 2>&1 | Tee-Object "$PileRun/pc-3d-through-m16.log"
if ($LASTEXITCODE -ne 0) { throw 'm16系列未完了。成功扱いにしない' }
```

`solve()` ごとのJSON行をログへflushする。全系列が終了した時だけ最終JSONを書き、`complete=true` を付ける。メモリー不足/時間切れ/キャンセル時の部分ログは保存し、未達・未完了とする。m16を越える系列、10⁵級の本体3D、別CAD、境界拡大/材料/すべりはPCで別途資源予算と入力出典を決める残件。Python内の上限はm=16であり勝手に増やさない。

各水準の沈下、最大辺と杭断面h、反力、真残差、elapsedを保存する。真残差は元の同一メッシュ検証の `<1e-10` を参考情報として併記できるが、Hex8空間収束の合格しきい値にはしない。元の6.687→4.717%を消したり、別量/外挿/高次要素の結果で置換しない。新しい工学的合否しきい値が正式に定義されるまでは `accepted=null`。同物理条件の高精度/誤差帯の裏付け、境界感度、断面力と先端/周面分担は未達のまま。

## CSV/JSON比較の動作確認（1秒程度）

```powershell
foreach ($PileKind in @('axial','beam')) {
    & $PilePython "$PileDir/export_reference.py" --verification-only --kind $PileKind --method exact `
        --output "$PileRun/$PileKind-exact.json"
    if ($LASTEXITCODE -ne 0) { throw 'exact export failed' }
    & $PilePython "$PileDir/export_reference.py" --verification-only --kind $PileKind --method fe --n 16 `
        --output "$PileRun/$PileKind-fe.json"
    if ($LASTEXITCODE -ne 0) { throw 'FE export failed' }
    & $PilePython "$PileDir/compare_pile.py" --reference "$PileRun/$PileKind-exact.json" `
        --actual "$PileRun/$PileKind-fe.json" --output "$PileRun/$PileKind-comparison.json"
    if ($LASTEXITCODE -ne 0) { throw 'comparison refused or failed' }
}
```

これは `verification-only` caseの比較で、旧3D杭の同定ばねを用いた結果ではない。`--metadata verified-case-metadata.json` を使う1D参照出力は、そのばね/先端剛性の出典とモデル対応が別途確認できてから。exporterは乾燥・線形載荷増分・明示した頭せん断/軸力・ゼロ頭モーメント・規定先端条件・固定基底だけを実装する。その他は拒否する。未知値を確認用パラメータで埋めて元3Dのcaseを名乗らせない。

## 将来の本体出力を比較する（契約が整った後）

本体の実CSV/JSONのschema、節点/断面/界面ID、座標・符号・単位、荷重フレーム、増分/総量、先端/周面の定義を出力仕様とともに確認する。未確認の間は未達。変換後も全量の完全なexpected_idsを保ち、欠損を隠す削除や交集合化をしない。参照と本体に共通の検証済みID/座標/断面標本がなければ比較不能。

```powershell
# 実ファイルから検証済みの明示変換を終えた後に使用。
# CSVのヘッダー: id,quantity,x_m,y_m,z_m,value
& $PilePython "$PileDir/compare_pile.py" --reference "$PileRun/verified-reference.json" `
    --actual "$PileRun/pile-fields.csv" --actual-metadata "$PileRun/pile-fields.metadata.json" `
    --output "$PileRun/actual-comparison.json"
if ($LASTEXITCODE -ne 0) { throw '比較拒否。単位/符号/基底/ID/物理/工程を確認' }

# JSON版は{metadata,rows}を一体化し、sidecarは渡さない。
& $PilePython "$PileDir/compare_pile.py" --reference "$PileRun/verified-reference.json" `
    --actual "$PileRun/pile-fields.json" --output "$PileRun/actual-json-comparison.json"
if ($LASTEXITCODE -ne 0) { throw '比較拒否' }
```

比較器は単位換算や座標回転を自動で行わない。元出力・変換入力・出力SHA、座標写像、節点変位か面積平均か、断面力/応力/周面摩擦/反力の意味をmetadataとprovenanceへ明記。正規化量・符号・座標・ID・条件が違う場合のexit codeは2で、合格JSONは生成しない。規約一致なら差分JSONを生成するが、しきい値は未定義なので `accepted=null`。新規44試験にCSV/JSONの正常比較と不一致拒否の実行例がある。

## 既存Python全件（本体の通常検証、独立参照とは別）

現CIと同じ `tests/`・`tools/` の収集対象を、以下の2コマンドで全件実行する（CI自体は4群に分割）。既存試験が使うCLIを対象コミットから通常手順で用意する。この工程は本体を検証するもので、参照生成の入力には使わない。時間の目安はCLIビルド数分〜15分、既存Python数分〜数十分（PC/既存キャッシュで変動、Windows未実測）。既存全件は約8 GBの保存証拠を読む試験を含むため、取得時間を別に見込み、展開/再監査用に数十GBのディスク余裕を確保する。

```powershell
cargo build --locked --release
if ($LASTEXITCODE -ne 0) { throw 'CLI build failed' }
& $PilePython -m unittest discover -s tests -p 'test_*.py' -v 2>&1 |
    Tee-Object "$PileRun/existing-python-tests.log"
$PileTestsExit=$LASTEXITCODE
& $PilePython -m unittest discover -s tools -p 'test_*.py' -v 2>&1 |
    Tee-Object "$PileRun/existing-tools-tests.log"
$PileToolsExit=$LASTEXITCODE
if ($PileTestsExit -ne 0 -or $PileToolsExit -ne 0) { throw 'existing Python failed。両ログを保持し未達とする' }
```

ログ末尾の実行件数・失敗/skipを確認し、本PRのLinux全件件数（`results/validation.json`）と対照する。既存テスト/CIの件数条件・収集設定を変更しない。Windows固有差で件数が変わる場合は理由を明記し、skipで通さない。GitHub Windows CIがタイムアウト/キャンセルなら成功とは書かない。A-01、NC-01、P6などの完成受入とは別に報告する。
