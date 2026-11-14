# P6-03 再現・ユーザーPC手順

保存先2026-12-07は依頼ラベル、実施日は2026-11-12 JST。PR #13の
`e3d830d8d1bdffd727c4ab7271755f52e556006d` を起点とする。PRはマージしない。
以下はリポジトリルートのPowerShell。出力は必ず未存在の場所を使う。
所要時間はPC未実測の予算目安で、上限・性能保証ではない。時間切れ、未収束、
キャンセル、infra失敗、部分出力は未完了としてログを残す。

## 共通準備と短い監査

```powershell
git switch codex/plastic-spatial-convergence
$P='interim/2026-12-07-p6-03-cloud-reconciliation'
$Run='runs/NEW-PC-p603'
New-Item -ItemType Directory -ErrorAction Stop $Run | Out-Null
py -3.13 -m venv .venv-p603
$Py='.venv-p603/Scripts/python.exe'
& $Py -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw 'dependency setup failed' }
$env:CARGO_BUILD_JOBS='2'; $env:RUST_TEST_THREADS='2'
$env:OPENBLAS_NUM_THREADS='1'; $env:OMP_NUM_THREADS='1'; $env:MKL_NUM_THREADS='1'
$env:PYTHONUTF8='1'; $env:PYTHONIOENCODING='utf-8'
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
cargo build --locked --release --bins --examples
if ($LASTEXITCODE -ne 0) { throw 'build failed' }
$Solver=(Resolve-Path target/release/soil-fem-3d.exe).Path
& $Py "$P/audit_status.py" > "$Run/audit.json"
if ($LASTEXITCODE -ne 0) { throw 'evidence/state audit failed' }
& $Py -m unittest discover -s $P -p test_audit.py -v 2>&1 | Tee-Object "$Run/audit-tests.log"
if ($LASTEXITCODE -ne 0) { throw 'audit tests failed' }
& $Py "$P/compare_mc_crest.py" "$Run/mc-crest" > "$Run/mc-crest-driver.json"
if ($LASTEXITCODE -ne 0) { throw 'MC crest comparison incomplete/failed' }
```

ビルド2〜15分、監査・天端比較は数秒〜数分。`audit.json`の`integrity_ok`と
`mc-crest/comparison.json`の3工程・元変位許容値`5e-8 m`（厳密な`<`）を見る。
これはMC柱の天端だけで、全場の再受入、MCC、Tet4の6.828%不合格の解消ではない。
1.034%は別の条件付き連続推定で、離散解を合格にしない。
本体が出す元5%・最大差10%・奥行き差1%の判定は別に保持する。

取り込んだ参照は `reference-sources.json` のPR/コミット/元パス/SHA/blobで固定した
未変更コピー。参照README内の相対リンク・旧実行時点は元のままである。
本手順のコマンドを入口にし、他PRへpushやマージをしない。

## MCC柱と既存MCC掘削空間差（PC）

**元.2m柱は非関連MC。別の1m MCC柱と混ぜない。** MCCはK0=1、
p0=pc0=100+20(1−z) kPa、γ=20、λ=.2、κ=.04、M=1.2、ν=.3、N=3.5、
水平変位固定、底uz=0、乾燥。130→100 kPaの載除荷、上半分撤去後は新上面100 kPa。
`prepare_pc_mcc.py` は既存準備器の入力定義を使うだけで求解しない。

```powershell
& $Py "$P/prepare_pc_mcc.py" "$Run/mcc-inputs"
if ($LASTEXITCODE -ne 0) { throw 'MCC input preparation failed' }
foreach ($c in 'space8','space16','space32','load16','load32','load64','mat1e4','mat1e5','mat1e6') {
  & $Solver solve-stages "$Run/mcc-inputs/$c/mesh.json" "$Run/mcc-inputs/ground.json" `
    "$Run/mcc-inputs/$c/stages.json" "$Run/mcc-inputs/$c/result.json" 2>&1 |
    Tee-Object "$Run/mcc-inputs/$c/solve.log"
  if ($LASTEXITCODE -ne 0) { throw "MCC $c incomplete; retain evidence" }
}
& $Py "$P/references/pr24/run_reference.py" --repo . --output "$Run/mcc-independent" `
  --column-spatial-cells 16 32 64 --column-material-levels 512 1024 2048 --column-spatial-material-n 4096
if ($LASTEXITCODE -ne 0) { throw 'independent MCC incomplete' }
# 元の非関連MC柱の別系列。MCCへ転用しない。
& $Solver plane-strain-2d user "$Run/original-mc-user.json"
if ($LASTEXITCODE -ne 0) { throw 'MC column incomplete' }
# 既存大規模掘削は旧JSONのhash付き読取。新求解ではない。
& $Py tools/summarize_mcc_spatial_records.py "$Run/mcc-legacy-inventory.json"
if ($LASTEXITCODE -ne 0) { throw 'legacy MCC bytes differ' }
```

MCC本体9経路30分〜数時間、独立1Dは数分〜30分、RAM数GB、非関連MCは数分〜30分を見込む。
共通物理点で変位、6応力、p/q/pc/e、体積ひずみ、受理刻み、撤去履歴凍結、全工程を確認する。
自動cutbackで意図した刻み列でなければ未確認。空間・荷重刻み・材料刻みを同時に変えて
空間次数を作らない。旧MCC柱のpc差`1.773e-3 kPa`不合格と元目標`1e-4 kPa`は保持する。
新しい独立1D材料/空間診断は許容値未定義で`accepted=null`。
本体の全場を共通標本へ根拠付き出力できた後だけ、PR24の`assess_series.py`/
`compare_history.py`を使う。共通点写像・不足GP量を推測補完しない。

既存の274〜14909要素の掘削は、この1D柱で置換できない。最後の固定壁の
1.463403911→1.527772593%非減少を保持する。数GB生結果・元の凍結入力/実行版を
PCに揃えた後、[MCC_MESH_OVERLAP](../../docs/MCC_MESH_OVERLAP.ja.md) と
[DEPTH_SPATIAL_DIAGNOSTICS](../../docs/DEPTH_SPATIAL_DIAGNOSTICS.ja.md)の元入口で全工程を監査する。
未確認の生結果パスや新規メッシュ生成規則をここで発明しない。この実斜面MCCの追加求解コマンドは
元の入力一式を取得・確認するまでPC手順待ちであり、クラウド未実行。

## 押出し応力・静的振動モニター・有限幅奥行き応力（PC）

旧MC立方体の材料/荷重/拘束/撤去位置は既存ドライバのまま。
物理監査は力残差≤1e−8、全体力誤差≤1e−8、線形許容1e−10、
正規化降伏≤1e−7、工学ひずみ≤.05。

```powershell
& $Py tools/plastic_spatial_convergence.py "$Run/plane-4-8-16" --divisions 4 8 16 --timeout 3600
if ($LASTEXITCODE -ne 0) { throw 'plane series incomplete' }
& $Py tools/plastic_spatial_convergence.py "$Run/plane-8-16-32" --divisions 8 16 32 --timeout 3600
if ($LASTEXITCODE -ne 0) { throw 'finer plane series incomplete' }
& $Py tools/plastic_spatial_convergence.py "$Run/plane-half" --divisions 4 8 16 `
  --increment 0.03125 --material-substep 0.00005 --control-baseline "$Run/plane-4-8-16" --timeout 3600
if ($LASTEXITCODE -ne 0) { throw 'half-step control incomplete' }
foreach ($xz in 2,4) {
  & $Py tools/plastic_spatial_convergence.py "$Run/finite-depth-xz$xz" --axis depth --finite-width `
    --in-plane-divisions $xz --divisions 4 8 16 --timeout 3600
  if ($LASTEXITCODE -ne 0) { throw 'depth series incomplete' }
}
& $Py tools/plastic_spatial_convergence.py "$Run/finite-inplane" --finite-width --divisions 4 8 16 --timeout 3600
if ($LASTEXITCODE -ne 0) { throw 'finite in-plane series incomplete' }
```

各列30分〜数時間以上、交差体積積分に数GB以上を見込む。`--timeout`は各求解だけの上限で
交差積分の上限ではない。32分割は押出し6144/有限幅12288 Tet4。
`report.json`の全テンソルRMS/相対L2、共通位置(.5,0,.5)mの全変位・初期値・掘削増分、
物理ゲート、`richardson`の拒否理由を見る。元の必要条件は差比`<.7`、次数[.5,8]、
スカラーは非零・同符号差、値の1e−8帯を超える変化。三水準だけの条件付き推定を
漸近収束や工学的応力精度の合格にしない。絶対応力許容値は未定義。
有限幅奥行き列は物理長1mを固定したメッシュ細分化で、物理長変更とは別。
静的モニターに物理秒・質量・減衰を後付けしない。

PR33の弾性箱/Lamé/動的検算法を詳しく見る場合だけ、別問題として次をPCで実行する。
数分〜60分以上、数GB〜16GB以上。元MCの合否には使わない。

```powershell
& $Py "$P/references/pr33/run_reference.py" "$Run/elastic-depth-reference" `
  --widths 1 --lengths 4 8 16 32 --levels 1 2 --ring-levels 4 8 16
if ($LASTEXITCODE -ne 0) { throw 'separate reference incomplete' }
```

## 杭（PC）

元10×10×8m地盤・1m角6m杭・完全付着・乾燥・側面ローラー・底uz=0、
自重後(+10,0,−100)kN増分を保持。旧Tet10は436/819/2279、
杭頭面積平均の追加沈下差6.687→4.717%。厳密解誤差ではない。

```powershell
& $Py "$P/references/pr28/solid_pile_reference.py" --pc-series --subdivisions 2 4 8 `
  --output "$Run/pile-hex8-through-m8.json" 2>&1 | Tee-Object "$Run/pile-m8.log"
if ($LASTEXITCODE -ne 0) { throw 'pile reference incomplete' }
# m8が完了しRAM/時間予算を確認してから別系列を実行
& $Py "$P/references/pr28/solid_pile_reference.py" --pc-series --subdivisions 4 8 16 `
  --output "$Run/pile-hex8-through-m16.json" 2>&1 | Tee-Object "$Run/pile-m16.log"
if ($LASTEXITCODE -ne 0) { throw 'pile m16 incomplete' }
```

m8:数分〜1時間/RAM数GB〜16GB、m16:20分〜4時間/RAM16〜64GB以上の未実測目安。
`complete=true`、体積794/6m³、杭頭面積1m²、荷重/反力、残差、追加沈下とhを見る。
新Hex8は近似参照で空間収束許容値未定義、`accepted=null`。
旧同一Tet10メッシュの基準は全節点変位相対差<1e−7、応力差<1e−7、e差<1e−8、
真残差<1e−10、合力atol=1e−5kN/rtol=1e−9（`verify_solid_pile.py`）。
これを異なるHex8空間差の基準へ転用しない。1D棒/梁のばね係数を3D地盤Eや旧沈下にfitしない。
旧本体のより細かいTet10/CAD・先端周面分担・塑性/すべりは元入力/凍結実行版と
[杭手順](../../docs/SOLID_PILE.ja.md)を揃えた後のPC残件。新しい支持力モデルは作らない。

本体Tet10も同じ有限地盤のまま、既存PowerShell入口を使用できる。
SolidCraft 0.20.0 CLIをPCへ既存配布手順で配置してから実行する。各30分〜数時間以上、
CAD/メッシュ生成/全場/再開検証の時間を含む。要素数はGmshで変わるので旧436等との一致を仮定しない。

```powershell
$CadCli='tools/solidcraft-0.20.0/solidcraft-0.20.0-windows-x64-portable/solidcraft-cli.exe'
pwsh -NoProfile -File scripts/Invoke-SolidPile.ps1 -SolverPath $Solver -CliPath $CadCli `
  -PythonPath $Py -OutputDirectory "$Run/pile-tet10-original-pair" -CoarseMeshSizeM 2 -FineMeshSizeM 1.4
if ($LASTEXITCODE -ne 0) { throw 'original-size pile pair incomplete' }
pwsh -NoProfile -File scripts/Invoke-SolidPile.ps1 -SolverPath $Solver -CliPath $CadCli `
  -PythonPath $Py -OutputDirectory "$Run/pile-tet10-finer-pair" -CoarseMeshSizeM 1.4 -FineMeshSizeM 1.0
if ($LASTEXITCODE -ne 0) { throw 'finer pile pair incomplete' }
```

`independent-verification.json`の全節点/GP・面積平均追加沈下、`history/verification.json`を見る。
上記同一メッシュ基準を保持する。異なる空間水準の工学的合格基準は未定義のまま。

## 降雨（PC）

中央0.25m²正方形雨・VG/Mualem・初期全水頭H=0（ψ=−z）・底H=0・側面無流束、
雨3e−6m/sを0〜3600秒、休止〜7200秒、再降雨〜10800秒。物性はPR32の元5行定義のまま。
1D/帯状2D/等面積円に置換しない。

```powershell
foreach ($n in 8,16,32) {
  & $Py "$P/references/pr32/pc_rainfall.py" --user-pc --dim 3 --n $n --dt 100 `
    --out "$Run/rain-3d-n$n.json"
  if ($LASTEXITCODE -ne 0) { throw "rain n=$n incomplete" }
}
foreach ($dt in 100,50,25) {
  & $Py "$P/references/pr32/pc_rainfall.py" --user-pc --dim 3 --n 16 --dt $dt `
    --out "$Run/rain-3d-n16-dt$dt.json"
  if ($LASTEXITCODE -ne 0) { throw "rain dt=$dt incomplete" }
}
```

本体の旧RT0系列は次で準備・実行する。元`rain-original/run.py`の入力定義だけを新規出力へ
転記する準備器で、旧記録や旧ソルバSHAを上書きしない。今回はPC専用・未実行。

```powershell
& $Py "$P/prepare_pc_rainfall.py" "$Run/rain-body-inputs"
if ($LASTEXITCODE -ne 0) { throw 'rain RT0 preparation failed' }
foreach ($c in 'rt0-4-dt100','rt0-8-dt100','rt0-16-dt100','rt0-8-dt50','rt0-8-dt25') {
  & $Solver solve-richards "$Run/rain-body-inputs/$c/mesh.json" "$Run/rain-body-inputs/$c/case.json" `
    "$Run/rain-body-inputs/$c/result.json" 2>&1 | Tee-Object "$Run/rain-body-inputs/$c/solve.log"
  if ($LASTEXITCODE -ne 0) { throw "RT0 $c incomplete" }
  & $Py tools/verify_richards.py "$Run/rain-body-inputs/$c/mesh.json" `
    "$Run/rain-body-inputs/$c/result.json" "$Run/rain-body-inputs/$c/audit.json"
  if ($LASTEXITCODE -ne 0) { throw "RT0 $c discrete audit failed" }
}
```

本体RT0は各数分〜数時間（24576 Tet4を含む）、RAM数GB以上。
独立FV8³:秒〜数分、16³:数分〜数十分、32³:数十分〜数時間/RAM数GB。
出力の`complete`、10800秒到達、貯留/累積雨量/底排水、最大収支欠損を確認。
独立FV累積収支の元基準は`<1e−9 m³`。旧RT0の6tetブロック平均、FVの上下配列、
元の三線形標本化とstorage offsetを監査した後だけ、本体4/8/16空間列、
本体8固定100/50/25秒の生出力と照合する。本体の元出力一式がない間は比較不能。
PC上で旧 `rain-original/compare.py` の条件を保持し、次の5値を再集計する。

|元量|旧値|厳密な元上限|
|---|---:|---:|
|3600秒・最細水頭RMS m|0.0046205475457005325|<.003|
|3600秒・本体空間差比|1.348244318335049|<.8|
|7200秒・本体空間差比|.998450823613834|<.8|
|10800秒・最細水頭RMS m|.0038863726982864173|<.003|
|10800秒・総貯留絶対差m³|.00014219441327791016|<.0001|

本体8/16/24・FV64の後続29水理ゲートは別の窓。湛水/流出、一般地質、塑性/湿潤接触は残る。
新独立3Dの完了だけで旧5値の合格を宣言しない。

## Linux / Windows全件確認とfixture不足

cargo全件20〜60分以上、Python10〜60分以上、旧NC archive約6.4GBの取得・展開時間は別。
全コマンドの終了値を記録し、失敗しても他群は実行する。

```powershell
cargo fmt --all --check 2>&1 | Tee-Object "$Run/fmt.log"
$fmt=$LASTEXITCODE
cargo clippy --locked --all-targets -- -D warnings 2>&1 | Tee-Object "$Run/clippy.log"
$clippy=$LASTEXITCODE
cargo test --locked --release --all-targets --no-fail-fast 2>&1 | Tee-Object "$Run/rust.log"
$rust=$LASTEXITCODE
cargo test --locked --release --doc 2>&1 | Tee-Object "$Run/rust-doc.log"
$doc=$LASTEXITCODE
foreach ($g in 'python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence') {
  & $Py tools/run_ci_python_tests.py --group $g --output "$Run/$g" 2>&1 | Tee-Object "$Run/$g-driver.log"
  Write-Host "$g exit=$LASTEXITCODE"
}
@{fmt=$fmt;clippy=$clippy;rust=$rust;doc=$doc} | ConvertTo-Json | Out-File "$Run/rust-exits.json"
```

疎配置で欠けたfixtureは初回ログを保存し、`git show HEAD:元パス`の元blobをバイナリで
復元する（PowerShellのリダイレクトでバイト変換しない）。復元後はSHA/blobを記録して該当群を再試験。
本記録の `fixture-recovery.json` を参照。6.4GBのNC archiveを取得できない場合は未実行と明記する。
NC mean/spatialの厳密setUpClass不一致はskip/許容差化/参照書換えで通さない。
各4メソッド未実行と準備エラー各1を別計上する。過去617/625・8未実行・準備エラー2は保持する。
新先端のWindows CIは別途headと全群の完了/成功を確認するまで未確認。
