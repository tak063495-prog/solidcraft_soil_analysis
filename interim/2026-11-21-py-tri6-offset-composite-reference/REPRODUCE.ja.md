# Windows PCでの再現

対象ブランチ `codex/py-tri6-offset-composite-reference`、起点
`f377e8eda41777c84c3bb849af556bf61683b4cb`。PRはdraftのまま、マージしない。
PowerShellをリポジトリー直下で使う。実行先は毎回新しい名前にし、保存済み結果・過去の失敗を上書きしない。

## 独立参照だけを実行

```powershell
git fetch origin codex/py-tri6-offset-composite-reference
git switch codex/py-tri6-offset-composite-reference
py -3.13 -m venv .venv-tri6-reference
$refPython = Join-Path $PWD '.venv-tri6-reference\Scripts\python.exe'
& $refPython -m pip install numpy==2.5.4 scipy==1.18.1
if ($LASTEXITCODE -ne 0) { throw 'dependency installation failed' }
$refDir = 'interim/2026-11-21-py-tri6-offset-composite-reference'
$refRun = 'runs/py-tri6-offset-pc-01'
& $refPython "$refDir/checks.py" --output "$refRun/checks"
if ($LASTEXITCODE -ne 0) { throw 'independent software checks failed' }
& $refPython "$refDir/run_reference.py" --output "$refRun/reference" --h-series 2 4 8
if ($LASTEXITCODE -ne 0) { throw 'reference incomplete' }
```

通常の8コア/16GB PCでの未実測目安：依存取得1～5分、checks数秒～1分、標準参照10秒～2分。
Linuxでは最終の標準参照を約4.5秒で実測したが、Windows時間の保証ではない。
独立側はRustのbuildも既存Pythonソルバも不要。

出力：`measurements.json`、`reference-output.json`、`reference-csv/`、
Winklerとアンカーのnごとの全節点・全要素重心JSON。
`checks/validation.json`は実行件数とfail/error/skip/xfail/xpassを保存する。
ソフトウェア検査の合格は全件実行・終了0・各失敗/skip等0。
参照生成の終了0は「完了した診断」であり、B-02受入ではない。`accepted=null`を保持する。
力釣合い・仕事・零モードの生値とh次数をREADMEの定義で確認する。
許容値・最低次数の適用定義がないため、PCの数値差から新しい合格値を作らない。

## 大きめのh系列はPCだけ

```powershell
& $refPython "$refDir/run_reference.py" --output 'runs/py-tri6-offset-pc-h-01' --h-series 8 16 32
if ($LASTEXITCODE -ne 0) { throw 'larger reference incomplete' }
```

未実測目安1～20分、RAM8～16GB以上を推奨（Python COO組立と疎LUのメモリーに依存）。
さらに`--h-series 16 32 64`は未実測10～90分、RAM32～64GB程度を見込み、所有PCで余裕がある場合だけ行う。
今回クラウドではn>8を求解していない。PCでもメモリー不足・時間切れ・キャンセルは未完了とし、部分結果から次数・成功を主張しない。
新しい系列でも同じ厚さ・材料・荷重・境界・GL5を保持する。厚さやばね定数を一緒に変えてh収束と呼ばない。

## 将来の本体出力と比較

現在の本体はオフセット板と水/アンカー等の複合を拒否する。入力拒否を外すパッチや架空の本体結果を本PRに含めない。
将来実装された本体から、READMEの契約に沿って同じ入力・n・局所基準の生出力をエクスポートする。
`actual_formulation`を実際の計算から記録し、旧方式を新方式と偽装しない。
例えば同じn=8のWinklerケースなら：

```powershell
& $refPython "$refDir/compare.py" `
  --reference "$refRun/reference/winkler-n8.json" `
  --candidate 'runs/future-native-winkler-n8/export.json' `
  --output 'runs/future-native-winkler-n8/comparison-01.json'
$compareStatus = $LASTEXITCODE
if ($compareStatus -eq 1) { throw 'output contract rejected; inspect reason' }
if ($compareStatus -ne 2) { throw 'unexpected comparison status' }
```

CSVの場合は`--candidate`に`metadata.json,nodes.csv,sections.csv`を持つディレクトリーを指定する。
候補の`source_kind`は`production-export`、単位はm/rad/kN/kN*m/kN/m/kN/kPa。
節点反力は支持→板、板応力は引張正、Mpは中立面、Msは付着面。
別のグローバル軸の本体出力は、出典付きのエクスポーターで同じ局所基準へ明示変換してから渡す。比較器は勝手に回転・反転しない。
本体が6評価点を出す場合、今回の重心r=s=1/3へ未監査の平均補間をせず、同じ点の生出力を追加するか、独立側で同じ評価点の参照を定義してから比較する。

終了1は符号・単位・座標・方式・点集合・完了状態等の拒否。
終了2は正常診断、`accepted=null`で**合否未定義**。
比較値が0でも既存限定受入の閾値を複合全体へ転用しない。
人工往復の動作確認には付属`reference-output.json`を両側へ指定できるが、本体比較と報告しない。

## 既存Python全件をPCで再実行

この部分だけは既存統合試験がRust CLIと過去の証拠を使う。独立参照器から呼ぶ経路ではない。
フルcheckoutで過去アーカイブを配置する。疎checkoutで欠けたまま試験をskipしない。
既存の開発環境の指示に従ってRust 1.99.1を準備する。

```powershell
. ./scripts/Enter-Dev.ps1
python -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw 'existing dependencies failed' }
cargo build --locked --release
if ($LASTEXITCODE -ne 0) { throw 'CLI build failed' }
cargo build --locked --release --examples
if ($LASTEXITCODE -ne 0) { throw 'existing examples build failed' }
$validationRun = 'runs/py-tri6-offset-existing-pc-01'
foreach ($group in @('python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence')) {
  python tools/run_ci_python_tests.py --group $group --output "$validationRun/$group"
  if ($LASTEXITCODE -ne 0) { throw "existing Python group failed: $group" }
}
```

未実測Windows目安：初回build5～20分、通常2群5～30分、歴史的証拠2群5～30分。
保存証拠の展開領域を含め空き20GB以上を確保する。新規の大規模NC求解は追加しない。
4群の発見ID集合が同じで、割当の重複・欠落0、実行件数の和が全発見数に一致することを確認する。
Linuxでの正確な件数は`results/validation/validation.json`を参照。
fail/error/skip/expected-failure/unexpected-successはすべて0を要求する。
CIの時間切れ・キャンセル・未完了は成功ではない。PCの成功も新先端Windows CIの成功を代用しない。

塑性土/MCC、JOINT全履歴、BIRTH/DEATH/MATSET複合の独立全場・状態移送・実checkpoint再開は別作業として残る。
NC-01・P6等を完了へ変更しない。
