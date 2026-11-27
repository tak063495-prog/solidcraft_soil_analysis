# Windows PCでの再現・長い系列

リポジトリの `codex/py-rainfall-infiltration-reference` ブランチのルートからPowerShellで実行する。
新規ファイルだけのdraftであり、マージしない。旧成果物を保存先に使わない。
以下の想定時間は目安。PCのCPU/メモリー依存で、未実測のPC時間を完了証拠にしない。

## 環境と短い参照（数秒〜1分）

Python 3.13を使用する。検証環境はPython 3.13.17、NumPy 2.5.4、SciPy 1.18.1。
新環境と出力はリポジトリ内の新しい`work/`に置く。

```powershell
$P = 'interim/2026-11-27-py-rainfall-infiltration-reference'
py -3.13 -m venv work/rainfall-venv
$Py = '.\work\rainfall-venv\Scripts\python.exe'
& $Py -m pip install -r "$P/requirements.lock.txt"
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$env:PYTHONUTF8 = '1'
$env:PYTHONDONTWRITEBYTECODE = '1'
& $Py "$P/checks.py"
if ($LASTEXITCODE -ne 0) { throw 'new checks failed' }
& $Py "$P/inspect_inputs.py" --root . --out work/rainfall-inputs.json
if ($LASTEXITCODE -ne 0) { throw 'input extraction failed' }
& $Py "$P/study.py" --out work/rainfall-small
if ($LASTEXITCODE -ne 0) { throw 'reference calculation failed' }
```

出力：`summary.json`、`reference-packet.json`、同CSVとmetadata JSON。
旧5項目は`rainfall-inputs.json`の`old_rows`。元記録からの抽出で新規本体求解ではない。
新規試験44件、skip/expected failureなし。次数・保存・解析式・契約拒否は参照実装の自己検査。
本体受入基準へ転用しない。`summary.json`ではCN時間半減差も併記する。
Linuxの最細空間で時間半減差/連続解差は約.000304（.0304%）。

元旧5項目の合否は、同条件の場合だけ次の**厳密な不等号**を用いる。
水頭RMS `<.003 m`、貯留絶対差 `<1e-4 m³`、隣接空間/時間差比 `<.8`、
独立FV累積水収支欠損 `<1e-9 m³`。
5量を比較する新しい契約の本体受入許容値は**未定義**で、accepted=null。
試験内の数値自己検査値をこの空欄へ代入しない。

## 長い1D・2D系列（PC専用、目安10秒〜数分）

```powershell
& $Py "$P/study.py" --pc-series --out work/rainfall-pc-1d
if ($LASTEXITCODE -ne 0) { throw 'PC 1D calculation failed' }
foreach ($n in 8,16,32) {
  & $Py "$P/pc_rainfall.py" --dim 2 --n $n --dt 100 --out "work/rainfall-pc-2d-n$n.json"
  if ($LASTEXITCODE -ne 0) { throw "2D n=$n failed" }
}
foreach ($dt in 100,50,25) {
  & $Py "$P/pc_rainfall.py" --dim 2 --n 16 --dt $dt --out "work/rainfall-pc-2d-n16-dt$dt.json"
  if ($LASTEXITCODE -ne 0) { throw "2D dt=$dt failed" }
}
```

1Dは64/128/256/512セルのGardnerと48/96/192セルVG類似柱。
2Dは帯状雨の類似問題で、**元の正方形3D雨ではない**。雨面積は奥行き1 mにつき.5 m²。
3期間の全水頭・貯留・流入出・全ステップの収支とNewton残差を保存する。
Newton失敗は`.failure.json`と非ゼロ終了。過去出力があるパスは拒否するので新しい保存先を使う。
次数・差は診断、最細値を厳密解にしない。未定義の受入値はnullを維持する。

## 元正方形雨の独立3D系列（PC専用、クラウド未実行）

`pc_rainfall.py` は自作のCartesian保存形FVであり、Rustや既存Pythonソルバを呼ばない。
3Dでは元の1 m立方体・中央.25 m²正方形・VG/Ss・H=0・3期間を扱う。
面の調和K、実貯留差BE、2d+1色の疎数値Jacobianと減衰Newtonを使う。
本体RT0と同一の離散化ではなく、独立な水理参照。
短い2Dで保存・反射対称性を検査したが、**3D実行・3D収束の検証は未確認**。
PCでの3D出力を得るまでは旧5行の新独立値は未定義のまま。

```powershell
foreach ($n in 8,16,32) {
  & $Py "$P/pc_rainfall.py" --user-pc --dim 3 --n $n --dt 100 --out "work/rainfall-pc-3d-n$n.json"
  if ($LASTEXITCODE -ne 0) { throw "3D n=$n failed or incomplete" }
}
foreach ($dt in 100,50,25) {
  & $Py "$P/pc_rainfall.py" --user-pc --dim 3 --n 16 --dt $dt --out "work/rainfall-pc-3d-n16-dt$dt.json"
  if ($LASTEXITCODE -ne 0) { throw "3D dt=$dt failed or incomplete" }
}
```

想定：8³は数秒〜数分、16³は数分〜数十分、32³は数十分〜数時間、メモリー数GBの余裕を推奨。
PC未実測の予算見積りで保証値ではない。長時間・メモリー不足なら中断し未達と記録する。
`complete=true`の新しい出力、実測経過時間、全10800 s到達を確認する。
Ctrl+Cや時間切れは成功ではなく、未出力・途中ファイルを合格にしない。
粗密の配列は上から下への深さ順で、元FV記録の底から上のz順とは異なる。
軸の並べ替え、元の明示三線形標本化、RT0の6tet平均、storage offset定義を
個別に監査してから元5量を再計測する必要がある。
本追加は旧RT0生結果の新しい取り込み/補間器を実装していない。
**独立3D参照だけの完了で、旧本体との差・元空間差比の合格は成立しない。**
本体4/8/16の元生出力との対応と元`compare.py`の定義を保った比較はPC残件。
本体の新規大規模求解や本体schemaの追加も本タスクでは行わない。

## 外部CSV/JSON比較

```powershell
& $Py "$P/compare.py" work/rainfall-small/reference-packet.json work/rainfall-small/reference-packet.csv --out work/rainfall-self-comparison.json
if ($LASTEXITCODE -ne 2) { throw 'expected undefined-acceptance diagnostic (exit 2)' }
```

これは参照の人工自己照合。`production_comparison=false`、差0でもaccepted=null。
将来の本体出力を外部で明示変換・監査した場合だけ、candidateを指定する。

```powershell
& $Py "$P/compare.py" work/rainfall-small/reference-packet.json work/production-audited.json --out work/rainfall-production-comparison.json
```

CSV列は `quantity,id,time_s,coordinate_m,value` の5列（余剰列も拒否）。
CSVと同じstemの`.metadata.json`、またはJSONの`metadata`と`rows`を必須とする。
例のmetadataを読み、同じ全パラメータ/モデル、時刻、量別ID/座標と照合する。
`source_kind=production_output`は**実本体値を監査済みの場合だけ**指定する。
例の値をそのまま本体値と名付けない。
head/fluxのcoordinateは深さ下向き、受け手effective_stressは標高上向きで、モデル欄に明示済み。
headは標高基準z=−xの全水頭、saturationはθ/θs、fluxは瞬時下向きm/s、
infiltrationは正味内向き累積m、effective_stressは圧縮正kPa。
0/300 sの表面流束は切替後の右極限、内部初期場は切替前の初期値。
座標・時刻の補間、自動単位変換、欠損補完、符号反転、クリップはしない。
終了2は未定義の正常診断、終了1は契約拒否。終了0による物理受入は実装しない。

## 既存Python全件（PC、目安10〜40分以上）

既存の固定依存一式も必要。証拠の部分取得状態ではアーカイブを元Git blobから
完全に取得する。署名/参照値/しきい値を変更しない。

```powershell
& $Py -m pip install -r requirements-reference.lock.txt
cargo build --locked --release --bins --examples
if ($LASTEXITCODE -ne 0) { throw 'existing CLI build failed' }
foreach ($group in 'python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence') {
  & $Py tools/run_ci_python_tests.py --group $group --output "work/rainfall-existing-$group"
  $code = $LASTEXITCODE
  Write-Host "$group exit=$code"
}
```

4群をすべて実行し、全IDの一意割当・件数・失敗/エラー/skipを確認する。
1群の失敗で後の群を省かない。Linuxで旧NC2群は元の厳密比較のsetUpClassエラーが既知。
各4件のメソッド未実行を成功/skipへ変換しない。Windowsのplatform差規則も元のまま。
既存全件が成功しなければ「全件成功」と書かない。
新先端Windows CIの成功は別途headと全群の実行・結論で確認する必要がある。
