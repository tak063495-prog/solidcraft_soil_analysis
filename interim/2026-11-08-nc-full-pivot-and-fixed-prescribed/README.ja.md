# 完全ピボットNC実験と固定水頭PRESCRIBED＋小部材の限定検証

2026-11-08。指定先端`9a15c83cd7de73e42e9344a8b4bfab31918e5d49`を取得し、そのコミットから同じ`codex/plastic-spatial-convergence`で続行した。PR #13はマージしない。**低プレロードNCは全撤去1.0未達、passed=false。固定水頭PRESCRIBEDは既存小ケース・弾性・閉鎖固着だけの限定合格。P3b/P4各行・構造/水/JOINT全体・P6全体は未完。**

[条件・実測・未達をまとめた日本語文書](../../docs/NC_FULL_PIVOT_FIXED_PRESCRIBED.ja.md)、[全計測JSON](measurements.json)、[固定貯水位の入力由来1D解析対照](analytic-reservoir.json)を参照する。

## 優先1：全撤去には届かない

既存の54節点・12 Tet10・324材料GP・32接触GP、`--mcc --pc 100 --preload-kpa 5 --force 1 .5 0`のみを使った。既存生成関数と3入力を完全一致で照合。初期撤去刻み1/128、最小1/4194304、Newton60、中央差分幅1e-9、力1e-8・真の線形残差1e-10を維持し、100 kPaプレロード合格を転用しない。

前回の未採用実験の部分ピボットLUを完全ピボットLUへ変更した。再分割時のNewton現在位置保持と最大8回の単調な残差補正は前回と同じである。元の行列積による真の残差を検査し、基準緩和や対角シフトは使わない。`numerical-experiment.patch`はこの数値差分と、詳細全場診断を無効にして進捗だけを残すexample差分を含む。補正回数や失敗Newtonの詳細全場を測定したとは主張しない。

|項目|今回の実測|
|---|---:|
|最高受理撤去進捗|0.9512972831726074|
|受理増分|262|
|前回最高との差|0.0000011920928955078125|
|最後の棄却先|0.9512975215911865|
|最後の真の線形残差|1.5814795527509818e-10|
|保存場の別再組立力残差|9.955027653861772e-11|
|残存域λ=0残差|1.5255013280638015e-12|
|保持荷重の残存割合|0.04870271682739258|

停止理由は切戻し枯渇であり、180秒の受理境界予算切れではない。全棄却で確定状態不変、除荷完了時から撤去済み材料/接触履歴は凍結。保存場の監査を独立塑性3D変位解にしない。全解放/その後の凍結、入力のみからの独立塑性3D全場、両実装の全経路刻み感度、撤去完了までの途中再開は未達である。

わずかな部分経路改善だけでは採用条件を満たさず、実験を生産CLIへ採用していない。生産Rustソース・exampleは指定起点と同じバイトへ戻し、release CLIを再ビルドした。新しい封印保存点を作っていない。旧封印・入力・実行版封印、WindowsへのLinux保存点封印付替えは行わない。破壊/物理的限界/離散接線特異点の結論は出さない。

## 優先2：既存固定貯水位仕様と小部材

`joint-water-column-tet10-mesh.json`、`ground-wet-plate-beam.json`、`stages-wet-plate-beam.json`の既存小ケースを使う。水圧ポリシーを既存PRESCRIBEDへ替え、`level_m`は元地盤入力の初期水位1 mから得る。固定貯水位に競合する既存rise（1.25 m）だけを除外し、その他の設置・載荷・定常Darcy・除荷・1 m水位復帰・部材撤去を維持する。形状・材料・構成則・接触法則・部材位置を変更しない。

`--water-mode prescribed-fixed`と`fixed_prescribed=True`の明示を要求する。既存のCONTINUOUS/PRESCRIBED追従水位/SOIL参照の既定範囲を広げない。固定値と初期水位/静水位の競合、塑性種状態は参照で拒否し、競合水位は本体でも入力拒否される。拒否試験は物理合格と区別する。

入力だけから地盤・JOINT・板・梁の自由並進/回転を直接3D連立系へ保持した独立参照を使う。全7フレームの全節点変位・工程増分・全GP応力/空隙比/水圧・接触力と物理水圧・両側荷重・全構造回転と設置参照・梁力・板生評価点とエネルギー・反力/モーメントを既存しきい値で照合した。GL5→7も照合。Rustの変位・応力・接線・設置参照を参照生成に使わない。弾性、閉鎖固着、全体ランク/正定値を独立解側で検査する。

定常Darcyは接合固定水頭1 m、外部水位0.75 m、bottom/top水頭2/3 mを区別し、水頭・全GP流速・圧力水頭・Dirichlet流出・水収支を独立FEMで照合した。さらに同じ入力から1D解析解を作り、下部の流速+0.0002 m/s、上部-0.0008 m/s、接合物理水圧5 kPaを照合した。これはPRESCRIBEDの非ゼロ流れで、前回SOILの非ゼロ浸透流を合格にしたものではない。

変位最大相対L2は2.6845491097390954e-11、増分3.471490722721888e-11、GP応力3.8384581232696116e-12、GL5→7変位差1.1630766646304642e-14。変位最大値は共通の初期フレームで旧CONTINUOUS保護値と同じになったが、別に実行した今回の結果であり、旧記録を上書きしない。各フレームは`measurements.json`に分離して保存した。

7地点（初期平衡と全6工程）の別CLI途中再開、意図した2棄却と同じ受理刻み列の棄却なし対照、VTU/生GP/JOINT/板CSVの再読込が通った。追加Python5試験には改変した流速・水頭・片側荷重・物理面圧・設置参照の検出も含む。予備試験のフィールドキー誤用2件は修正し、エラーログと修正後ログを生証拠に分けて保存した。予備失敗を最終成功として隠さない。

## 生証拠と再計測

`evidence.zip`は124エントリー。`archive-sha256.json`にZIPと全エントリーのSHA256を記録した。NCの3入力・進捗JSONL・Run・実験ソース・元ソース・差分・ビルドログ・実行版SHA、固定PRESCRIBEDの入力/Run/独立全場/7再開/棄却なし対照/出力/入力拒否証拠を含む。実行ファイル本体は同梱しない。実験差分を別の新規ディレクトリの元ソースに適用し、保存実験ソースと完全一致することを確認した。

新規展開先で全エントリーSHAを照合し、`remeasure.py`によるNC保存場監査・入力のみの独立全場/GL差・保存再開照合・VTU水圧/全応力の再検査を行った。`remeasurement.json`は保存証拠の再計測であり、求解再実行とは区別する。

```powershell
.venv/Scripts/python.exe -m zipfile -e interim/2026-11-08-nc-full-pivot-and-fixed-prescribed/evidence.zip runs/NEW-fixed-evidence
.venv/Scripts/python.exe interim/2026-11-08-nc-full-pivot-and-fixed-prescribed/remeasure.py --evidence runs/NEW-fixed-evidence --output runs/NEW-fixed-remeasured.json
.venv/Scripts/python.exe interim/2026-11-08-nc-full-pivot-and-fixed-prescribed/audit-fixed-reservoir.py --evidence runs/NEW-fixed-evidence/fixed --output runs/NEW-fixed-analytic.json
```

## Linux全試験

最終検査の件数・終了コード・実測時間は`validation.json`、全文は同名ログへ保存する。fmt、clippy（`--locked --all-targets -- -D warnings`）、locked release buildと、全Rust（`--all-targets`）/Python・tools試験を実行する。試験前の722実行/試験/入力/文書/CIソースを`source-snapshot.json`/ZIPに固定し、終了後の同一バイトを検査する。本interim内の証拠処理スクリプトは別に実行確認し最終ファイル指紋へ含める。後続の本interim内文書/JSON追記は、固定ソースとは区別する。

## WindowsユーザーPC手順

新先端Windows CIは未確認。CIの追加やLinux成功だけでWindows成功にしない。OS/先端SHA/件数をそれぞれ保存し、過去先端CIやWindows固有試験件数を混用しない。

リポジトリの固定Rustと、`requirements-reference.lock.txt`を入れたWindows `.venv`を用意し、未存在の出力先を指定する。

```powershell
git switch codex/plastic-spatial-convergence
git pull --ff-only
powershell -ExecutionPolicy Bypass -File interim/2026-11-08-nc-full-pivot-and-fixed-prescribed/windows.ps1 -OutputDir runs/NEW-fixed-windows
```

NC実験は元checkoutと別の使い捨てworktreeで行う。Linux保存点を使わず、新しいWindows実行版から同じ入力を新規求解する。既存出力を上書きしない。

```powershell
$Repo = (Get-Location).Path
$Scratch = Join-Path $Repo 'work/NEW-nc-full-pivot'
$Python = Join-Path $Repo '.venv/Scripts/python.exe'
$Patch = Join-Path $Repo 'interim/2026-11-08-nc-full-pivot-and-fixed-prescribed/numerical-experiment.patch'
git worktree add --detach $Scratch HEAD
if ($LASTEXITCODE -ne 0) { throw 'worktree failed' }
git -C $Scratch apply $Patch
if ($LASTEXITCODE -ne 0) { throw 'patch failed' }
Push-Location $Scratch
try {
  cargo build --locked --release --example nc_newton_trace
  if ($LASTEXITCODE -ne 0) { throw 'build failed' }
  New-Item -ItemType Directory work/nc-input | Out-Null
  @'
import sys
from pathlib import Path
sys.path.insert(0, 'tools')
import probe_low_preload_nc as low
for name, value in zip(low.NAMES, low.expected_inputs(1/128, minimum_increment=1/4194304, dense_lu=True, tangent_perturbation=1e-9)):
    low.write(Path('work/nc-input')/(name+'.json'), value)
'@ | & $Python -
  if ($LASTEXITCODE -ne 0) { throw 'inputs failed' }
  Get-FileHash target/release/examples/nc_newton_trace.exe -Algorithm SHA256
  & target/release/examples/nc_newton_trace.exe work/nc-input work/NEW-nc-result
  if ($LASTEXITCODE -ne 0) { throw 'process failed; no acceptance' }
} finally { Pop-Location }
```

exampleは180秒の受理境界予算を持つ。終了コード0でもNC合格ではない。`complete`、1.0到達、元の残差基準、保持荷重ゼロと全必要証拠を読む。未収束・時間切れ・入力拒否を破壊/合格にしない。1.0到達時だけ、全解放とその後の凍結、独立塑性3D全場、両実装の全経路刻み感度、完了までの途中再開を揃える。

## 保護と未達

起点2874追跡ファイル中2872が同一バイト。変更した既存2件は独立参照PythonとCIのみ。旧interim 2137件、全既存docs、Rust、examples、全保護記録は不変（`base-preservation.json`）。Fs両区間、Griffiths中央面L2全4値の不合格、柱Tet4最細不合格/連続推定合格、既存MCC空間差、押出し応力/振動モニター/有限幅奥行き応力の判定不能、杭、旧降雨5項目未達、旧CONTINUOUS湿潤JOINT値を変更しない。

SOIL非ゼロ浸透流、一般塑性/開離/すべり、Richards、接触BIRTH、部分亀裂先端、高次板、弱結合、一般空間収束は未達。部分亀裂先端/BIRTH/高次板/弱結合は仕様と独立証拠が足りないまま受け入れない。形状や構成則を捏造せず、ユーザーPCで検証可能な入力仕様と独立参照が揃うまで未完とする。P3b/P4各行、構造・水・JOINT全体、P6全体を完了にしない。指定された大規模再求解・実CAD本計算・長時間Richards・実断面SRM・Native UI・別PC/GPU・完成配布は行っていない。

NCの受理済み収束ゲートを`nc-convergence-gates.json`に保存した。力残差最大8.921514504787992e-9、受理済み641線形求解報告の真の残差最大7.732861802999355e-11は、それぞれ元基準1e-8/1e-10以内である。最後の棄却をこの最大値から除外して合格へ見せるものではなく、棄却残差と`passed=false`は同JSONに併記した。

## 最終Linux実測

fmt、clippy（`--locked --all-targets -- -D warnings`）、locked release build、全Rust **408件**・Python **299件**＋tools **2件**が成功した。失敗・無視・Rustフィルター除外・Pythonスキップは0。追加はPython5件で、Rust数値ソースは変更していない。

Rust全試験は1131.181秒（ビルドを含む）、Python全試験は247.651秒。試験前後の722固定ソースと生産CLIは完全一致した。生証拠ZIPの新規展開と再計測も成功した。新先端Windows CIは未確認であり、Linux全試験成功をNC全撤去やWindows成功へ読み替えない。
