# 低プレロードNCの反復位置保持実験とSOIL水圧＋既存部材の限定検証

2026-11-08。PR #13の指定先端 `9cd1b719e8dbb65f6d4dafc2f53df2df90679075` をfetchし、そのコミットから同じ `codex/plastic-spatial-convergence` に積んだ。mainから実装していない。PRはマージしない。過去のinterimと保護記録は変更しない。

**優先1の5 kPa低プレロードNCは撤去1.0未達、`passed=false`。今回の数値手順は0.9512958526611328・受理261増分で停止し、前回最高0.9512960910797119を改善しなかった。優先2は既存小ケースのSOIL＋部材について、弾性・閉鎖固着・既存の不透水構造側と流れ境界に限定して独立3D比較を通した。構造・水・JOINT全体の完了ではない。**

## 優先1：再分割時の位置保持と線形反復改良

前回 [Newton計測](../2026-11-08-nc-newton-and-prescribed-joint/README.ja.md) の結論を踏まえ、次の2操作を同時に一時適用した。

- 分割数を増やすとき、確定状態の構成則積分の起点はそのままで、全体Newtonの現在の反復変位を保持する。次反復で新しい固定分割の力残差と接線を再評価する。旧残差と新接線を混用しない。候補を未評価のまま受理する変更でもない。
- dense LUは同じ非対称行列の分解を使い、元の演算で再計算した真の線形残差が1e-10を超える場合、最大8回の残差補正を試す。改善しない補正は捨て、改善した場合だけ更新し、最後にも元基準を検査する。対角シフト、接線の対称化、収束基準の緩和を使っていない。

入力は既存 `make_free_contact.py --mcc --pc 100 --preload-kpa 5 --force 1 .5 0` の54節点・12 Tet10・324材料GP・32接触GP。前回標準経路と同じdense LU、Newton60、撤去初期1/128、最小1/4194304、中央差分幅1e-9である。力基準1e-8・真の線形基準1e-10、構成則・局所積分・形状・荷重を維持した。100 kPa追加プレロードの結果を転用していない。

`numerical-experiment.patch` に数値差分と、既存exampleの詳細観測を無効にして進捗だけ保存する差分を含めた。実験版の2ソース、ビルドログ、実行版SHA、入力、全進捗JSONLとRunを保存した。実行後はRustソースとexampleを起点と同じバイトへ戻し、生産CLIを改めてrelease buildした。今回この実験を生産実装へ採用しない。実験版の既存 `LinearReport.iterations` は反復改良回数を計測するよう変更しておらず、その回数の実測には使わない。失敗Newtonの詳細全場や局所反復全履歴を今回採取したとは主張しない。

最後の棄却先0.9512960910797119は真の線形残差2.2719550405006375e-10で停止した。切戻し枯渇であり、180秒の受理境界予算を理由とする停止ではない。未収束・切戻し枯渇を破壊・安定限界・合格にしない。今回の差分では再分割による位置喪失を避けても1.0へ進まず、悪条件化を解決したとはいえない。2変更を組み合わせた一経路なので、個々の寄与を切り分けたともいえない。

最高受理点の保存応力から別再組立した力残差は1.2420529790311404e-9、残存域λ=0は1.5244191234979123e-12。撤去済み材料・接触履歴は除荷終了時から凍結し、全棄却で確定状態が保たれている。これらは保存場の監査で、入力だけからの独立塑性変位解ではない。保持荷重の約4.8704%が残り、全解放と全解放後の凍結、独立塑性3D全場、両実装の全経路刻み感度、完了までの途中再開は未達。1.0未到達のためそれらの受入へ進めない。過去の部分経路・100 kPa経路を完了証拠にしない。

## 優先2：既存SOILの役割と流れ境界を引き継いだ独立3D

形状・材料・部材・工程は既存 `joint-water-column-tet10-mesh.json`、`ground-wet-plate-beam.json`、`stages-wet-plate-beam.json` の同じ54節点・12 Tet10・板2枚・梁1本を使う。既存 `stages-joint-water-soil.json` から次を一組として引き継ぐ。

- JOINTの `water={kind:soil, soil_side:a}` と不透水連続体側 `structure_bodies=['UpperClay']`。
- 既存SOILの定常Darcy入力全体（材料1、bottom水頭2 m、外部水位0.75 m）。構造側topの水頭3 mを残さない。構造側を帯水層にしない。

このSOIL流れ境界は土側bottom水頭2 mのみで、源項・流量境界もないため、独立解は土側の定水頭2 m・流速ゼロである。非ゼロSOIL浸透流の受入証拠にしない。

部材の取付位置・回転・設置参照、設置/載荷/除荷/撤去、静水圧上昇/復帰、初期のMC種状態は元の小ケースのまま。`UpperClay` の材料定数を架空の鋼材値等へ置き換えない。不透水連続体は既存実装どおり全応力基準であり、構造側の間隙水圧を全時点でゼロにする。地盤側水頭から接触面の物理水圧を独立に面積分し、有効接触ばねと分ける。

独立参照は入力のみから地盤・JOINT・板・梁の自由並進と構造回転を同じ直接3D連立系に保持する既存Python実装を拡張したもの。Rustの変位・応力・接線・設置参照で参照解を生成しない。MCの弾性範囲、接触の閉鎖固着、行列のランクと正定値を独立に検査する。定常Darcyも独立に組み立て、全水頭・定義域・流速・圧力水頭・Dirichlet流出と水収支を照合する。部材の全回転/設置参照、梁端力/断面力、板生評価点/エネルギー、全材料GP応力/空隙比/水圧、接触力/水圧/片側荷重、反力/モーメントを照合した。

| 計測 | 結果 |
|---|---:|
| 全節点3成分変位・最大相対L2 | 2.521224214760812e-11 |
| 工程変位増分・最大相対L2 | 3.3529909900857026e-11 |
| 全材料GP応力・最大相対L2 | 3.6343001709696405e-12 |
| 独立GL5→7変位相対差 | 1.1311530276725194e-14 |
| Darcy最大水頭差 | 0 m |
| Darcy最大流速差 | 0 m/s |
| Darcy最大圧力水頭差 | 2.220446049250313e-16 m |
| 独立Darcy水収支誤差 | 0 |

8フレームと、初期平衡・全7工程の途中8地点の別CLI再開を実行した。確定状態、全フレーム、試行とflow/excavation/transition監査は実測時間だけを除き一致した。意図した棄却2回と同じ受理刻み列の棄却なし経路の状態・フレーム・flow監査も一致。既存しきい値を変更せず、新しい合格数値を測定値から作っていない。独立全場は生証拠中の `soil/independent-3d-fields.json`、全比較は `soil-joint.json` を参照する。

既定CONTINUOUS参照はSOILを引き続き拒否する。明示 `--water-mode soil` が必要。構造側の宣言欠落・土側の逆指定・構造側top水頭境界、構造側への水圧混入、片側水圧荷重/物理水圧/部材設置参照/地盤変位/水頭定義域の改変、塑性種状態を合格にしない4試験を追加した。塑性、開離/滑り、Richards、固定水頭PRESCRIBED、接触BIRTH、部分亀裂先端、高次板、弱結合、一般空間収束や一般地質を受入対象へ拡大しない。

## 出力で発見したschema /14の扱い漏れ

最初のSOIL比較で生産求解の全場は参照と一致したが、既存 `high_order_io.py` は不透水構造の水圧ゼロと全応力基準フラグの処理をschema /13だけに適用していた。部材を含むschema /14の初期フレームでは、構造側GPに最大4.964701844302226 kPaの誤った水圧を出力し、`material_stress_reference` がなかった。これを合格にせず、同じ処理を/14へ適用する1行を修正した。生産求解の材料更新・荷重・接線は変更していない。

修正後、全8フレームのVTUを再読込し、入力から得た構造側フラグ・全GP水圧・全応力を独立参照と照合した。既存のGP/JOINT/板CSVとVTUの往復照合も通した。修正前GPファイルは生証拠 `regression/` に保存し、旧ソースのSHAと新ソースを区別する。

## 試験と生証拠

Linuxのfmt、clippy（`--locked --all-targets -- -D warnings`）、locked release build、元の診断exampleのrelease build、全Rust408件・Python294件＋tools 2件が成功した。追加はPython4件、失敗・無視・Rustフィルター除外・スキップ0。Rust全試験は1154.223秒、Python全試験は243.319秒。元ソースの診断実行版を再ビルドし、実験版と区別した。

全検査の件数・終了コード・実測時間は `validation.json`、全文は同名ログ。`source-snapshot.json` とZIPは検査前に720実行/試験/入力/文書/CIソースを固定した。終了後に同一バイトを照合する。収集範囲は既存 `tools/freeze_sources.py` の対象ディレクトリである。本interimの再計測・記録補助スクリプトは証拠の新規展開検査で別に実行し、最終マニフェストにも記録する。後続のこのinterimの文書・計測追記は別扱いであり、検証済みソースを変更しない。

`base-preservation.json` で、起点の2844追跡ファイル中2841のバイト不変、旧interim 2108件すべてのバイト不変を確認した。変更した既存3件は独立参照Python・出力Python・CIだけで、Rust数値ソース・examples・旧docsは不変。実験差分から保存ソースを完全再構成できることも確認した。

`evidence.zip` はNCの実験とSOILの入力・出力・全再開・独立全場・修正前出力を含む。実行ファイル本体は同梱せず `binary-sha256.json` に記録する。`archive-sha256.json` はZIP本体と全130エントリーのSHA256。Python `zipfile` で**新しいディレクトリ**へ展開し、`remeasure.py` で保存場監査・SOIL独立3D/積分差・保存された再開一致・出力を再検査できる。再計測は求解の再実行とは区別する。

```powershell
.venv/Scripts/python.exe -m zipfile -e interim/2026-11-08-nc-retained-iterate-and-soil/evidence.zip runs/NEW-soil-evidence
.venv/Scripts/python.exe interim/2026-11-08-nc-retained-iterate-and-soil/remeasure.py --evidence runs/NEW-soil-evidence --output runs/NEW-soil-remeasured.json
```

## ユーザーPC（Windows）手順

今回の新先端Windows CIは**未確認**。CIにはSOILの小規模ドライバと生証拠の保存を追加したが、追加しただけでは成功とはしない。過去先端のWindows成功やLinux結果を新先端へ転用しない。

`rust-toolchain.toml` のRustと、`.venv` の `requirements-reference.lock.txt` 固定依存を用意し、最新ブランチを取得してリポジトリルートで次を実行する。既存の出力・保存点を上書きしない。

```powershell
git switch codex/plastic-spatial-convergence
git pull --ff-only
powershell -ExecutionPolicy Bypass -File interim/2026-11-08-nc-retained-iterate-and-soil/windows.ps1 -OutputDir runs/NEW-soil-windows
```

この手順は指定起点の祖先確認、fmt/clippy/release build、全Rust/Python試験とSOIL独立3D/8再開を行う。Windows固有Rust試験をLinux件数へ混ぜず、手元の実測件数とOS/先端SHAを保存する。未収束・中断・入力拒否は合格にしない。

NCの数値手順を追試する場合は、元checkoutと別の**使い捨てworktree**へ `numerical-experiment.patch` を適用する。検証済み生産ソースを編集したまま使わない。以下でPythonは元checkoutのvenvを使う。Scratchと出力はすべて未存在のパスを選ぶ。

```powershell
$Repo = (Get-Location).Path
$Scratch = Join-Path $Repo 'work/NEW-nc-experiment'
$Python = Join-Path $Repo '.venv/Scripts/python.exe'
$Patch = Join-Path $Repo 'interim/2026-11-08-nc-retained-iterate-and-soil/numerical-experiment.patch'
git worktree add --detach $Scratch HEAD
if ($LASTEXITCODE -ne 0) { throw 'worktree failed' }
git -C $Scratch apply $Patch
if ($LASTEXITCODE -ne 0) { throw 'patch failed' }
Push-Location $Scratch
try {
  cargo build --locked --release --example nc_newton_trace
  if ($LASTEXITCODE -ne 0) { throw 'build failed' }
  New-Item -ItemType Directory -Force work/nc-input | Out-Null
  @'
import sys
from pathlib import Path
sys.path.insert(0, 'tools')
import probe_low_preload_nc as low
for name, value in zip(low.NAMES, low.expected_inputs(1/128, minimum_increment=1/4194304, dense_lu=True, tangent_perturbation=1e-9)):
    low.write(Path('work/nc-input')/(name+'.json'), value)
'@ | & $Python -
  if ($LASTEXITCODE -ne 0) { throw 'inputs failed' }
  & target/release/examples/nc_newton_trace.exe work/nc-input work/NEW-nc-result
  if ($LASTEXITCODE -ne 0) { throw 'experiment process failed; no acceptance' }
} finally { Pop-Location }
```

exampleは180秒の受理境界予算を持つ。終了コード0でもNC合格ではなく、`run.complete`、全撤去進捗、元の残差基準と全必要証拠を読む。低速PCで時間切れなら `passed=false` のまま記録する。WindowsへLinux保存点の実行版封印を付け替えない。元保存点再開は [既存NC文書](../../docs/LOW_PRELOAD_NC.ja.md) の契約に従う。1.0に達した場合に限り、撤去中全履歴、残存域λ=0、全解放/凍結、入力のみからの独立塑性3D全場、両実装の刻み感度、完了までの途中再開一致を揃える。

## 保護記録と残件

Fs区間 `[1.3250000000000002, 1.3531250000000001]` と `[1.3320312500000002, 1.3337890625000002]`、Griffiths中央面L2 15.582／15.886／15.593／15.672％不合格、柱Tet4最細6.828％不合格と連続推定1.034％合格、既存MCC空間差、押出し応力・振動モニター・有限幅奥行き応力の判定不能、杭6.687→4.717％、旧降雨5項目未達、旧CONTINUOUS湿潤JOINT変位相対L2 2.6845491097390954e-11は元ファイルのまま保持する。今回のSOIL数値で旧湿潤JOINT数値を置換しない。

P3b/P4各行、構造・水・JOINT全体、P6全体は未完。Griffiths Tet10細分、既存大規模MCC再求解、斜面SRM、82944降雨、105456性能、実CAD、長時間Richards、実断面SRM、Native UI、GPUは今回実行していない。全ソフトウェア試験の成功をこれらの受入へ読み替えない。
