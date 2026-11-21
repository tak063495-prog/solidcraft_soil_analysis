# NC行列平衡化・既存陰的接線実験と非ゼロSOIL源項流

2026-11-08。起点はPR #13の未マージブランチ `codex/plastic-spatial-convergence` の `b020b492ea091288dbfffac0514c570898f8c408`。このコミットを直接の親にして同じブランチへ積む。mainからの再実装ではない。PRをマージしない。追加クレジットは購入していない。

**低プレロードNCは最高0.9512977600097656で全撤去1.0未達、passed=false、production_adopted=false。SOIL体積源項による非ゼロ流れは、既存小ケース・弾性・閉鎖固着だけの限定合格。P3b/P4各行・構造/水/JOINT全体・P6全体は未完。**

生産Rust・既存封印保存点は変更せず、隔離した使い捨てworktreeで元の5 kPa入力から新規求解した。力1e-8・真の線形残差1e-10は緩めていない。100 kPaプレロード合格を転用していない。全解放・独立塑性3D全場・刻み感度・完了までの再開は未達で、合格や生産採用にしない。

追加の水検証は同じ54節点・12 Tet10・324材料GP・32接触GP、小さい湿潤柱・板2枚・梁1本だけを使う。`stages-joint-water-soil.json`のSOIL側a・構造側UpperClay・底面水頭2 m・外部水位0.75 mを引き継ぎ、既存Darcy材料入力 `volume_source_per_s` に、既存 `tests/darcy.rs::quadratic_source_exact_solution_and_outflow_sign` と同じ関係 `s=2Kzz` を与える。これは既存入力仕様による新しい小規模試験であり、新しい形状・構成則・接触仕様ではない。

下部土領域はz=0..0.5 m、Kzz=1e-4 m/s、s=2e-4 /s。底面H=2 m、接合面と側面は流量ゼロ。入力から解析すると、t=z-z0、L=0.5 mに対し `H=H0+s(Lt-t²/2)/Kzz`、`qz=s(t-L)`。したがって接合H=2.25 m、物理水圧17.5 kPa、底面流出1e-4 m³/s、接合流速0。流速は土領域内部で非ゼロである。源項がない前回のSOIL定水頭・流速ゼロや、PRESCRIBED固定貯水位の結果と区別する。

全節点変位・全材料GP応力/間隙水圧・JOINT物理面圧/有効接触力/両側荷重・板/梁の回転/設置参照/力/エネルギー・反力/モーメントを、入力のみから組み立てる独立3D連立系と元しきい値で照合する。Darcy水頭・全GP流速/圧力水頭・Dirichlet反力流出・源項/境界流出/水収支は独立FEMに加えて解析解と照合する。GL5→7、8工程途中の別CLI再開、意図した2棄却と棄却なしの同じ受理刻み列、VTU/CSV再読込を検査する。

限定範囲は飽和Darcy、MC弾性枝、閉鎖固着、同じ小部材のみ。一般塑性・開離/滑り・Richards・部分亀裂先端・接触BIRTH・高次板・弱結合・一般空間収束は未検証。P3b/P4各行、構造・水・JOINT全体、P6全体を完了にしない。

## NC実験の実測と未達

物理入力は既存生成関数 `make_free_contact.py` の `--mcc --pc 100 --preload-kpa 5 --force 1 .5 0` と一致する。数値制御は初期撤去1/128・最小1/4194304・Newton60・中央差分幅1e-9・dense LU。全実験は起点から分離したworktreeで初期状態から新規求解した。再分割で現在のNewton位置を保つ前回実験を共通に使う。陰的接線は既存 `implicit_stress_tangent_partition` を呼ぶもので、主値の戻り写像・構成則を置換しない。

|隔離実験|最高受理進捗|受理増分|最後の棄却残差|
|---|---:|---:|---:|
|行/列平衡化＋完全ピボット＋補償残差補正|0.9512968063354492|261|2.2485878254795473e-10|
|上記＋既存MCC陰的接線|0.9512972831726074|262|1.1157101737922449e-10|
|上記＋double-double完全ピボット、通常積/補償積の両ゲート|0.9512972831726074|262|1.3244034880452305e-10|
|double-double＋補償付き真の残差ゲート|0.9512977600097656|260|3.3328194891721366e-10|

**全実験で全撤去未達、passed=false、production_adopted=false。** 最後の経路でも前回最高からの差はわずか `4.76837158203125e-7` で、保持荷重は `0.048702239990234375`（約4.870224％）残る。停止理由は切戻し枯渇。180秒の受理境界予算による停止ではなく、破壊/安定限界/特異点の確定ではない。

通常積と補償積の差を隠さない。最後の実験は閾値を変えず、FMA積誤差＋補償加算で元の行列・右辺から残差を計算した。保存した**全1,213線形解**を、各floatを正確にDecimalへ移し70桁で独立再計算した。線形ゲートを通った1,206件の真の残差最大は `9.843326286208951e-11`。このうち67件は通常のf64行列積だけでは1e-10を超えたが、70桁の再計算は補償積と一致し、元の真の残差閾値以内だった。これは生産の通常積ゲートを変更したものではなく、未採用の隔離実験である。最後の棄却は70桁でも `3.3328194891721356e-10`、通常積は `6.904211393893965e-10`。実際の基準超過であり、内積丸め誤差だけを除けば完了するという結論は出さない。

最後の経路の受理済み線形報告の最大は `5.7368249997964867e-11`、受理済み力残差最大は `8.25858795864388e-9`。最終保存応力からの別再組立力残差は `2.9649446549173733e-10`、残存域λ=0は `7.710426636321882e-15`。撤去済み材料162 GPと接触32 GPの履歴凍結、全棄却での確定状態保持を照合した。保存場の監査は独立塑性変位解ではない。

double-doubleの2経路では終盤0.95129..0.95131の詳細診断も保存した。通常積/補償積両ゲートの経路は308 Newton場/803線探索候補、補償積ゲートの経路は322場/1,170候補を保存。点ごとの実分割数、全324材料GP状態、変位、実際の離散接線行列、線探索の受否を含む。保存された反復と候補の撤去済み材料凍結、保存された確定状態の接触凍結を検査した。全経路の全材料/接触履歴保存や、接触候補場を全部保存したという主張ではない。全解放と解放後の凍結は未到達。

入力のみからの既存独立3D MCC率式＋RK4参照も128/256刻みで再実行した。到達は `0.9509514227042946` / `0.9509514099202983`、最後の線形残差 `4.490642967739375e-10` / `3.582143726547472e-10` で停止。実測11.526/16.529秒で、300秒予算切れではない。0.9375での独立128→256変位相対L2差は `7.144051321864514e-6`。これは部分経路の刻み感度であり、新しい合格値や全経路収束ではない。全撤去までの本体/独立両実装の刻み感度、独立塑性3D全場、完了までの途中再開は引き続き未達。

## SOIL源項流の限定実測

8フレーム最大変位相対L2 `2.521224214760812e-11`、最大応力相対差 `3.6343001709696405e-12`、独立GL5→7変位差 `1.1311530276725194e-14`。最大値は共通初期フレームで前回SOILと同値になるが、旧記録の更新ではない。非ゼロ流れのflowフレームの変位差は `2.9584711241819295e-12`。全フレーム別の結果を `measurements.json` に分けて保存した。

解析解に対する最大節点水頭誤差 `6.661338147750939e-13 m`、最大GP流速誤差 `2.506404372243365e-16 m/s`。独立FEMとの流速差は `2.506810948058047e-16 m/s`。土側の水頭と水圧を構造側の間隙水圧ゼロと区別し、物理面圧/両側力/部材配置を検査した。8地点の別CLI再開、意図した2棄却と同一受理刻み列の棄却なし経路、VTU/CSV検査が通った。新規Python5試験は源項ゼロや流速/水頭/物理面圧/水収支/片側荷重/配置の改変、塑性・無指定SOIL・構造側水頭の不正入力を合格にしない。

初回追加試験の部材キー誤用1件と、解析監査で正規化後の境界キーを無視した比較1件は修正した。初回エラーログと修正後ログを生証拠で分離した。予備失敗を成功として隠さない。

## 証拠と再計測

`evidence.zip` には全4NC経路の入力・Run・進捗、実験ソース全体とビルドログ/実行版SHA、2経路の終盤詳細Newton場、最後の経路の全線形行列/右辺/解、独立RK4参照2経路の全到達場、SOIL源項の独立3D全場・8再開・棄却なし対照・全出力を保存する。`archive-sha256.json` にZIPと全エントリーSHAを記録する。巨大な2本のtraceはZIPのLZMA方式を使うため、Python標準 `zipfile` で展開する。実行ファイル本体は含めない。

新規展開先で全エントリーSHAを照合し、`remeasure.py` で入力照合・NC保存場/全線形70桁残差/保存反復履歴・独立部分刻み感度・SOIL独立全場/GL差/解析解/保存再開/VTUを再計算する。`measurements.json` と `remeasurement.json` は完全一致する。保存証拠の再計測であり、新しい求解や再開ではない。

```powershell
.venv/Scripts/python.exe -m zipfile -e interim/2026-11-08-nc-equilibration-and-soil-source/evidence.zip runs/NEW-soil-source-evidence
.venv/Scripts/python.exe interim/2026-11-08-nc-equilibration-and-soil-source/remeasure.py --evidence runs/NEW-soil-source-evidence --output runs/NEW-soil-source-remeasured.json
```

`base-preservation.json` に起点2,904追跡ファイルの全件照合を保存。うち2,902件は同一バイトで、変更した既存2件は独立検証PythonとCIだけ。旧interim 2,165件、全既存docs/examples/srcは不変。生産CLIのSHA256 `598a0d42a13bb00bdaffef52e82c45debb2a974c0f54c94cce9926b27fc4bfd6` も起点記録と一致する。

## WindowsユーザーPC手順

新先端Windows CIは未確認。Linux結果や旧先端の成功を転用しない。固定Rust toolchainと `requirements-reference.lock.txt` を満たしたWindows環境で、未存在の出力先へ実行する。

```powershell
git switch codex/plastic-spatial-convergence
git pull --ff-only
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release
cargo test --locked --release --all-targets
.venv/Scripts/python.exe -m unittest discover -s tests -p 'test_*.py' -v
.venv/Scripts/python.exe -m unittest discover -s tools -p 'test_*.py' -v
.venv/Scripts/python.exe tools/verify_wet_structures.py --solver target/release/soil-fem-3d.exe --output runs/NEW-soil-source-windows --water-mode soil-source
.venv/Scripts/python.exe tools/verify_soil_source_column.py --evidence runs/NEW-soil-source-windows --output runs/NEW-soil-source-analytic-windows.json
```

終了コードと先端SHA/OS/件数を保存し、コマンドごとに失敗を確認する。CIには同じSOIL源項ドライバと出力保存を追加したが、それ自体はWindows成功の証拠ではない。NC実験は生産経路へ採用しない。Linux保存点のWindows実行版への封印付替えや、新しい封印保存点の付替えをしない。

NCの数値実験をWindowsで再現する場合も、Linux保存点を再封印せず、以下のように保存ソースから別のWindows実行版をビルドして初期入力から新規求解する。コピー元の証拠は上の展開先を使い、出力先は未存在にする。これは合格手順ではなく、未達の実験の追試である。

```powershell
$Repo = (Get-Location).Path
$Evidence = Join-Path $Repo 'runs/NEW-soil-source-evidence'
$Scratch = Join-Path $Repo 'work/NEW-nc-accurate-windows'
git worktree add --detach $Scratch b020b492ea091288dbfffac0514c570898f8c408
if ($LASTEXITCODE -ne 0) { throw 'worktree failed' }
Copy-Item "$Evidence/nc/sources/accurate/stages.rs" "$Scratch/src/stages.rs"
Copy-Item "$Evidence/nc/sources/accurate/soil.rs" "$Scratch/src/soil.rs"
Copy-Item "$Evidence/nc/sources/accurate/nc_newton_trace.rs" "$Scratch/examples/nc_newton_trace.rs"
Push-Location $Scratch
try {
  cargo build --locked --release --example nc_newton_trace
  if ($LASTEXITCODE -ne 0) { throw 'build failed' }
  Get-FileHash target/release/examples/nc_newton_trace.exe -Algorithm SHA256
  New-Item -ItemType Directory work/NEW-linear-dumps | Out-Null
  $env:NC_LINEAR_DUMP = Join-Path (Get-Location) 'work/NEW-linear-dumps'
  & target/release/examples/nc_newton_trace.exe "$Evidence/nc/input" work/NEW-result
  if ($LASTEXITCODE -ne 0) { throw 'solver process failed' }
} finally {
  Remove-Item Env:NC_LINEAR_DUMP -ErrorAction SilentlyContinue
  Pop-Location
}
```

終了コード0だけでNCを合格にしない。`run.json` のcomplete/失敗理由・1.0到達・保持荷重・全残差と履歴を確認する。1.0に届いた場合にも、入力のみの独立塑性3D全経路と刻み感度、同一Windows実行版で作った未変更保存点からの完了までの再開一致、全解放後の凍結を揃えるまでは合格にしない。

## 保護記録

既存記録を変更しない。Fs区間 `[1.3250000000000002, 1.3531250000000001]` / `[1.3320312500000002, 1.3337890625000002]`、Griffiths中央面L2 15.582/15.886/15.593/15.672％の不合格、柱Tet4最細6.828％不合格/連続推定1.034％合格、既存MCC空間差、押出し応力/振動モニター/有限幅奥行き応力の判定不能、杭6.687→4.717％、旧降雨5項目未達、旧CONTINUOUS湿潤JOINT変位相対L2 `2.6845491097390954e-11` は元ファイルで保持する。前回SOIL定水頭の `2.521224214760812e-11`、固定PRESCRIBEDの `2.6845491097390954e-11` も別記録のまま保持する。

## 最終Linux全試験

fmt、clippy（`--locked --all-targets -- -D warnings`）、locked release build、全Rust **408件**・Python **304件**＋tools **2件**が成功した。新規はPython5件。失敗・無視・Rustフィルター除外・Pythonスキップは0。Rust全試験はビルドを含め1184.775秒、Python全試験313.502秒。全文ログ・終了コード・件数を `validation.json` と同名ログへ保存した。

試験前の724実行/試験/入力/既存文書/CIソースを `source-snapshot.json`/ZIPで固定し、終了後も全件同一だった。本interim内の記録・証拠処理スクリプトは別指紋 `final-files-sha256.json` に含める。生産CLIも起点記録と一致。全1,400生証拠を持つZIPの全エントリーSHAを新規展開で照合し、再計測JSONは完全一致した。

新先端Windows CIは未確認。CIの追加やLinux全試験成功をWindows成功、NC全撤去合格、P3b/P4/P6全体の完了へ読み替えない。
