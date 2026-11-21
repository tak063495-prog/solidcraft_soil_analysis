# 5 kPa NC全上部撤去：保存点再開と停止原因の計測

2026-11-08。指定先端 `8e772568f5b9f773398b984289463b0905156cae` を取得・チェックアウトし、同じ `codex/plastic-spatial-convergence` に追加した。mainからの作業ではない。PR #13は未マージのまま。

**最高受理進捗は 0.9512960910797119（受理266増分）。全上部撤去1.0は未達で、`measurements.json` は `passed=false`。破壊判定はしない。** 初期平衡・載荷・除荷の3フレームだけが完了している。100 kPa追加プレロードの合格を転用していない。

## 元の0.951保存点からの再開

前回 `native-evidence.zip` の `dense-1-128-min16384-v2/job` を新規作業領域へ展開した。最高受理点の封印保存点SHA256は `71f57e58ac1cd3da606a6f15ab85243392af9413c03c95b0faf25a505f648100`。再ビルドした生産CLIは元の封印と同じ `afc5b3f542bd2bd68b9520f5a392780cb23fbf8cdc630f689f71a6ca57dd260e` になった。同一実行版を展開コピーのsolverへ戻し、既存 `prepare-stages-recovery` → `run-stages-job` を実行した。

新しい受理は0件、0.95123291015625・合計259増分のまま、0.381秒で前回と同じNewton上限・切戻し枯渇となった。最後の棄却残差も `2.8009435309853623e-5`。時間予算の延長で進むとはしない。入力・保存点・封印・バイナリ結合・継続情報を付け替えていない。

## 別入力での既存数値制御の細分化

54節点・12 Tet10・324材料GP・32接触GP、p=pc=100 kPa、追加プレロード5 kPa、水平1/0.5 kNを固定した。形状・拘束・構成則・材料/接触局所刻み・摩擦・許容値・分割上限は元の生成値である。数値制御を変える既存保存点編集APIはないため、以下は同一の物理入力を最初から実行した別経路であり、元保存点から設定変更して進んだ結果とは呼ばない。

生産CLI・Rust数値ソースは起点と同一。既存dense LU、Newton60、撤去初期1/128を使用。元から許可された接線中央差分幅1e-7〜1e-9と最小刻みだけを細分した。`probe_low_preload_nc.py` に差分幅の明示オプションと入力照合を追加した。力基準1e-8・真の線形基準1e-10は緩めていない。

|最小刻み|中央差分幅|最高受理進捗|受理数|実測秒|最後の停止|
|---|---:|---:|---:|---:|---|
|1/262144|1e-7|0.9512939453125|263|65.460|Newton60、力残差2.4257404604674747e-6|
|1/262144|1e-9|0.9512939453125|262|63.064|Newton60、力残差1.956854326139585e-6|
|1/4194304|1e-9|0.9512960910797119|266|63.368|真の線形残差1.194127284855861e-10が元基準1e-10を超過|

最後の経路で棄却した試行先は0.951296329498291。その線形失敗直前に記録されたNewton力残差は `1.1383847722886682e-7` で、こちらも元基準1e-8を超える。これを線形失敗時の最終収束残差と混同しない。全経路は300秒の観測予算より先に数値停止し、時間切れではない。

## 切り分けで分かった範囲

- **接触**：撤去中の全32接触点は非活動。保存されたstatusはstickだが、現在の有効接触力の意味ではない。全接触履歴と撤去済み162材料GPは除荷完了時と完全一致。接触切替はこの終盤の直接の駆動要因として確認されなかった。
- **撤去荷重の配分**：残存域の保存荷重を、元の終端物理荷重＋(1−進捗)×元保持荷重と別計算で比較した。自由度相対誤差は全経路約2.18e-16。元のλ=0釣合い・接触消失項を含む保持荷重も既存監査で照合した。配分不整合は検出されなかった。
- **材料戻り写像**：前回の最後の受理増分0.951171875→0.95123291015625の全162活動点で、観測された最大分割10以内の全候補（1620条件）を既存Rust MCCで再評価した。局所エラー0件。点ごとの最も近い候補と保存状態の最大成分差は5.684341886080802e-14。点ごとの実分割数は保存されていないため、実際の分割列を完全復元したとは主張しない。
- **材料接線**：その候補で中央差分1e-7/1e-8/1e-9と既存の暗黙微分接線を比較した。最大相対差は1.2767494033124116e-10 / 4.5329994287435416e-10 / 4.0215881551489865e-9。全候補選択点で±摂動の塑性分割数切替は0。暗黙微分と元の積分の状態は完全一致し、入力状態も不変だった。この接線を生産CLIへ採用したわけではない。
- **全体接線の悪条件化**：保存場に対して、既存の独立連続体速度式で撤去方向の載荷/除荷枝を選び、残存域54自由度の接線を別組立した。156/162点が塑性載荷枝。最小固有値は0.431134→0.064702/0.083992→0.026503 kN/m、条件数は8.09万→53.74万/41.40万→131.15万へ悪化した。負の局所連続体接線を持つ点は80→81点。保存場での全体最小固有値は正であり、離散Newton接線の特異点や物理的な限界荷重を確定した証拠ではない。

悪条件化は実測で支持される。一方、失敗したNewton全反復の材料場・点別分割数・線探索候補は保存されていない。接線の枝混合・分割変更による再始動・降伏近傍の更新を失敗反復で完全に分離できておらず、停止原因を一つに断定しない。材料写像の成功を全体釣合いの成功とも扱わない。

これらは**本体保存場の診断**であり、入力だけから積分する独立塑性3D解ではない。`nc_material_probe` は同じ戻り写像を観測する読み取り専用exampleであり、独立参照・新構成則・新ソルバーではない。

## 最高受理点の履歴・釣合いと途中再開

最高点の保存応力からの別再組立力残差は `1.6058573425517926e-9`、同じ経路の残存域λ=0は `1.5244191234979123e-12`。元の力基準1e-8を満たす。棄却時の確定状態不変も検査した。ただし保持荷重の残り割合は `0.048703908920288086` で、全解放・全解放後の凍結は未確認。

同じ経路の封印済み265増分保存点（進捗0.9512958526611328）から、既存 `resume-stages` を別プロセスで実行した。266増分・最高点と後続の数値停止まで、初めからの経路の**Run全体**と完全一致した。除外は再帰的な実測時間 `elapsed_seconds` だけ。終了コード1も記録し、失敗の一致を合格にしない。`highest-partial-restart.json` を参照。撤去1.0までの再開一致ではない。

全撤去に未到達のため、今回の経路を独立塑性3D全場・刻み収束の受入判定へ進めていない。前回のRK4部分経路・場差は元ファイルのまま保持し、今回の最高点の参照や合格へ流用しない。全解放、独立全経路・両実装の刻み感度、完了までの途中再開が残る。P3b-04、P3b/P4全体、P6全体は未完。

## JOINT残件で追加した拒否試験

既存の同じ小さい湿潤柱・板2枚・梁1本だけを使用した。SOIL/PRESCRIBEDの水圧指定は既存の小ケース入力から取り、CONTINUOUS用の独立弾性統合参照が両条件を拒否することを確認した。塑性履歴を持つ種状態も弾性参照で合格にならない。Python試験2件（SOIL/PRESCRIBEDの2サブケースと塑性種状態）と `joint-unverified-gates.json` に記録した。

これは「証拠のない条件を既存参照で合格にしない」検査だけであり、SOIL/PRESCRIBED＋部材や塑性JOINT統合の完成証拠ではない。部分亀裂先端、接触BIRTH、Richards、2D互換、板補間、新形状・新構成則は作っていない。

## 生証拠とLinux検証

`measurements.json` に全経路の設定・到達点・全Newton残差列・線形残差・分割上限、保存場/接線診断、材料写像の点別比較を保存した。`evidence.zip` は入力・封印保存点・全試行・結果・ログ・材料写像の生データ・途中再開を含み、実行版は除外してSHAを保存した。`archive-sha256.json` でZIPと全エントリーを照合する。

Linuxのfmt・clippy（`--locked --all-targets -- -D warnings`）・locked release buildは合格。`cargo test --locked --release --all-targets` はRust407件が合格し、失敗・無視・フィルター除外0件。Python285件＋tools 2件もすべて合格し、失敗・スキップ0件。追加はRust1件・Python5件。全コマンドの終了コード0と実測時間を `validation.json` と各ログへ保存した。新規観測コードの開発途中に発生したPython試験の属性名衝突は修正済みで、最終全試験の結果と区別する。

試験前のソーススナップショットと最終照合を保存する。日本語文書は検証後の追補を含むため、実行/試験/入力/CIソースの不変照合とは区別する。`base-preservation.json` で指定起点の全追跡ファイルを照合し、旧interim・Rust数値ソース・保護記録を保持する。

## ユーザーPC手順

すべて新規出力先を使う。Windows実行版へLinux保存点の封印を付け替えない。Windowsは同じ5 kPa物理条件から再計算する。

```powershell
git merge-base --is-ancestor 8e772568f5b9f773398b984289463b0905156cae HEAD
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release
cargo build --locked --release --example nc_material_probe
cargo test --locked --release --all-targets
.venv/Scripts/python.exe -m unittest discover -s tests -p 'test_*.py' -v
.venv/Scripts/python.exe -m unittest discover -s tools -p 'test_*.py' -v

.venv/Scripts/python.exe tools/probe_low_preload_nc.py --solver target/release/soil-fem-3d.exe --output runs/NEW-nc-h9-fine --budget-seconds 300 --grace-seconds 5 --excavation-increment .0078125 --minimum-increment .0000002384185791015625 --dense-lu --tangent-perturbation 1e-9
# 上が終了した後、同一exe/入力/数値制御で保存点から再開。
.venv/Scripts/python.exe tools/probe_low_preload_nc.py --solver target/release/soil-fem-3d.exe --output runs/NEW-nc-h9-fine-recovery --resume-from runs/NEW-nc-h9-fine/job --budget-seconds 300 --grace-seconds 5 --excavation-increment .0078125 --minimum-increment .0000002384185791015625 --dense-lu --tangent-perturbation 1e-9
```

Linuxで元の0.951保存点を再開する場合は、前回の `native-evidence.zip` を新規フォルダへ展開し、`dense-1-128-min16384-v2/job/solver/soil-fem-3d.exe` に、封印と同じSHA256のLinux実行版を復元する。SHAが一致しなければ再開せず、新規入力から計算する。

```bash
python tools/probe_low_preload_nc.py --solver target/release/soil-fem-3d --output runs/NEW-original-resume --resume-from runs/NEW-previous/dense-1-128-min16384-v2/job --budget-seconds 180 --grace-seconds 5 --excavation-increment .0078125 --minimum-increment .00006103515625 --dense-lu
```

保存場診断は今回の `evidence.zip` を新規フォルダへ展開し、次で再計測できる。これは封印を変えて求解する処理ではない。

```bash
python tools/diagnose_low_preload_nc.py --previous-job runs/NEW-evidence/previous/dense-1-128-min16384-v2/job --profiles runs/NEW-evidence/original-resume runs/NEW-evidence/h7-min262144 runs/NEW-evidence/h9-min262144 runs/NEW-evidence/h9-min4194304 --material-probe target/release/examples/nc_material_probe --output runs/NEW-remeasured
```

数値停止後、同じ設定で時間だけ延ばしても到達の保証はない。次の改善では失敗反復の点別分割/降伏枝・材料残差・線探索と全体接線を計測し、元の力/線形基準を満たすか確認する。1.0に到達した場合だけ、独立塑性3Dの別積分・別刻み全場、全解放/凍結、λ=0、全工程履歴、完了までの再開一致を揃えて判定する。現在の独立参照コマンドと全工程途中再開コマンドは [既存手順](../../docs/LOW_PRELOAD_NC.ja.md) を参照する。

SOIL/PRESCRIBED＋部材の受入には、既存の各水圧モードについて物理水圧と有効接触力の別積分・両側荷重・部材回転/設置基準を含む独立3D全場、積分細分、保存再開、出力再読込を揃える必要がある。現在の `verify_wet_structures.py` は拒否するので、その終了を成功扱いにしない。塑性統合には入力のみからの塑性全経路が必要。新しい亀裂/BIRTH仕様は別途必要なままである。

## 保護記録

Fs区間 `[1.3250000000000002, 1.3531250000000001]` と `[1.3320312500000002, 1.3337890625000002]`、Griffiths中央面L2 15.582／15.886／15.593／15.672％不合格、柱Tet4最細6.828％不合格と連続推定1.034％合格、既存MCC空間差、押出し応力・振動モニター・有限幅奥行き応力の判定不能、杭6.687→4.717％、旧降雨5項目未達、湿潤JOINT変位相対L2 `2.6845491097390954e-11` は起点の元ファイルを変更しない。実断面SRM・大規模求解・GUI/GPU/Native UI・別PC実測は今回未実行。

今回の新先端に対するWindows CIは未確認。起点8e77256のWindows成功を新先端の試験結果へ転用しない。ZIP全88エントリーのSHAと再展開後の全計測JSON完全一致を確認済み。指定起点2781追跡ファイル中2778はバイト不変で、変更した既存3件は診断Python・拒否試験・元本文を接頭部保持した日本語文書だけである。
