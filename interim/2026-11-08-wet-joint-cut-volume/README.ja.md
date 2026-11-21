# 既存湿潤JOINTの部分水没体積積分の切り分け

2026-11-08。未マージPR #13の指定先端 `8397d04b569ce264bb939cad731ff373a0d930a3` を起点とする。同じ `codex/plastic-spatial-convergence` に直接の子コミットを1件追加する作業であり、mainから開始せず、PRはマージしない。

**今回の前進は、部分水没の差を水圧体積項へ定量的に切り分け、GL7にも残る誤差を解析積分で測れたこと。生産経路の修正・採用と湿潤JOINT全履歴の収束は未達で、`passed=false` を維持する。** 新しい形状・構成則・仕様を作らず、NCも再実行していない。

## 対象と再実行の範囲

既存の `examples/joint-water-column-tet10-mesh.json`、`ground-joint-water-mc.json`、`stages-joint-water-cycle.json` をバイト保持して使用した。54節点、12 Tet10、324地盤GP、32接触GP、MC弾性枝と既存Coulomb接触、CONTINUOUS静水圧、水位0.75 m、接合面圧2.5 kPaである。初期平衡とshear/open/recontact/reverse/reload/returnの7フレームを対象とする。

起点の一様64/128/256分割の本体・独立参照は再求解していない。前回の `../2026-11-08-wet-joint-separated-convergence/evidence.zip` と `archive-sha256.json` を読み取り専用の対照にし、アーカイブ全体と使用エントリーのSHA256、入力の一致を検査した。新しい求解は水圧体積項だけを解析積分へ置換した独立256分割・接触内8分割の1経路である。剛性・有効初期応力には従来の独立GL7、接触・地盤の出力GPには従来の演算子を使う。旧参照・旧生結果・旧計測JSONを上書きしない。

ソフトウェア試験では既存の小規模テストも実行するが、記録済み64/128/256収束系列やNCの再求解と区別する。開発中の同じ新積分の確認計算は `work/` に分け、ここに保存した最終新経路と重複加算しない。

## 解析積分の方法と適用限界

`tools/hydrostatic_cut_volume_reference.py` は検証用であり、Rust本体に接続していない。既存要素の中間節点が辺の中点であることを検査し、曲がったTet10は拒否する。水平水位面で四面体を切り、元の四面と水位面の多角形を内点から小四面体へ分割する。これは積分領域の内部的分割で、FEMの形状・節点・材料GP・接触GPを増やすメッシュではない。

各小四面体で `p=γw(H−z)` は一次、直線Tet10の `grad(N)` は一次なので、内部水圧力 `−∫grad(N)p dV` は二次多項式の積分になる。一次場a,bの積の恒等式 `∫ab dV=V[(Σa)(Σb)+Σ(ab)]/20` を用い、固定高次求積で水位の折れを跨がない。別方式として、分割領域上のGL3/GL5と既存形状関数微分でも照合した。元の未分割GL7を厳密解とは呼ばない。

全12要素で解析式対分割GL3/GL5の力成分差は最大 `8.881784197001252e-16 kN`。既存1 m角柱の解析量は湿潤体積0.75 m³、`∫p dV=2.8125 kPa m³`、`∫xp dV=∫yp dV=1.40625 kPa m⁴`、`∫zp dV=0.703125 kPa m⁴`であり、丸め誤差内で一致する。剛体並進・回転への仕事0と全アフィン変位場への仕事も試験した。乾燥、頂点/面が水位上、部分/全水没とTet4/Tet10の検査は、同じ既存柱を使った積分単体試験である。

切断領域そのものは、上下2高さの頂点を持つ直線四面体の中間水位で、下側頂点が1/2/3個なら湿潤体積比が1/8・1/2・7/8になる解析式でも確認した。下部の全水没要素を含む12要素別の照合を `element-domain-check.json` に保存し、全体量の一致だけに依存しない。

この多項式の正確さは直線要素・水平静水位の範囲である。曲面要素、Richards、塑性材料の任意経路、一般地形へ拡張していない。GP水圧・GP応力表示の値を解析力に合うよう校正していない。

## 水圧項の差と因果対照

以下の相対L2の分母は解析積分による組立済み水圧体積力で、測定値であり新しい合格基準ではない。GL11は力の積分だけであり、全経路の追加求解ではない。

|未分割の体積則|水圧力相対L2差|水圧力最大成分差 kN|
|---|---:|---:|
|入力GJ3（本体と同じ則）|0.002378863195959428|0.0033571074564643505|
|独立GL7|0.00032684856921443673|0.0004952190528565437|
|独立GL11|0.00004719499084616943|0.00006564137050463215|

一方、GJ3対GL7の組立剛性相対差は `7.452780472658479e-16`、有効初期応力力の相対差は `1.6691285703652517e-15`。部分水没でない下部要素の水圧力差も丸め誤差内である。剛性・有効初期応力・接合の面水圧を主原因にせず、上部の水位を跨ぐ水圧体積項へ局在できた。各要素の力ベクトルと全組立行列を生証拠へ保存した。

初期平衡についてのみ、入力から閉鎖固着の剛性を独立に組み立て、`(Kbulk+Kcontact)δu=−δIpore` を拘束空間で解いた。旧独立解と新独立解の初期接触が全点閉鎖・摩擦円内であることを別に検査してから、この線形補正を比較した。旧変位を補正量の生成に使っていない。

予測補正と、全非線形参照を実際に解いた初期変位差の食い違いは、GJ3対解析で最大 `7.569364349350246e-19 m`、GL7対解析で最大 `8.207621045899225e-19 m`。したがって初期変位差は、水圧体積項の差で説明できる。この対照は初期の閉鎖固着フレームだけのもので、開離/すべりの全履歴を線形解と認定しない。

新解析積分参照対旧独立GJ3の初期変位相対L2は `4.872971800166699e-4`、対旧独立GL7は `5.280761157772165e-5`。この表の分母は各旧参照変位である。旧GJ3対GL7の約 `4.947e-4` という記録は保持する。GL7側にも誤差があると分かったので、旧差をそのまま連続積分誤差とみなさない。

新経路の全1792小刻みは各工程 `i/256` を満たし、独立力残差は最大 `9.621468806620752e-11`、従来の `1e-10` 内である。全7フレーム・全変位・応力・反力・接触すべり/仕事/圧力と差を `measurements.json` に残す。MC弾性範囲と仕事恒等式も従来参照の条件を通過した。新経路の完了は全履歴収束の合格ではない。

## 収束・保護記録と残件

全体128→256の累積すべり1.70％・ベクトル22.4％という既存未達は変更しない。今回の解析体積積分参照について全体刻み感度を検証しておらず、既存履歴判定による実際の収束は確認していない。`history_converged=false`、`passed=false`、`production_adopted=false`。意図した棄却や再開一致を履歴合格の根拠にしていない。

生産側へ採用するには、既存の材料GP・静水位/定常場の移行・保持荷重・掘削・保存再開と矛盾しない水圧力の扱いを別途検証する必要がある。今回、RustのGP値を変更したり、この参照を直接生産関数へ差し込んだりしていない。全体刻み依存は別の残件として残り、単に体積項を直せば履歴収束するとは結論しない。

5 kPa NCの中点積分180秒・0.951171875は未完了であり成功/失敗判定へ変えない。Backward Eulerの0.9512960910797119と線形残差1.194e-10による切戻し枯渇も保持する。Newton/LU/許容/切戻しの再調整・100 kPaの結果の転用は行っていない。

Fs区間 `[1.3250000000000002, 1.3531250000000001]` と `[1.3320312500000002, 1.3337890625000002]`、Griffiths中央面L2 15.582／15.886／15.593／15.672％不合格、柱Tet4最細6.828％不合格と連続推定1.034％合格、既存MCC空間差、押出し応力・振動モニター・有限幅奥行き応力の判定不能、杭6.687→4.717％、旧降雨5項目未達、湿潤JOINT統合の変位相対L2 `2.6845491097390954e-11` を元ファイルのまま保持する。`base-preservation.json` に起点との差分と保存範囲を記録する。

今回(2)を切り分けて進められたため、代替のJOINT残件を追加で試していない。亀裂先端、接触BIRTH、高次板/弱結合、未存在入力の生成は行わない。大規模塑性4/8/16・8/16/32、MCC274〜14909、Griffiths1200/4800/19200、実断面/1200要素SRM、実CAD Richardsの再求解も行っていない。P3b/P4各行、構造・水・JOINT全体、P6全体は未完了のままである。

## 証拠とユーザーPC手順

`evidence.zip` は新入力のバイトコピー、新しい解析積分参照、項別の組立力/行列、所要時間の6エントリー。`archive-sha256.json` に全エントリーとZIPのSHA256を保存する。旧対照は起点にある前回ZIPへ依存し、今回ZIPへ二重収録していない。新規展開とSHA照合後、求解を呼ばない `--remeasure` による再計測が `measurements.json` と完全一致した（`archive-validation.json`）。これは新規求解や新しい履歴収束の成功ではない。

固定toolchainと `requirements-reference.lock.txt` を満たし、未存在の出力先を指定する。以下のWindows用手順は未実行であり、Windows CI成功を述べるものではない。

```powershell
git switch codex/plastic-spatial-convergence
git pull --ff-only
git merge-base --is-ancestor 8397d04b569ce264bb939cad731ff373a0d930a3 HEAD
if ($LASTEXITCODE -ne 0) { throw 'wrong history' }
git rev-parse HEAD
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release
cargo test --locked --release --all-targets
.venv/Scripts/python.exe -m unittest discover -s tests -p 'test_*.py' -v
.venv/Scripts/python.exe -m unittest discover -s tools -p 'test_*.py' -v
# 新しい解析積分経路だけ。記録済みGJ3/GL7の64/128/256を再求解しない。
.venv/Scripts/python.exe tools/diagnose_wet_joint_cut_volume.py --output runs/NEW-cut-volume
# 保存証拠だけの再計測。旧対照ZIPはリポジトリの元位置に必要。
.venv/Scripts/python.exe -m zipfile -e interim/2026-11-08-wet-joint-cut-volume/evidence.zip runs/NEW-cut-evidence
.venv/Scripts/python.exe tools/diagnose_wet_joint_cut_volume.py --remeasure runs/NEW-cut-evidence --output runs/NEW-cut-remeasured.json
```

Linuxでは環境の既存有効化後に同じコマンドを `python` で実行できる。時間切れ・未収束・欠落/非有限場・変わった刻み列は成功へ読み替えない。診断終了コード0と `diagnostic_completed=true` は計測完了を示すだけで、`passed=false` を真に変更しない。次にPCで履歴感度を調べる場合も、今回の参照と本体の体積項を区別し、既存判定規則・全場・実際の受理列を保持することが必要である。

Linuxの最終試験件数・ログ・固定ソース照合は本フォルダーの `validation.json`、各ログ、`source-snapshot.json` と `source-check.json` に保存する。


## 最終Linux検証

`cargo fmt --all --check`、`cargo clippy --locked --all-targets -- -D warnings`、`cargo build --locked --release`、`cargo test --locked --release --all-targets` がすべて終了コード0。Rust全408件、Python全324件（追加5件を含む）、tools全2件が成功した。失敗・Rust無視/除外・Pythonスキップは0件。Pythonの報告時間は363.940秒。コマンドと実測時間、件数と全文ログは `validation.json` と各ログに対応する。

最終検証の開始前に504ソースを固定し、全検証終了後に全件一致を確認した（`source-snapshot.json`、`source-check.json`）。生産CLIのSHA256は `598a0d42a13bb00bdaffef52e82c45debb2a974c0f54c94cce9926b27fc4bfd6` で、起点の記録と同じ。Rust64ファイルと既存入力138ファイルは不変である。

起点の追跡済み3002件中3001件を保持し、変更した既存ファイルは `tools/wet_joint_cycle_reference.py` のみ。旧interim2255件、既存docs90件、CI設定も保持した。新規追加は検証用積分器・診断器・5試験と本フォルダーの記録であり、生産側の変更はない。

`validation.json` の成功はLinuxソフトウェア検証についてのみ。数値受入の `measurements.json` は `passed=false` のまま。Windows CIは未確認、P3b/P4各行、構造・水・JOINT全体、P6全体も未完了である。
