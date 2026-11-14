# 構造・水・JOINTの限定追補：証拠と保存範囲

起点 `4b6c06bd9fb2a0ec0f87ae531bffb13bbf081b7e`、実装 `0151b24d94bc501d2a0d3729d80088ce8c805d10`、監査欠落拒否の補強 `5cbad531da9209f143533cfef3100b766823f6cc`。
同じPR #13・`codex/plastic-spatial-convergence`を更新した。mainだけを起点にせず、過去の変更を保持する。
[棚卸し・実装範囲・ユーザーPC手順](../../docs/STRUCTURE_WATER_JOINT_REMAINDERS.ja.md)を参照。

## 新しい限定証拠

項目名：**静水圧・定常Darcy CONTINUOUS JOINTと12 Tet10・板2枚・梁1本の独立3D弾性統合**。
54節点・12 Tet10・324地盤GP・32接触GP、8フレーム、施工の受理36/意図した棄却2（別に初期平衡の受理4）。
初期平衡・設置・水位上昇・部材載荷・定常Darcy・部材除荷・水位復帰・撤去を対象とする。
地盤はMC弾性枝、部材は線形、JOINTは閉鎖固着である。水頭だけを共有し、有効接触力と物理水圧を分ける。
旧部材例の幾何・荷重と既存湿潤柱を組み合わせた検証用入力で、現地の形状/強度同定ではない。

|指標|全8フレームの最大値（相対量は比率）|
|---|---:|
|全地盤節点3成分変位の独立相対L2差|2.6845491097390954e−11|
|工程増分自身で正規化した変位差|3.471490722721888e−11|
|全GP有効応力の独立相対L2差|3.8384581232696116e−12|
|独立GL5→7の変位相対差|1.233125503200593e−14|
|独立全体直接系の真残差|1.4670168596564533e−14|
|本体の全物理力相対残差|4.9732768816830085e−12|

最小有効接触圧は全工程を通じ85.29678194321733 kPa。独立解側でも摩擦円内部・MC降伏面内部と全体行列の全ランク/正定値を検査する。
全GPの応力/間隙比/水圧、全構造並進/回転/設置参照、全板生評価点・エネルギー、全梁端作用力/断面力、反力・力/モーメント、水頭を照合した。
入力から独立に生成した解であり、保存応力からの力再組立だけを独立3D変位解と呼んでいない。

8地点の別CLI再開は、初期平衡λ=0.25、設置0.125、水位上昇0.25、載荷0.125、Darcy移行0.25、除荷0.25、水位復帰0.25、撤去0.25。
最終確定状態・全フレーム・水理/施工監査・試行経路は測定時間`elapsed_seconds`だけを除いて完全一致する。
棄却のない同じ受理刻み列の新計算も全状態/フレーム/水理監査が一致する。
VTU/PVD/CSVの地盤GP・JOINT GP・板生評価点は既存の出力器で再読込して一致し、撤去時の空フレームを含む。

`small-verification.json`が数値要約と入力・ソルバー・参照コードのSHAを持つ。
`small-evidence.zip`は同じ検証の113ファイル（入力・結果・各ログ・出力・8再開保存点/結果・正常刻み列）を含む。
Linux実行版SHA256は `afc5b3f542bd2bd68b9520f5a392780cb23fbf8cdc630f689f71a6ca57dd260e`。
既存releaseバイナリの数値ソースは変更していない。再build後も同じSHAで、追加計測の実行版と一致する。
`source-snapshot.json/zip`は最終補強後の704ファイルを固定し、ZIP内の全ファイルと現行版をバイト照合する。
`small-verification.json`の一部パスは実行時の`work/structure-water/accepted-final`を指す。アーカイブ内のルートがそのフォルダー内容に対応し、`stage-results.json`は明示したSHAで照合できる。

## 既存記録の保存と限界

`inventory.json`は指定13文書＋P6一覧の起点版SHAと過去の`runs/`参照84件の配置状況を記録する。
今回のクローンに当該過去生ファイルは配置されておらず、歴史的な限定受入の内容は文書から棚卸しした。今回その生ファイルを再監査したとは主張しない。
起点に含まれるinterim内の証拠は削除・再計算・変更していない。
`preservation.json`は起点の2688追跡ファイル中2683ファイルが変更なしであることをGit差分で検査する。
変更した既存ファイルはCI定義と4文書のみ。4文書は元の本文バイトを保持し追補だけを加えた。
本体`src/`、既存試験、既存例、既存ツール、過去の参照値/失敗/数値証跡は変更なし。新規の例・検証器・試験・文書を追加している。

元Fs区間と6 Tet10の別Fs区間、Griffiths中央面4値の不合格と連続側判定不能、柱Tet4最細6.828％不合格/連続推定1.034％合格、MCC掘削空間差、押出し応力/振動モニター/有限幅奥行き応力の判定不能、杭6.687→4.717％、旧降雨5項目未達と元のしきい値を保持する。
実断面SRMは未実行。大規模・数時間級の指定求解、GUI/GPU、Native UIは行っていない。
P3b/P4のどの行も完了へ変えず、P6全体も未完了である。

## 試験とWindows CI

新規はRust2件・Python4件。実CLI・独立全体3D・8別プロセス再開に加え、水圧/接触力/地盤場/板設置参照の改変と水理監査欠落を拒否する。
未実装BIRTH/Richards・部分体積水没・摩擦円外の候補を、この限定参照器で合格にしない。これらの拒否を機能完成とは扱わない。
Linuxの最終全試験は、**Rust406件・Python269件＋tools 2件**が合格した。新規Rust2件/Python4件を含み、失敗・無視・スキップ0。
fmt、clippy（`--locked --all-targets -- -D warnings`）、locked release buildも合格した。
`local-validation.json`、`rust-tests.log`（60ターゲット、doc testの0件を含む）、`rust-test-list.log`、`python-tests.log`、`python-tools-tests.log`、`fmt.log`、`clippy.log`、`release-build.log`を対応付けた。
全Pythonは監査欠落拒否を含む最終補強後に新しく全実行し、固定704ソースとのバイト一致も確認した。途中版の全試験ログを最終版の全試験へ流用していない。
Windowsの[CI 42423617612](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/42423617612)は成功した。
Rust408件（Windows固有2件を含む）、Python269件＋tools 2件、fmt・clippy・release buildと小規模統合ドライバが合格し、失敗・無視・スキップ0。
実装は`5cbad531da9209f143533cfef3100b766823f6cc`、CIマージコミットは`0c6533dff888a42b1688de9ec0fc514fcdb37085`で、両Gitツリーは差分0。
`windows-verification.json`、`windows-job-status.json`、取得したログ全文`windows-job.log`、件数/測定の抽出`windows-test-summary.log`、ログ内JSONを復元した`windows-small-verification.json`へ対応付けた。
ログ全文はUTF-8、BOM除去・CRLF→LFで保存し、そのSHAを記録した。Windowsの追加ケースも全8フレーム・8再開地点・切戻し/出力が合格し、最大変位相対差は2.6845435035127625e−11だった。

artifact13595586491のメタデータ（450279 bytes、digest `sha256:35e7a1bf049f016c9da1993674e66acbcb85e38c7f8d21a0d5a07584b6fba152`）は`windows-artifacts.json`へ保存した。
コネクターはファイル参照を返したが、そのダウンロードURLがHTTP 403となったため、artifact内部のソース指紋/生結果バイトを直接照合したとは主張しない。
Gitツリーの一致、ログの704ソース凍結・全試験件数・小規模検証レポートは別に確認した。既存Fs追補のartifact取得限界も変更しない。
初稿0151b24のCI42423445660は後続pushでcancelledとなり、合格証拠に使わない（`windows-superseded-run.json`）。
最終の証拠追記はこのinterim内だけで、検証対象704ソースは実装版から変更していない。

実行コマンドはfmt、clippy `--locked --all-targets -- -D warnings`、`cargo build --locked --release`、`cargo test --locked --release`、`python -m unittest discover -s tests -p 'test_*.py' -v`、`python -m unittest discover -s tools -p 'test_*.py' -v`。
activeなRust全ターゲットと`tests/`/`tools/`のPython全試験が対象であり、`interim/`に固定保存した旧失敗版/研究試作を現行試験へ混ぜていない。

LinuxはRust1.99.1、Python3.13.17、固定9依存。詳細は`runtime.json`。
`retained-first.log`/`retained-second.log`/`retained-third.log`は検証器の初回実装での、未設置板の省略項目・地盤側へ記録する部材外力・水理監査の計測時間比較に対する対応漏れを保存する。
`retained-focused-rust.log`は新試験からenum材料状態のフィールドを直接参照できずコンパイルエラーになった記録で、試験側をシリアライズされた既存状態で照合するよう修正した。
これらを合格証拠へ流用していない。本体の数値経路や受入許容値の修正ではない。

## 継続する残件

低プレロードNC全上部撤去、部分亀裂先端、接触BIRTH、板/地盤の補間整合、Richards湿潤接触、2D互換JOINT、一般塑性の独立統合、SOIL/PRESCRIBED＋部材、実CAD堤体/掘削排水、一般地質自由水面、湛水/流出、杭列/群杭、梁塑性/座屈、大規模疎ブロック・山留め掘削は未完了。
ユーザーPCで揃えるべき出典/入力・独立証拠・再開/細分化・不合格の扱いは上記手順書に列挙した。形状や参照値を捏造しない。
一般空間収束、斜面5％、斜面の正確なFs、実断面SRM、別PC/GPU、完成配布、P6全体の未完了も継続する。
