# 有限幅の二層MCC掘削モデル（途中成果物）

SolidCraftの地盤10ボディと局所掘削用1ボディから、2種類のMCC材料を持つ274 Tet10・7398積分点・128板要素の5工程を解析しました。領域は10×4×6 m、掘削はx=3〜7、y=1〜3、z=3〜6 mです。底面8 m²、掘削量24 m³、四周の壁面積72 m²です。仮想材料による乾燥・排水条件の小ひずみ検証モデルです。

初期K0=0.6の応力場、pcと間隙比、プレロード、壁設置、除荷、掘削を渡します。材料属性の `pc_over_mean_effective` はpc/p′の比で、鉛直応力で定義する通常のOCRとは区別します。MCCの硬化・除荷履歴は全工程で保持します。

面積平均の掘削増分は、底面+5.312085 mm、周辺地表−0.032572 mmでした。施工開始からの表示累積変位は底面−0.774747 mm、周辺地表−17.248414 mmです。隆起増分と累積沈下を区別してください。粗い試験メッシュの数値であり、現地の設計値や空間収束を認定した結果ではありません。

掘削面の保持力を段階解放し、元の地表だけに上載圧を与えるため掘削底面の指定面圧はゼロです。3試行の棄却履歴復元、648撤去積分点の履歴凍結、掘削進行率0.03125からの保存再開を検査しました。別プロセス再開と連続実行は、時間記録を除く8群の全物理データが完全一致しました。新旧MCCソルバーも同じ未解析入力から全物理履歴が完全一致しています。

実SolidCraftの公式結果ビューアで変位・四周壁の曲げ・地盤18 m²断面・壁12 m²断面と生の節点/板値を確認しました。CAD変更による結果失効とUndoでの復帰、所有アプリの正常終了も確認済みです。PNGは実画面で、節点補間・板の生のDKT値等の表示規則は既存の結果ビューア文書を参照してください。

`finite-width-pit.solidcraft.json` を既存の `scripts/Start-SolidCraftSoil.ps1` のDesignへ指定して開きます。固定SolidCraftへの `integrations/solidcraft-projected/soil-fem.patch` を適用・ビルドした専用アプリと、本リポジトリのRustソルバーが必要です。公開名・別フォルダーへコピーした設計も別CLIで再表示し、11個のSTEPとground/stages/入力スナップショットが一致することを検査しました。アプリのSave Asで設計名を変えた場合は、解析属性のCAD対応を再確認して保存してください。

既にメッシュ化した未解析入力は、リポジトリ直下から以下で再計算できます。出力先は未作成のフォルダーを指定してください。

```powershell
cargo build --release --locked
python interim/2026-10-27-finite-width-pit/unpack_inputs.py runs/NEW-finite-pit-inputs
.\target\release\soil-fem-3d.exe prepare-stages-job runs/NEW-finite-pit-inputs/mesh.json runs/NEW-finite-pit-inputs/ground.json runs/NEW-finite-pit-inputs/case.json runs/NEW-finite-pit-job
.\target\release\soil-fem-3d.exe run-stages-job runs/NEW-finite-pit-job
```

`inputs.zip` には上記3入力とNative出力の11 STEPを収録しています。保存済み解・チェックポイントから初期状態を移したものではありません。中断は `stop-stages-job`、保存履歴からの別ジョブの準備は `prepare-stages-resume PREVIOUS_JOB NEW_JOB`、続いて `run-stages-job NEW_JOB` を使用します。出力・監査・独立参照の各CLIはリポジトリの既存文書を参照してください。

本版はRust263件・Python128件、fmt/clippy/buildと、小さい実CAD斜面の独立200→400分割による厳密な比較に合格しています。この有限幅モデル自体の独立100/200分割は計算中で、合格とは扱っていません。一般の空間収束も未完です。同一62ファイルのSolidCraft総合CIは合格済みで、無視試験とChrome未検出によるWebオフライン試験の省略は別の報告に保持します。

細分945要素モデルでは、詳細な板属性がNativeの64 KiB保存上限を超え、外部経路でも構造結合の記号行列サイズ上限により増分開始前に停止しました。失敗を残し、別作業用ソースで行列の全結合を保持する修正を検証中です。本公開版で解決済みとは扱いません。独立比較、細分化・新規PC配布、水/JOINT/アンカー等の全P0〜P6要件を引き続き進めます。
