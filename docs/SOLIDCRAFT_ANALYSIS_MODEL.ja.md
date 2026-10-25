# SolidCraft内の解析属性とネイティブSTEP受渡し

2026-10-24。専用作業ツリー`external/solidcraft-soil-fem`で実装中。**保存/改名/Undo・実GUI出力・実STEP→MCCまで限定検証したが、P5全体の完了を認定した記録ではない。** [実ネイティブ検証と修正版](SOLIDCRAFT_NATIVE_VERIFICATION.ja.md)を参照。engine6試験、取込の当初5試験と静的品質検査は[初回時点の記録](../runs/solidcraft-analysis-quality-progress-20261024-v2.json)を保持する。現在は取込6件・表示受渡し7件に加え、UI3試験・clippy・書式・WASMと実画面の小数表示を確認済み。[外部ジョブワーカー](SOLIDCRAFT_JOBS.ja.md)は5試験と実ジョブの成功/停止を限定検証した。WASMの未変更3ファイルには警告がある。nativeのジョブ/結果UIは引き続き未完。

## 保存するものと実行境界

`metadata.org.solidcraft.soil-fem`へ、版付き地盤入力・施工入力、解析ボディ名とCADボディの対応、作成フィーチャー・コンポーネント・出力ファミリー、確認したCADのSHA256を保存する。材料の密度・色・CAD表示材から土の物性を推定しない。K0、OCR/pc、間隙比などの意味は既存FEM入力の定義を保持する。

engineのコマンドは次のとおり。プログラム呼出しでダイアログを開かない。

|コマンド|動作|
|---|---|
|`analysis.soil.model`|属性・CAD候補・対応不成立・変更検出と版ハッシュを読む|
|`analysis.soil.read_inputs`|許可された地盤/施工JSONを読み、受渡し契約を検査する。設計は変更しない|
|`analysis.soil.configure`|完全なボディ対応と読込時CAD/属性ハッシュを照合し、1回のUndoで保存する|
|`analysis.soil.snapshot`|対応が成立し、CAD変更を確認済みの場合にだけ入力スナップショットを取得する|
|`analysis.soil.export`|新規フォルダーへ各ボディのSTEP、入力、元設計、出典ハッシュを出力する。既存フォルダーを上書きしない|

保存後のCAD改名は、同じUndo操作の中でCAD名の対応を更新する。解析ボディ名は保持する。作成フィーチャーの出力が増減・変更した場合、名前が同じ別フィーチャーへ置換された場合や、ボディが消えた場合は未対応となる。境界箱や体積の似たボディへ自動で割り当てない。

編集中にCADまたは保存属性が変わると、古いハッシュによる保存を拒否する。CADハッシュは解析メタデータと表示だけの可視性を除き、他の設計内容を保守的に含める。したがって、幾何以外の設計属性変更でも再確認が必要な場合がある。未確定の移動、編集中のスケッチ、評価エラーがある状態では出力しない。コンポーネントの配置インスタンスを持つアセンブリは、この版の取込対象外である。

属性保存はメタデータの64 KiB上限を使い、解析ボディ/材料/工程は各256まで。大規模な積分点初期状態を外部参照として保持する仕組みは後続項目で、上限を超えたデータを切り捨てて保存しない。

## 画面の構成

`Soil Analysis`コマンドの対話開始で、Ground / Construction / CAD bodiesの編集窓を開く。既存の地盤・施工JSONを読み、物性/初期状態と施工の値を編集し、各解析ボディへCADボディを明示する。保存するまで設計属性は変えない。ファイル選択は既存の非同期ホストサービスを使い、選択中もイベントループへ戻る。

この最初の編集窓では、読み込んだ入力の値を編集する。全材料・工程の新規作成、列挙方式の切替、任意CAD面の永続対応、ジョブの起動・監視・中断と結果コンター表示は後続の実装である。UIの作成だけをP5完了とは扱わない。

## 新しいSolidCraftから既存FEMへ渡す

新しい設計属性には専用のpending format noteを追加した。旧portable版で新形式の元設計を読み直してはいけない。ネイティブSolidCraftが各ボディを直接STEPへ出力し、FEM側では次のファイルを検査する。

- `ground.json`、`stages.json`：明示した解析入力。
- `body-XXXX.step`：ボディ別のmm単位形状。論理名は対応表から取得する。
- `input-snapshot.json`、`snapshot-canonical.json`：同じ入力の表示用/ハッシュ用表現。Rust/Pythonの浮動小数点の文字列表現差を避けるため、SHA256はRustが保存した正規化バイト列に対して検査する。
- `source.solidcraft.json`：元設計の保存記録。取込アダプターは実行しない。
- `capture.json`：全ファイル・入力・出力元実行版のハッシュと完了判定。最後に書き、部分出力と区別する。

取込時はファイル名、フォルダー外参照、全ハッシュ、入力/対応表の一致と数値型を検査する。Gmshで各STEPが単一の閉ソリッドであること、体積・mm→mの換算を独立に調べる。FEMの工程パラメータや材料の妥当性は、メッシュ作成後のRustソルバーが判定する。画面で属性を保存しただけでは解析の成立を認定しない。

```powershell
# SolidCraftのExport analysis inputsで作った新規captureフォルダーを指定する
pwsh -NoProfile -File scripts/Invoke-SolidCraftCapture.ps1 `
  -CaptureFolder D:\Work\MY-NEW-CAPTURE -MeshSizeM 2 -ElementOrder 2
# JOINTがある場合は、変位節点分離の対応を明示する
pwsh -NoProfile -File scripts/Invoke-SolidCraftCapture.ps1 `
  -CaptureFolder D:\Work\MY-OTHER-CAPTURE `
  -SplitFile examples/joint-column-split.json -MeshSizeM 4
```

FEMは既存の共形メッシュ、履歴JOINT、施工解放、保存再開、独立再組立、VTU/PVD再読込へ進む。Python/Gmshは外部FEMプロジェクトのメッシャー境界で使い、SolidCraftアプリへリンク/組込みしない。メッシュサイズ等はこの外部ジョブの指定で、CAD属性スナップショットとは別にメッシュ出典へ記録する。

現在の外部コントローラーは、その後に[工程別表示データとVTU照合](SOLIDCRAFT_RESULTS.ja.md)も作成する。これは外部受渡し基盤で、SolidCraft内への描画は未接続。取込側には指数オーバーフローの拒否を追加し、現在の試験数は6件（`runs/native-capture-import-tests-20261024-v3.log`）。冒頭の品質検査記録はその前の固定時点の5件である。

## 検証の区別

取込アダプターの5試験は、別のCAD表示名と解析名の対応、400 m³の既知形状、単位換算、非有限JSON、数値1とboolのすり替え、STEP改変、不正パス、未完了出力、既存出力保全を対象とする。実SolidCraftの画面/コマンドからFEMまでの統合試験は、この単体試験とは別に行う。

既存の[湿潤JOINT固定版206試験](../runs/wet-joint-acceptance-20261024-v1/final-validation.json)は、その凍結ソースと解析バイナリの証拠である。今回追加したSolidCraft拡張・UI・取込経路の受入を代用しない。P6全体の未達要件は[P6一覧](P6_STATUS.ja.md)を保持する。
