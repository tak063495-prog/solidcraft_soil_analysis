# CSR疎求解と施工途中の保存再開

MCC全体Newtonの接線をCSRへ直接組み立て、非対称GMRES＋ILU(0)で解く。施工の受理増分境界をチェックポイントとして保存し、同じ入力・ソルバーバイナリから別プロセスで再開できる。

## 実行

```powershell
cd D:\Work\soil-fem-3d
.\scripts\Invoke-LocalExcavation.ps1 -MeshSizeM 0.75
.\scripts\Test-SparseAndRestart.ps1
```

局所掘削スクリプトの既定は`examples/stages-local-sparse.json`。`examples/stages-local.json`は密LU比較用。以前の段階解析入力は、`linear_solver`を省略した場合に密LUを使う。

```json
"maximum_free_dofs": 100000,
"linear_solver": {
  "method": "sparse_gmres_ilu0",
  "options": {"restart":60,"maximum_iterations":3000,"pivot_tolerance":1e-14}
}
```

設定は`settings`内へ置く。密LUは`{"method":"dense_lu"}`。CSRの設定上限は300000自由度、密LUは2000自由度。これは入力上限であり、全材料・全工程でその規模まで検証済みという意味ではない。

## 組立と収束判定

節点支持・平面拘束・有効要素を反映した自由度の接続から、CSRの行ポインタ・列番号を作る。数値組立は同じ疎パターンへ加算し、密行列を経由しない。平面拘束で同じ自由度へ集約する成分も加算する。

行ごとの最大絶対値で方程式をスケーリングし、ILU(0)を構成する。右前処理付き再始動GMRES、2回の修正Gram–Schmidtを使う。Hessenbergの推定残差は検査のきっかけにだけ使い、合格は必ず**元のAとbによる真の相対残差**で判定する。接線を対称化せず、対角シフト・人工安定化・密LUへの自動切替はしない。

ILUのゼロ/小ピボット、Arnoldi破綻、反復上限、入力不正、中断を別の線形エラーとする。ILU(0)はピボット交換を行わないため、正則な行列でも前処理が作れない場合がある。その失敗を物理的な地盤破壊と判定しない。

受理および検証用棄却の試行に、線形反復回数・再始動回数・真残差・CSR格納量・作業配列量・求解秒数を`attempts[].linear_solves`として保存する。線形エラーなどで試行全体が失敗した場合は停止理由を記録するが、それ以前の部分的な線形診断をすべて保存する実装にはなっていない。CSR格納量/作業配列量とOSのプロセスピークメモリーは別の数値。

材料点は従来と同じ不変の確定状態から更新する。材料接線とサブステップの規則、荷重切戻し、掘削面の保持力解放は疎求解への移行で変更していない。境界面は一度作った共有面トポロジーから判定し、各境界面で全要素を走査しない。

## チェックポイント

入力済みのメッシュ・地盤・工程を使い、受理増分数を指定する。次は最初の受理増分で止める例。出力名は未作成のものを使う。

```powershell
$groundRun = 'D:\Work\soil-fem-3d\runs\sparse-restart-test-20261023-214652-380e2009\cad\ground'
$stageFile = 'D:\Work\soil-fem-3d\examples\stages-local-sparse.json'
& .\target\release\soil-fem-3d.exe checkpoint-stages "$groundRun\mesh.json" "$groundRun\ground.json" $stageFile 1 .\new-checkpoint.json
& .\target\release\soil-fem-3d.exe resume-stages "$groundRun\mesh.json" "$groundRun\ground.json" $stageFile .\new-checkpoint.json .\new-resumed-result.json
```

`checkpoint-stages`は、指定数の増分を受理した後に`soil-fem-checkpoint/1`を書き、終了コード0で停止する。試行の途中を保存しない。指定数が全受理増分数を超えた場合は、完了結果をチェックポイントとして扱わずエラーを返す。`solve-stages`は従来どおり最後まで実行する。

保存内容は全節点変位、全材料点の有効応力/e/pc/塑性ひずみ、撤去済み要素を含む履歴、有効要素、外力・仮保持荷重、工程進行度、完了工程結果、試行記録。さらに**現在工程の開始確定状態、終端物理荷重、初期保持力、次の増分、試行数、棄却注入の消化状態**を保持する。掘削途中から再開するとき、撤去済み要素を再び掘削対象として扱わず、元の工程開始状態から解放経路を再構成する。

メッシュ・地盤・工程の型付きJSONを正規化してSHA256を計算し、全保存履歴のチェックサムと、実行中EXEのSHA256/パッケージ版も保存する。JSONの空白等は影響しないが、入力値・ソルバー設定・メッシュ内のCADソースハッシュ・バイナリが変われば再開を拒否する。チェックサムはファイル破損検出であり、第三者による署名や認証ではない。

再開時はチェックサムに加え、初期K0状態、工程開始状態、有効要素・保持力・荷重経路の再構成一致、材料状態、支持/同一平面変位、現在の力残差を確認する。以前に棄却した試行や次の刻みも復元し、試行上限を再開でリセットしない。

CLIのJSONは同じ出力先ディレクトリ内の新規一時ファイルへバッファ付きで書き、内容を同期後、既存パスを置換しないhard linkで原子的に公開する。通常のNTFSで検証した。未公開の一時ファイルは後処理で除去する。既存の結果・チェックポイントは上書きしない。

保存要求で止めた内側のRunは`complete=false`、理由`checkpoint_requested`。完了後の`resume-stages`出力だけを従来のVTUエクスポータへ渡す。突然のプロセス終了で自動的に最新チェックポイントを作る機能はまだない。保存済みチェックポイントからの再開と、SolidCraft画面内の中断/監視は別で、後者はP5の未達要件。

固定静水圧は保存された地盤入力から再構成する。未実装の構造部材・接触・過渡水理の履歴まで保存できるとは扱わない。また、再開時に起動中SolidCraftへ自動照会して変更を失効表示するP5機能も未実装。

## ベンチマークと検証

`benchmark-tet4 CASE NEW_REPORT`は実Tet4の線形剛性を組み立て、既知の離散変位から作った荷重を解く。`verify_sparse.py`は小モデルの行列と解をNumPyの独立組立/密求解で比較し、約10³/10⁴/10⁵要素を測定する。これは全非線形解析の10万要素性能を保証しない。

`verify_checkpoint.py`は、塑性プレロード→除荷→局所掘削の途中で保存し、別CLIプロセスで再開する。確定履歴、工程結果、増分列の完全一致、入力/破損拒否、既存ファイルの保持を検証する。単体試験では各受理増分でJSON往復と再開を繰り返す。

[実測と証跡](SPARSE_RESTART_VERIFICATION.ja.md)／[P6全体の未達要件](P6_STATUS.ja.md)。数値解法の設計参考：[Netlib Templates](https://netlib.org/templates/templates.html)。ハッシュAPI：[RustCrypto sha2 0.10.9](https://docs.rs/sha2/0.10.9/sha2/)。
