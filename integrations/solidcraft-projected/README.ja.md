# 最新の途中差分（2026-10-26 更新3）

現在は55ファイルです。実GUIで使った54ファイル版アプリはビルド済みで、実CADの開始・掘削途中の停止・保存再開・パケットv2の手動読込と3D表示を確認しました。ジョブ完了時の自動読込は旧view/1固定条件で失敗しました。掲載する最新差分には候補提示と完了スキーマ判定の修正を含み、後者の関連試験は採取時点で実行中です。最新55ファイル全体のアプリ再ビルド・全CI・自動完了の受入は未完了です。

固定上流への適用方法・ライセンスは下記を参照してください。下記の旧件数とビルド状況は前回記録です。

# 新P1体積MCCのSolidCraft隔離実装・途中差分

固定した上流コミット `edf57cec293dcea301201f1d5173b33e16351466` へ適用する途中差分です。隣の差分と重ねず、クリーンな別チェックアウトへ適用します。上流のライセンス・著作権通知を保持してください。

```powershell
git clone https://github.com/tangency-labs/solidcraft.git external/solidcraft
git -C external/solidcraft checkout --detach edf57cec293dcea301201f1d5173b33e16351466
git -C external/solidcraft apply --check ../../integrations/solidcraft-projected/soil-fem.patch
git -C external/solidcraft apply ../../integrations/solidcraft-projected/soil-fem.patch
```

クリーンな上流Gitインデックスへ適用し、変更48ファイルのGit正規化後SHA256一致を検査しました。通常のフルクローンで上流手順に従いビルドしてください。疎チェックアウトの場合は上流試験の `skills/` も必要です。

起動時に `SOLIDCRAFT_SOIL_FEM_ROOT` をこのFEMリポジトリ、`SOLIDCRAFT_SOIL_FEM_SOLVER` を使用するビルド済みソルバーへ設定します。[ジョブ手順](../../docs/SOLIDCRAFT_JOBS.ja.md)と[新方式の意味・範囲](../../docs/SOLIDCRAFT_PROJECTED_FIELDS.ja.md)を参照してください。

Tet10変位・材料別P1圧力/体積方式の選択、乾燥MCC入力検査、パケットv2読込、材料/釣合い応力の区別、P1補間、節点/材料の組と生積分点の照会を実装しました。関連engine試験29件が合格し、実データを必要とする既存1件は未実行です。手で値を設定した転送試験は実FEMの精度や実画面の受入を証明しません。アプリ本体ビルド・実SolidCraftモデルによる通し検証・新方式の全体CIは未完了です。従来版の4ファイルのCI修正はこの隔離差分には含みません。完成版ではありません。
