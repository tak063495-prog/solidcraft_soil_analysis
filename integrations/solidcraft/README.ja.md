# 従来方式のSolidCraft統合・CI修正の途中差分

固定した上流コミット `edf57cec293dcea301201f1d5173b33e16351466` へ適用する途中差分です。隣の差分と重ねず、クリーンな別チェックアウトへ適用します。上流のライセンス・著作権通知を保持してください。

```powershell
git clone https://github.com/tangency-labs/solidcraft.git external/solidcraft
git -C external/solidcraft checkout --detach edf57cec293dcea301201f1d5173b33e16351466
git -C external/solidcraft apply --check ../../integrations/solidcraft/soil-fem.patch
git -C external/solidcraft apply ../../integrations/solidcraft/soil-fem.patch
```

クリーンな上流Gitインデックスへ適用し、変更51ファイルのGit正規化後SHA256一致を検査しました。通常のフルクローンで上流手順に従いビルドしてください。疎チェックアウトの場合は上流試験の `skills/` も必要です。

起動時に `SOLIDCRAFT_SOIL_FEM_ROOT` をこのFEMリポジトリ、`SOLIDCRAFT_SOIL_FEM_SOLVER` を使用するビルド済みソルバーへ設定します。[ジョブ手順](../../docs/SOLIDCRAFT_JOBS.ja.md)と[新方式の意味・範囲](../../docs/SOLIDCRAFT_PROJECTED_FIELDS.ja.md)を参照してください。

実GUIの開始・協調停止・保存再開・再起動後再開は旧変位型MCCで限定検証済みです。Windows診断パス、CRLF、解析パネル幅の試験を修正し、全試験バイナリの実行を既定最大4並列にしました。全件実行のキュー試験を含むxtask33件は合格。総合CI v5は書式・資産・全ターゲットclippy・WASM・docを通過し、全試験は採取時点で実行中です。全体CIの合格ではありません。新P1方式は隣の隔離実装を参照してください。
