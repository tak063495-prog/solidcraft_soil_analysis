# SolidCraft統合の途中差分

上流の固定コミット `edf57cec293dcea301201f1d5173b33e16351466` に適用する差分です。解析属性、ジョブの開始・協調停止・再開、結果表示を含みます。上流のライセンスと著作権通知を同梱しています。

```powershell
git clone https://github.com/tangency-labs/solidcraft.git external/solidcraft-soil-fem
git -C external/solidcraft-soil-fem checkout --detach edf57cec293dcea301201f1d5173b33e16351466
git -C external/solidcraft-soil-fem apply --check ../../integrations/solidcraft/soil-fem.patch
git -C external/solidcraft-soil-fem apply ../../integrations/solidcraft/soil-fem.patch
```

固定したクリーンな上流のGitインデックスへ実際に適用し、変更47ファイルの内容一致を確認しました。通常のフルクローンでビルドしてください。疎チェックアウトでは上流MCP試験の参照に `skills/` も必要です。

上流手順でアプリをビルドし、起動時に `SOLIDCRAFT_SOIL_FEM_ROOT` をこのFEMリポジトリ、`SOLIDCRAFT_SOIL_FEM_SOLVER` を使用するビルド済みソルバーへ設定します。操作・信頼する実行ファイル・適用範囲は [ジョブ手順](../../docs/SOLIDCRAFT_JOBS.ja.md) を参照してください。

実GUIの開始・停止・再開・アプリ再起動後の再開を限定検証しました。P1体積の新しいMCC形式は外部CLIの検証段階で、native画面は未対応です。全体CIは採取時点で実行中、P5/P6全体は未完です。
