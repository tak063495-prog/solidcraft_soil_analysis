# SolidCraft統合の途中差分

上流SolidCraftの固定コミット `edf57cec293dcea301201f1d5173b33e16351466` に適用する差分です。解析属性、ジョブ操作、力学結果表示を含みます。上流のライセンス・著作権通知を同梱しています。

```powershell
git clone https://github.com/tangency-labs/solidcraft.git external/solidcraft-soil-fem
git -C external/solidcraft-soil-fem checkout --detach edf57cec293dcea301201f1d5173b33e16351466
git -C external/solidcraft-soil-fem apply --check ../../integrations/solidcraft/soil-fem.patch
git -C external/solidcraft-soil-fem apply ../../integrations/solidcraft/soil-fem.patch
```

差分は固定したクリーンな上流ソースへ実際に適用し、変更ファイルの内容一致を確認しました。上流のビルド手順に従ってアプリをビルドし、起動時に `SOLIDCRAFT_SOIL_FEM_ROOT` をこのFEMリポジトリへ設定します。詳細は `docs/SOLIDCRAFT_JOBS.ja.md` を参照してください。

P1体積の新しいMCC形式は外部CLIの限定検証段階です。native画面の新形式対応と協調停止の接続は継続中です。
