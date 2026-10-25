# 新規環境CIと追加診断の途中成果（更新7）

FEM専用のGitHub Actionsを追加しました。固定したActionsとRust、Python、依存ファイルを使い、Windowsでfmt・clippy・release試験・CLI build・Python全試験を実行します。初回のリモート成功はこのスナップショットでは未確認です。SolidCraft全体CIやP6全体の合格を示しません。

公開前のローカルPython全67件は合格。Rust83ファイルとCargo定義は全247件合格版と一致しています。SolidCraft全体CIのseed 23失敗を同じバイナリで単独再現しました。上流engineによる別ビルド比較は進行中です。基準値は変更していません。

MCCの新しいLUごとの残差補正を含む別順序の参照でも、元の厳密な成分別比較は未達でした。過去の失敗と全非線形経路の未受入を保持しています。

今回のsource-sha256.jsonに挙げた全ファイルは、公開先の同名ファイルとバイト単位で一致します。過去の更新6では公開用READMEの装飾をローカル指紋採取後に加えていたため、そのREADMEだけは当時のソース指紋と異なります。実行ソースの差ではありません。

[CIの範囲](../../docs/CI.ja.md)、[全残件](../../docs/P6_STATUS.ja.md)、[直前の途中成果](../2026-10-26-plate-reference/README.ja.md)。検証要約中のローカルrunパスは由来を示し、全runデータの同梱を意味しません。
