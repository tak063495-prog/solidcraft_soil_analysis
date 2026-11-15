# filtered取得・実バイト・キャッシュの記録

固定先端23c7b0eをblob:none/depth=1で取得し、初回全履歴転送は行わない。
初回一括checkoutは中断した。config.worktreeのsparse=falseがlocal設定を上書きしており、
no-lazy取得の失敗も保持した。worktreeの設定を訂正して本文から先に復元した。
必要な固定OIDは本文32個/アーカイブ1個でfiltered取得した。
10,330必要オブジェクトは取得時点で全て存在し、初期不足7,852の回復後の不足0。
旧P6固定コミット49件とP6-03文書blob29件も必要分だけ取得した。
全実ファイルのGit blob/実行属性一致をbase-preservation-before-registration.jsonへ保存した。
既存CIと同じpostfold74パスだけはGit参照を保持して実体化対象外。必要証拠をskipしない。

63GBディスクで二重保持を避けるため、実バイトとSHA/Git blobを再照合した後で
再取得可能な単一blob promisor packだけを削減した。旧作業ファイル/コミットを削除しない。
cache-pruning-*.jsonはpack SHA、元OID、検証済み実ファイルを列挙する。
最終的なGitキャッシュの欠落と必要実ファイルの欠落を区別する。
後のcommitは検証済みindexからwrite-tree --missing-ok / commit-tree / 条件付きupdate-refで単親作成する。
検証済みの旧tree参照を保持し、改変や全履歴の捏造ではない。
生の各取得ログは封印のb02follow-preparation-sealed/fetches-v3/にも保持する。
HOME/NC設定は変更せず、portable PowerShellのXDGだけworkに置く。
