# 取得・実バイト・キャッシュの監査

初回は `fetch --filter=blob:none --no-tags` で指定ブランチを取得した。
一括checkoutの2回の未index packを中断し、A-04のrecover.pyと同じ方式で完全prefixだけを回収した。
455/251オブジェクトをstrict indexし、原header・切断位置・元pack/未完suffix SHAをrecovery.jsonへ保存する。
未完suffixの生バイトはwork/b01follow-git-recovery/{,second/}incomplete-suffix.binにも保持する。
このフォルダにはheader、SHAと再現スクリプトを保存し、未完suffix自体を新しい数値証拠とはしない。

その後は必要な固定OIDだけを32個の本文群/1個の大blobに分け、filtered取得した。
fetch-result-v3.jsonの全10,205必要オブジェクトは取得時点で存在し、残件0。
既存CIと同じpostfoldの74パスだけは実体化しない。新旧の必要アーカイブをskipしない。

実体化後、全12,403起点パス中12,327実ファイル・29,406,169,025 bytesはGit blob/実行属性と完全一致。
残る2実ファイルは旧台帳の元バイトを接頭部に保つ追記と監査器pinだけ。上記74パスのGit参照は保持する。
base-preservation.jsonの全SHA/パスを参照。旧の数値・生証拠・失敗・封印・拒否・NC・生産Rustの変更0。

63GBの一時ディスクでは実ファイルとGitキャッシュの二重保持が問題になる。
実ファイルのGit blobとSHAを再照合してから、再取得可能な単一アーカイブpromisor packだけを削減した。
作業ファイル/コミット/参照値は削除しない。削減したpackのSHA/OID/元実ファイルはcache-pruning*.jsonへ記録。
一部の削減済みpackを進行中Gitが開いていたため、終了後に約31GBの空きを確保した。

通常commitが無関係な旧アーカイブblobをlazy fetchしようとしたため中断し、再取得を止めると
`could not fetch b11d01d...` で不成功になった。旧実ファイルのSHAは先に全て確認済みである。
このpartial cloneでは `git write-tree --missing-ok`、`commit-tree -p`、旧SHAを条件とするupdate-refで
同じ検証済みindexから単親コミットを作る。履歴・ファイルを捏造せず、欠けたキャッシュを数値未完了やskipへ変換しない。
プランの行列は初回plan.jsonと同じround-trip定数、分類器の訂正は付記だけであり、数値を変更していない。
新しい重いデータをpushする際も、旧リモート起点を境界とし、必要な新規オブジェクトだけを送る。

初回commit後の自動repackも中断して元オブジェクトを保持し、この作業cloneだけgc.auto=0/maintenance.auto=falseにした。
HOMEやNC固有設定を変更せず、提供activate.shのUV_CACHE_DIRとwork配下XDG格納先を使う。

初回全件のP6-03は過去コミットの文書blob3個の不足で1 errorとなった。
通常fetchはbad revisionで不成功、v3と同じfetch.negotiationAlgorithm=noopを指定したfiltered取得で復旧した。
p603-fixed-blobs-recovery.jsonに固定Git blobとSHA256の全照合を保存する。
P6-03単独19件の最初の再試行も1 error、取得完了後は19成功、さらにcore-0全426件成功。
固定出典・許容値・試験を変更せず、失敗ログはvalidationに残す。
