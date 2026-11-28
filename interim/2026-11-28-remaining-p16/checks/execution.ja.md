# 実行記録

起点 SHA の疎チェックアウトを計算用にも用いた。使い捨てバイナリの初期ビルドは、過去の include と p8 の予測器入力が疎チェックアウトに含まれず停止した。入力ファイルを同じ SHA から復元し、最終ビルドは終了コード0。NC の初期化も以前の入力証拠を復元してから実行した。これらの停止時点で全体試行の追加は0だった。[include](build-sparse-stop.log)・[再試行](build-sparse-retry.log)・[予測器入力](build-missing-predictor.log)・[NC 入力](nc-missing-input.log)・[最終ビルド](build.log)

C-02 の最初の呼び出しは、旧スクリプト内の変位変数とモジュール名が重なり、831の面積分に入る前に停止した。モジュール参照を修正し、正負を開始位置から再実行した。最終方法で再計算した28件の JSON・NPZ は全件バイト一致した。[停止ログ](face-plus-name-stop.log)・[全28件の SHA-256](face-replay.json)

p16の部分テストの最初の実行では p15 の比較用 JSON が疎チェックアウトに含まれず1件が停止した。同じ SHA の比較用 JSON を復元すると4件とも通過した。その後、全結果を対象に9件すべてが通過した。[初回](partial-tests.log)・[復元後](partial-tests-restored.log)・[全件](p16-tests.log)
