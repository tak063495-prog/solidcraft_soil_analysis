# 不合格・拒否の保持

- outside-offset-scope.log / failed-studies.zip の offset-positive：t=0.1 mにe=0.1 mを与えた範囲外入力の明示拒否。受理範囲を拡張せず、最終のt/2入力と区別する。
- different-accepted-path.log / failed-studies.zip の offset-positive-2：棄却なし対照で除荷・撤去の刻み列も変わってしまい、全フレーム厳密一致に失敗。同じ受理列に固定した後の結果へ読み替えない。
- verification-flag-comparison.log：全物理場・受理Attemptは一致したが、全Run比較に意図的棄却済みフラグを残して失敗。最終では棄却操作の入力/フラグ/そのハッシュだけを区別する。物理量・受理履歴の許容値を緩めない。
- development-compile-error.log：開発中のテストコードが非Optionのlast_committedへas_mutを呼び、コンパイル失敗。修正後の全件実行と区別する。
- scope-refusals（evidence.zip内）：曲面二次辺・面外中間節点・円筒/球面形状、範囲外e、未検証組合せの拒否入力と終了コード/stderr。古典殻ベンチマークの合格へ転用しない。

Python全件の初回には疎チェックアウトの過去証拠ファイル不足も含む。元のGit blobと照合して配置した後の完全再実行を最終判定に使い、初回ログを保持する。

`rust-missing-frozen-fixtures.log` は、過去のGriffiths不合格JSONが未配置だった全Rust初回の実エラーを保存する。該当5ファイルは元Git blobと照合して復元し、内容・判定は変更していない。Pythonの未完了診断3実行は、各JSONでpassed=false/completed=falseを明示する。最終の無フィルター全件実行と区別する。
