# Issue #36 診断ランナーの Linux プロセスグループ終了

PR #13 の起点は `ad683b30c6cda1f2567a61157d977184731e8fd8`。使い捨ての Python 子・孫で、ランナーの PID だけに SIGINT を送った。修正前は直接の子が終了した後も孫が `S` で生存し、ランナーは `failed` を書いて終了 1 を返した。タイムアウト対照は終了 2、`incomplete`、`timed_out=true` で、孫は `Z` になった。このクラウドの PID 1 は回収を行わないため、再現用の親で回収して消失を確認した。修正前の割込みで残った孫も、再現用の親から停止・回収した。 [Issue #36](https://github.com/tak063495-prog/solidcraft_soil_analysis/issues/36)・[再現手順と実測](REPRODUCE.ja.md)

## 修正

Linux では `start_new_session=True` で直接の子を所有グループのリーダーにする。`waitid(WNOWAIT)` で終了を観測し、回収前に `killpg(SIGKILL)` を実行する。グループ停止までリーダーの PID を保持し、回収後は `killpg(pgid, 0)` による存在確認だけを行う。既に消えたグループの `ProcessLookupError` は許容する。セッション外のプロセスは対象に含まれない。

待機終了、タイムアウト、待機中の例外は同じ後始末を通る。ランナーを実行中だけ Linux の subreaper に設定し、所有グループから引き取った孫をグループ指定で回収する。直接の子と引き取った子孫を回収し、グループの消失を確認してから結果を記録する。起動時と後始末中はランナーで SIGINT を保留し、起動する子には元のシグナルマスクを渡す。

直接の子が先に終了しても同じグループの孫を回収する。`done`／`failed`／`incomplete` の判定は従来の直接の子の終了コードと待機結果に従う。SIGINT や待機例外は `failed` と終了 1、タイムアウトと子の終了 2 は `incomplete` と終了 2。全ケースで `numerical_acceptance=false`。Windows の停止状態での起動、Job への割当・再開、タイムアウト時の Job 終了、Job のクローズは従来どおり。

旧 2 本のソース指紋は保存時点のランナーを指す。指定起点の Git blob と旧 2 本を照合し、現行ランナーは [今回の SHA256](acceptance.json) と個別に照合する。P6-01 の追補台帳は修正コミット `ea2542403b98c334e957299dcec3dcf4fa965167` の `acceptance.json` を 40 桁コミット、Git blob ID、SHA256、型付き状態値で固定する。既存台帳のバイト列は追補後も接頭部に保持する。

## 検証

独立再現の修正後は SIGINT とタイムアウトの両方で、ランナー終了時に直接の子・孫・グループの 3 つが消失した。pytest 9.2.0 の軽量回帰試験は 7 成功、失敗 0、skip 0。親 PID だけへの SIGINT、待機例外、タイムアウト、通常終了、子の終了 7、子の終了 2、直接の子の先行終了を検査した。試験の後始末は、各試験自身が起動し回収権を持つプロセスに限る。

Linux 固有試験は `tests/linux` に置き、pytest の通常探索または明示パスで実行する。Windows 上でこの pytest ファイルを指定すると 7 件が skip になる。Windows CI は従来の unittest 探索を使い、Linux 固有のサブディレクトリを探索対象に含めない。既存の CI ランナーは実行中に発生した skip を検出すると終了 1 を返す。

## Linux 全件の実測

| 対象 | passed | failed・error | skipped | 未実行 |
|---|---:|---:|---:|---:|
| Rust 全ターゲット | 467 | 0 | 0 | 0 |
| Rust doc-test | 0 | 0 | 0 | 0 |
| Python CI と同じ 7 群（tests・tools） | 981 | 0 | 2 | 0 |
| Linux 専用 pytest | 7 | 0 | 0 | 0 |
| Python 合計 | 988 | 0 | 2 | 0 |

Python の 2 skip は PowerShell が必要なゲート試験。Linux では skip として集計し、Windows CI で確認する。Rust の ignored と filtered out はともに 0。P6-01 の 372 出典で整合性を確認した。旧追補台帳の接頭部を保持し、追補後の SHA256 は `95e46cbfd6ac22cd42b9d881cc5f64b39be8b8cc631d62a05e0604cd909beb36`。 [全件の内訳](REPRODUCE.ja.md)

A-01〜A-04、B-01〜B-03、C-01〜C-04、NC-01、P6 は未達、passed=false。
