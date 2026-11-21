# 全件検証の件数と初回失敗の保持

[validation.json](validation.json)はソフトウェア検証の確定記録で、数値受入とは別。
A/B/C各項、NC-01、P6は全体未達。Windows CIはpush後も未確定。

|Linuxの最終選択|発見|実行・成功|失敗 / error|未実行|skip / ignore / filter除外|
|---|---:|---:|---:|---:|---:|
|Rust全ターゲット＋doc|467|467|0|0|0|
|Python tests/tools全発見|979|979|0|0|0|

Pythonの新規13件とtools 2件は979件の内数。Rust doc内は0件。
Rustは実行後の独立`--list`の全467 IDと実行467 IDが完全一致した。
Pythonは初期成功分と標準のクラス共有初期化での完走分を選び、最終発見979 IDと完全一致。
ID欠落・最終選択内の重複・expected failure・unexpected successは0。
fmt、locked全ターゲットclippy (`-D warnings`)、locked release全ターゲットbuildは終了0。

## 成功へ換算しない初回記録

初期発見は977件（当時の新規試験11件を含む）。開発終了後の新規試験は13件で、全件を最終選択で実行した。
最初の実行器はクラス初期化を各試験で繰り返していたため、完了546 IDを保持してSIGINTで中断し、
残りを標準のクラス単位へ戻した。中断した1 IDに成功/失敗を捏造せず、後の完走を選んだ。

|全件実行の生の試行（固有ID数ではない）|test開始数|成功|test失敗|test error|fixture error（test未開始等）|skip|
|---|---:|---:|---:|---:|---:|---:|
|最初の部分実行|771|545|1|225|1|0|
|クラス単位での継続・再実行|442|439|0|3|73|0|
|計|1213|984|1|228|74|0|

上表と別に中断1件を保持。最初の977 IDのうち546 IDだけがその実行で最終化され、
残る431 IDは中断1＋未着手430である。継続側は436 IDを割り当て、全436を成功させた。
修正したCIクラス3 IDは初回結果を置き換えた。最初から選ぶ543 ID＋継続436 ID＝979固有ID。
リトライ中の成功重複もあるので984を固有成功件数として報告しない。

欠落archiveのFileNotFoundErrorはログを残して、必要な固定blobだけを取得してから再試験した。
欠落のまま成功/skipにしたものではない。fixture errorはtest開始前のエラー等で、test開始数へ加算しない。
唯一のtest失敗はCI sparseテンプレートの固定行数10に対し、明示例外追加後11となったもの。
11へ更新し、旧profileの全選択契約と新4partのPOSIX/NUL・LF/CRLF契約を実Gitで検査した。

別枠の準備・単独回帰は5試行61件＝59成功/2error/0skip。
親ディレクトリも消えるsparse復元と、固定順残差関数へlistを渡した試験fixtureの誤りを直した。
初回10/11、12/13のログを残し、これらを979へ二重加算しない。
登録後のP6/CI/新試験の結果も[post-registration.json](post-registration.json)へ別記する。

## 封印・ソース・実行環境

科学記録は61,109,258B・4part・663entryを全SHA照合。
全件ログ・全試行/ID・環境・追加因子観測・補助スクリプトは[full-suite-index.json](full-suite-index.json)の
ZIP（12,135,145B・1224entry）へ封印し、全entryを独立に再読込してSHA照合した。
初期ログ、中断、依存blob取得/キャッシュ解放、PowerShell公式配布SHAの照合も収録する。
実行版本体・観測器・入力・丸め行列・Cholesky因子/三角求解途中値のSHAを分けて記録した。

科学段階と検証段階のコード1653件のSHAを別々に固定した。
検証用bindingは継続実行中の保存で、最初の実行より前の固定とは主張しない。
検証binding保存後は全1653コード/設定SHAが終了時にも一致。
最初の実行開始後の変更はCI行数試験とP6追記/監査器pinだけで、影響範囲は再試験した。
追加の因子解析はwork内の明示観測版のみで、生産Rustには変更0。

全件に必要な旧NCの大きい封印もクラウドで選択取得し、過去保存場を再計測した。
新規NC求解、PC計算、旧記録の書換え、29GBのフル作業ツリー復元ではない。
PowerShellは公式7.5.4配布のSHAを照合して作業用領域へ導入し、一般のXDG格納先だけを作業領域へ指定した。
HOME・NCコード/設定は変更せず、PowerShell試験をskipにしなかった。

P6-01の科学出典は固定commit `18d54b140bf44b6ca20da2745c3b4dd367437277`。
確定validationも別の固定commitへ保存して登録し、旧台帳の全バイト接頭部と監査器pin以外の同一性を確認する。
最終登録のcommit/各SHAと監査結果はscientific-registration.json、validation-registration.json、
post-registration.json、final-p601-audit.json、preservation.jsonを参照する。
