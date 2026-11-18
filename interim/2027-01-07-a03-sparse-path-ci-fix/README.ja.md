# A-03 疎選択試験：Windowsパス引用の修正

PR #13、codex/plastic-spatial-convergence、起点6ee03451c49b86bc9f238ba0bed2f21b15f4299cへの追記。マージしない。
2027-01-07は保存先ラベル、作業日は2026-11-18 JST。
**A-03は未達、passed=false。数値・許容値・参照値・合否・拒否域・生産Rust・NCを変更しない。**
科学診断としての新しいFEM求解は0（既存全件試験内の求解は別）。12 Tet10の旧科学記録・生証拠・封印・失敗・数値を全て保持する。

## 失敗と原因

[Windows run 44234435017 / job 134567021528](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/44234435017/job/134567021528)
のpython-core-1は492件実行・491成功・1失敗、error/skip 0。
失敗は`test_ci_sparse_selects_new_archive`の91行目。
要求されたブランチ先端は6ee0345、ログが示す実チェックアウトはPRの合成merge ref
93a443ebfc3d9b34e73005b9c9f3e91012fb211cである。合成refの実行をブランチSHAそのものの実行とは表記しない。
[失敗ログ](retained-failure/windows-job-134567021528.log)はconnectorで取得した復号済みテキストを
UTF-8、BOM/CRLF保持で保存した。サービスの元ZIPとは区別し、SHAを[出典情報](retained-failure/windows-ci.json)へ固定した。
元失敗は成功へ書き換えない。

`str(Path.relative_to(ROOT))`はWindowsでバックスラッシュ区切りになる。
旧試験はこの値をGitの行区切りstdinへ渡し、stdoutを`.splitlines()`だけで比較していた。
Gitはバックスラッシュ等を含む行出力をC形式で引用し、二重引用符とエスケープを付ける。
Pythonの失敗表示はさらにreprでバックスラッシュを表示するため、ログには二重化が重なって見える。

親CIステップはPowerShell 7だが、問題の子Gitは`subprocess.run`の引数配列・既定の`shell=False`で直接起動している。
PowerShell/cmdによる再引用や改行だけが原因ではない。
LinuxでもPureWindowsPathで同じ入力を作り、直接Gitを呼ぶと同じ引用不一致を再現した。
さらに[実際の試験本体の対照](actual-case-probe/result.json)では、Path.relative_toの返値だけをWindows形式にして実Gitを呼んだ。
起点の旧試験は期待どおり1失敗、修正版は両core規則で成功した。旧失敗トレースも保持する。
これはLinux上のWindowsパス対照であり、Windows CI実行そのものの成功とは区別する。
LF/CRLF入力、元stdin/stdoutの生bytes、コマンドと結果を[reproduction](reproduction/)へ保持した。
LinuxへCRLFを直接渡す対照をWindows CRTのtext変換そのものとは扱わない。

## 修正と試験の強さ

対象試験だけで、リポジトリ相対パスを`.as_posix()`によりGitの`/`区切りへ統一する。
`git sparse-checkout check-rules -z`を使い、UTF-8＋NUL区切りのbytesをstdin/stdoutで直接比較する。
`text=True`、`.splitlines()`、引用符除去、後付けのエスケープ解除は使わない。
CRLFへのtext変換や行の分割も比較経路へ入れない。

全10アーカイブpartについて、実際のpython-core-0/1疎選択規則を使う検査を維持する。
比較は順序・件数・重複・末尾NULまでbytes完全一致であり、集合化や包含だけの検査へ弱めない。
追加の軽い2試験でPureWindowsPath→POSIXと、空白・引用符・バックスラッシュ・タブ・改行・非ASCII・重複を確認する。
後者はcore.quotePath=true/falseの双方で実Gitを呼び、除外対象が出力されないことも完全一致で検査する。
特殊名は論理的なGitパスであり、Windowsで作成できないファイル名を実ファイルとして作らない。
新しいRust統合試験や重いFEM再計算は追加しない。単独12件は成功、失敗/error/skip 0。

## P6-01と保全

e541fe1 / 2026-12-10の方式を継承する。旧台帳は厳密なバイト接頭部として保持し、
新フォルダのacceptanceと確定validationをそれぞれ先行する不変コミットへ保存してから登録する。
各JSONの作業ツリー一致、Git blob SHA-1、SHA-256、型付きpassed/a03_complete/diagnostic_onlyを検査する。
監査器はSUPPLEMENTAL_SHA256定数だけ更新し、SOURCE_NOT_INVENTORIEDの検査を弱めない。
旧A-03フォルダ全ファイルのSHAは[prior-a03-binding.json](prior-a03-binding.json)に保持する。
旧validationのLinux成功記録と「当該push後Windows CI未確定」の歴史も書き換えず、今回のWindows失敗を別記する。

全件検証の実数は確定後のvalidation.jsonに保存し、このREADMEへ追記する。
push後の新しいWindows CIは未確定。Linuxの成功・エミュレーションをWindows実行成功へ転用しない。


## Linux全件検証の確定記録

fmt、locked/all-targets clippy、release全ターゲットbuild、doctestは終了0。
**Rust467件、Python957件が成功**。failure/error・未実行・skip・ignore・filtered-out・expected-failure・unexpected-successは0。
全発見/割当/実行IDが一致し、欠落/重複0。全件中の1630コード/設定ソースは前後でSHA一致。
単独12件、実試験本体の旧版失敗/新版成功対照、登録後の再確認はこの固有件数に加算しない。
以前のLinux955件成功・Windows491成功/1失敗の記録はそのまま保持した。

|Python群|成功|実測秒（起動/discovery込み）|
|---|---:|---:|
|python-core-0|439|436.769|
|python-core-1|494|361.837|
|python-mean-evidence|4|862.848|
|python-spatial-evidence|4|130.422|
|python-c01-switch|4|47.948|
|python-c01-path|4|77.372|
|python-c01-cause|8|65.597|

全件ログ・ID一覧・process JSON・前後ソースSHAをfull-suite.zipへ封印し、full-suite-index.jsonでSHA固定した。
修正と失敗記録の固定コミットは428e6d82b1f5af925ada61945c919dd15cee0981。
確定validationを759d8b39a659e1fe379fe6ed87ec282b6fbfc900へ固定してP6へ登録した。
最終監査は362出典、integrity_ok=true、SOURCE_NOT_INVENTORIED=0。旧台帳は厳密なバイト接頭部を維持。
登録後のP6関連76件＋A-03の12件、計88件も成功し、全件957件へ重複加算していない。
起点の配置対象13004ファイルのうち13001件はバイト/Gitモード一致。
残る3件は対象試験（固定428e6d8のパッチと一致）、台帳の追記、監査器のSHA定数だけ。
旧A-03フォルダの全ファイルSHAも一致し、旧数値・封印・失敗・NCを保持した。
全件証拠と保全証拠のZIP/member SHAも再読込で確認した。
push後の新しいWindows CIは未確定。全件ソフトウェア成功をA-03の全体受入へ転用しない。

A-03は未達、passed=false。
