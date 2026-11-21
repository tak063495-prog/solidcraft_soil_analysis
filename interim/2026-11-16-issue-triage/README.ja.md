# 全issue照合と残件整理：範囲拡張(d)

依頼の保存先ラベルは2026-11-16、実施日は2026-11-11（日本時間）。PR #13、`codex/plastic-spatial-convergence` の指定先端 `d8616cefb526e70f422ec8c164afd2353c8f8149` から作業した。検証用Pythonの最小修正と回帰は `bfec86f0ff119f44a8046853258e3daf10a992f7`。PRはマージしない。

**起点では#15は解決済み。closedの#5/#14には追記の未解決経路があり、今回修正した。修正後は全9件の対象不具合を解決済みと判定する。** 元の閉鎖根拠の回帰と、以前から残る別経路を混同しない。NC・湿潤JOINT・曲面板・P6の数値受入を完了にしたものではない。

## 全issueの一覧と判定

GitHub検索 `is:issue`、リポジトリを限定、上限100で全9件。PRをissue件数へ含めない。各本文と全コメントも取得し、検索では省略されるstateを各issueの取得結果で照合した。取得時はopen #15、closed #1〜#5/#10/#11/#14。原文は[取得snapshot](github/issue-snapshot.json)、機械可読判定・全試験ID・修正SHA・行番号は[issues.json](issues.json)。

|issue|取得時|d8616ceでの判定|今回修正後|対象試験（全件内）|根拠と範囲|
|---|---|---|---|---|---|
|[#1](https://github.com/tak063495-prog/solidcraft_soil_analysis/issues/1)|closed|解決済み|解決済み|Python test_p2_resume_ownership 9件|所有Popenの終了確認後だけunlockし、確認不能時は所有情報付きlockを保持。通常・起動中signal・待機/status例外と安全な再起動を検証。|
|[#2](https://github.com/tak063495-prog/solidcraft_soil_analysis/issues/2)|closed|解決済み|解決済み|Python test_p2_optimized_validation 2件 / Python test_p2_acoustic 1件|構造照査・入力ZIP・exact診断の外部整合条件は最適化で消えない明示例外。normal/-O/PYTHONOPTIMIZE=1で不一致と正常対照を実行。|
|[#3](https://github.com/tak063495-prog/solidcraft_soil_analysis/issues/3)|closed|解決済み|解決済み|Python test_p2_acoustic 4件|各frameの確定element_material_idsをsampled/exact共通で解決。旧/新Mの候補、定数変更、legacy、非活動を含む不正inventoryを検証。|
|[#4](https://github.com/tak063495-prog/solidcraft_soil_analysis/issues/4)|closed|解決済み|解決済み|Python test_p2_epoch_budget 2件|実際に保存するindent=2の同じバイト列で個別/累積予算を検査。geometryとmanifestも含め、超過時は成功マーカーなし。|
|[#5](https://github.com/tak063495-prog/solidcraft_soil_analysis/issues/5)|closed|未解決（元不具合は解消、追記経路が残存）|解決済み|Python test_p2_safe_pause 4件 / Python test_issue_triage_followups 2件|元の親PID再利用は親子の本人情報と作成時刻のedge検査で解消済み。追記の最適化時handle照合消失は起点に残存し、今回は一切の処理より前に最適化モードを明示拒否。|
|[#10](https://github.com/tak063495-prog/solidcraft_soil_analysis/issues/10)|closed|解決済み|解決済み|Python test_p2_richardson 2件|単調Richardsonを適用する前に連続差分の符号と有限性を確認。正負の振動系列は実CLIでIndeterminate、単調既知次数の対照は維持。|
|[#11](https://github.com/tak063495-prog/solidcraft_soil_analysis/issues/11)|closed|解決済み|解決済み|Rust epoch_srm_gate 2件|全完了recordにも対象ground/caseの未低減工程prefixをverify_originで共通照合。正しく再封印したA/B混在を全phaseと保存検証/復旧で拒否。|
|[#14](https://github.com/tak063495-prog/solidcraft_soil_analysis/issues/14)|closed|未解決（元不具合は解消、追記経路が残存）|解決済み|Python test_p2_soil_source 2件 / Python test_issue_triage_followups 3件|元のGP非有限値拒否は有効。追記の後続no_flow流量NaNは起点に残存。さらに負/負無限大の相対水収支誤差を実API/CLIで誤受理したため、全境界流量と縮約値の有限scalar検査、相対誤差の有限・非負条件を追加。既存上限は不変。|
|[#15](https://github.com/tak063495-prog/solidcraft_soil_analysis/issues/15)|open|解決済み|解決済み|Python test_nc_mean_stress_evidence 4件 / Python test_nc_spatial_evidence 4件 / Python test_historical_utf8_isolation 2件 / Python test_spatial_replay_isolation 3件 / Python test_nc_replay_contract 2件|旧archive全part/entryのSHA・size・inventoryと可逆JSONLの復元検査を保持し、検証済み旧snapshotを-I -X utf8の別Pythonで再計測。現行checkoutのsource差は診断欄へ分離し、旧実行を新sourceへ結び替えない。|

## 今回の最小修正と修正前後

修正した既存コードは2ファイルだけ。`src/`、物理モデル・接触則・既定求解、許容値・参照値は変更しない。

- #5：停止候補の最適化実行を入口で明示拒否する。通常モードのsnapshot親子検査と保持handle照合は保持。新規`SafePauseEntry`2メソッドは通常のread-only helpと、`-O`/`PYTHONOPTIMIZE=1`の入口拒否を実別Pythonで検証する。子のCIM/native API/書込はstubで禁止し、実processを終了しない。
- #14：全境界面の流量を集計前に有限scalar検査し、縮約した合計/最大値も検査する。相対水収支誤差は有限・非負scalarを必須にする。元の`1e-9`、no_flowの`1e-12`等の上限は不変。新規`SoilBoundaryFlow`3メソッドは実外面12面、先頭/後続、NaN/±Infinity/不正形状/有限不一致/負の誤差をAPIと全3モードの実CLIで検証する。
- #15：既存の実証拠8試験に加え、新規`HistoricalReplayContract`2メソッドで小型の封印archiveを両replayへ渡す。異なる現行sourceを実行せず、CRLF/LFを別記する正常対照と、part・source SHA・過剰/欠落inventory・source欠落・離散値・strict float不一致の拒否を確認。計16 subTest条件であり、16試験と水増ししない。replay本体・例外・float条件は変更しない。

[修正前後の実API/CLIと偽handle検証](baseline-followups.json)、[修正前の5メソッド・14 subTest失敗記録](baseline-followup-tests.json)、`baseline-followup-tests.log`を保存した。失敗記録を成功件数へ加算しない。修正前の後続no_flow NaN、相対誤差-Infinity/-1は全3モードでexit0/passed=true。修正後は全9対照でexit1・reportなし。正常reportのcanonical SHAは修正前後とも `f337d08cf689ef038c8f6f910f5e688634b1866db69a41edf43080ccb6df1cdb`。この正常fixtureの一致と、既存の大規模数値結果の再求解は区別する。

## #15：歴史的証拠と現行HEADの分離

meanの隔離再計測導入は`113d83b`、明示UTF-8は`b235824`、spatialの隔離と不足依存の旧SHA限定補完は`c0addfb`。成功した参考先端8d62e7b/cecbb39/bc6eb1cはこれらを継承する。

|歴史的証拠|元archive SHA256|サイズ/part/entry|封印source数|
|---|---|---|---:|
|mean|`87434944b291e5958f8a98e3d36ff485a64adbe95be67f42d3fbbe2cee25f6f9`|5,490,818,201 bytes / 328 / 1,314|783|
|spatial|`cc268fb3b4b0f45a2774ef13125ff0a05431629e3de45acffca22a05ca45dcd4`|920,331,265 bytes / 55 / 411|32|

meanの783 source指紋はarchiveのsource entryそのもの（別のsource-fingerprints.jsonはない）。spatialには元のsource-fingerprints.jsonが623件あり、不足依存の補完もこの旧SHAで制限する。part、entryのsize/SHA/inventory、meanの可逆JSONL復元SHAを照合してから、検証済み旧snapshotを`-I -X utf8`の別Python・旧sourceをcwdとして実行する。新しいcheckerへの差は診断欄に記録する。spatialの不足依存は元のfingerprintと完全に同じbytesだけを補完する。現行sourceを旧sourceと同一だと要求し続けることと、旧証拠の真正性を検査することを分離した実装である。

両archiveは実行バイナリを配布せず、元の実行版SHA/実プロセス・要求応答SHAを記録する仕様。今回はその記録とsource snapshotを不変に保ち、旧バイナリの再ビルド一致・旧Rust求解の再実行を認定しない。現行epoch sourceのSHAが旧manifestと異なるのは予期する差であり、旧実行のhashを新sourceへ更新していない。[全source事前照合](source-binding-audit.json)と[保存照合](preservation.json)を参照。

## Windows CIを試験ログまで確認

[ci-summary.json](ci-summary.json)と`github/job-*.log`を保存した。ログはconnectorが返したdecoded textをLFへ正規化したもので、artifact ZIPのバイト取得ではない。

|head|run|mean実行|spatial実行|判定|
|---|---|---:|---:|---|
|8d62e7b|[42890061150](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/42890061150)|4件ok|4件ok|全job success|
|cecbb39|[42925852344](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/42925852344)|4件ok|4件ok|全job success|
|bc6eb1c|[42938912424](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/42938912424)|4件ok|4件ok|全job success|
|d8616ce|[42951479845](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/42951479845)|4件ok、job130296554677|4件ok、job130296554588|全job success、Rust453/Python505|

d8616ceのCI merge `5510d3b7ebc940e324b675a544c10dea90b02b2c` とheadのtreeはともに`9354363a57139ab28334f1e587b5a479ef5cfa52`。各NC4試験の名前・`... ok`・`Ran 4 tests`・`OK`がログにあり、クラスsetup停止や試験skipではない。Python505 IDの開始と集計を照合した。行の途中に子出力を挟む20試験と、同一行でokを出す485試験を区別し、後者だけを全件数にしない。最終`OK`とrunnerのskip禁止条件も確認した。matrixで非選択のstepがskippedなのは、選択された試験自体のskipとは別である。

WindowsのNCは離散証拠一致、受入falseを維持し、計測JSONはバイト一致ではない。d8616ceではmean157箇所、spatial18箇所の派生float差を記録する。旧成功、キャンセルや時間切れ、CIが緑であることだけを数値収束・現行FEMの認証へ転用しない。

## Linux全件検証

fmt、clippy（`--locked --all-targets -- -D warnings`）、locked release build、**Rust451件・Python512件（tests 510 + tools 2）**が成功。fail/error/ignored/skipped/filtered-out/expected-failure/unexpected-successは0。Windows固有cfgの2件はLinuxではコンパイルされず、d8616ce Windowsの453件に含まれる。ignored/skipで除外したものではない。新規7メソッドは全件に含み、集中再実行は加算しない。[validation.json](validation.json)、[Python発見ID/実行4群の照合](python-coverage.json)、`validation/`が根拠。

Rustのチェックは取得待ちの間にd8616ceのGit blobから作った隔離build checkoutで実行した。全144 Rust/config入力のGit blobとSHAを検査し、最終ブランチとも同一bytesであることを確認した（`rust-build-source-binding.json`）。Python全4群は指定ブランチのcheckout上で実行し、mean/spatialの実archiveも再計測した。ビルド場所を変えた現行CLIを旧実行版と同一バイナリとは扱わない。初回は隔離コピーに過去の斜面不合格JSONを配置し忘れ、lib 59成功/1失敗で後続が未実行だった。`retained-failures/initial-rust-*`へ保存し、元Git bytesの5記録を配置して全件を再実行した。初回59件を最終件数へ加算しない。

## 残件・issue操作・再現

[todo.ja.md](todo.ja.md)に(a)(b)(c)とP6/残件表の未達・拒否を集約した。既存2文書へ本整理への追記だけを行い、P3b/P4/P6の状態行と過去の本文は保持する。低プレロードNC全撤去、湿潤JOINT塑性すべり全履歴収束、曲面Tri6、大規模/実CAD/GUI/GPUは完了にしない。既存のFs両区間、Griffiths、柱Tet4、杭、降雨の不合格/判定不能は変更しない。

コメント文案は[comments/](comments/)。既存closedは維持し、#5/#14には追記経路の修正根拠、#15には解決根拠を投稿しcloseする方針。取得時のstateと実操作後のstateは区別し、push後の実操作結果はGitHub issue timelineと最終報告を参照する。投稿/更新権限が拒否された場合は文案を残す。今回push後のWindows CI成功は本書で確定させない。PR #13のマージ操作は行わない。

[REPRODUCE.ja.md](REPRODUCE.ja.md)に全件実行、小型反例、旧source隔離、ユーザーPC専用の大規模・実CAD・実停止・GUI/GPU手順を記載する。
