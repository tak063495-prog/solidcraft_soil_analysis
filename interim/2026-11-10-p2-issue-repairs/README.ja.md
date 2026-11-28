# PR #13: P2 Issue 8件の修正と回帰記録

起点は `f19ce642647207a87b9dda6cc31b006b85b5936c`、ブランチは `codex/plastic-spatial-convergence`。mainから修正を始めず、この起点を直接の親とするコミットにまとめる。PR #13はマージしない。日付はユーザーのAsia/Tokyo基準。

8件とも起点に不具合が残っていた。本文・全コメントを取得して確認し、修正前の実行失敗を `baseline-*.log/json` に保存した。「すでに修正済み」として省いたIssueはない。ここでいう修正は検証・表示・プロセス管理の不具合の修正であり、一般FEM・P6の受入ではない。

| Issue | 起点での確認 | 修正と回帰範囲 |
|---|---|---|
| #14 | 正常合成fixtureは受理。1 GP/全162 GPのhead/pressure/flux NaNをAPIが誤受理。起点ファイルを別名で再実行した実CLIも6条件すべてexit 0・passed=true（`baseline-issue14-cli.json`）。 | 全GPの座標(3)、水頭scalar、圧力水頭scalar、流速(3)を誤差集計より前に形状・有限性検査。誤差最大値も有限検査。既存許容値を維持し、従来のassert条件も最適化で消えない明示例外へ変更。正常API/CLI、4フィールド×1/全GP×NaN/±Infinity/有限不一致、形状不正を照査。拒否CLIは非ゼロ・reportなし。 |
| #11 | 本文のA=荷重1、B=0.5、要求1.6/最大2でBのFinish(None) payloadをAの外側identityへ接続すると受理。正しい型のpayloadバイト列からsealを計算した2回帰が起点で失敗（`baseline-issue11-sealed.*`）。 | `verify_records(mesh,ground,case,spec)` から各完了originに `verify_origin` を実施。本文A/Bおよびgroundのみ変更をverify/advance/保存ジョブ検証/復旧で拒否。複数gateでConstruct/Search/Finish(None/Some)/Doneを確認。active Searchの起点・保存constructionは対象Aの正しいものを保持し、完了recordだけをBへ差し替える。Doneのreport整合・ディスク外側hashも再計算。既存の中断/保存/復旧と連続結果の一致も試験。 |
| #10 | 実 `plane-strain-2d judge` が振動系列をpassにした（`baseline-issue10.*`）。 | 連続差分の符号と有限性を確認してから既存の正倍率Richardsonを適用。正/負の振動系列はIndeterminate、既知次数1/2/3の単調系列は正しい外挿。差分ゼロ、solver floor、発散、非有限、2点を実CLIで照査。既存5%とsolver floor条件は変更しない。 |
| #5 | 起点の実while-loopをAST抽出し、親23時/子22時で所有集合が{100,200}へ誤拡張する純データ再現（`baseline-issue5.json`）。 | 親・子それぞれのsnapshot identity（PID一致・image・timezone付きCreated）を検査し、各上下edgeに親Created<=子Createdを必須化。不明・欠落・parse不能を所有集合へ入れない。正しい子孫・excludedも試験。既存の保持handleによる本人image/生成時刻照合は保持。Windows CIMや実終了は実行せず、無関係processへsignalを送っていない。 |
| #4 | 実exporterの小型transport fixtureで、実ファイル/累積上限の1バイト超過を成功扱い（`baseline-issue4.*`）。 | indent=2の実保存バイトを一度生成し、同じバイトを検査・保存。geometry/frame/points/manifestを累積。manifest予算確認後にのみcompleteと成功receiptを公開。配列主体の各種payload、個別/累積の上限ちょうど・1 byte超、失敗時manifest/receiptなしを試験。力学受理部分のみ既存fixtureでmock、実FEM大規模表示とは区別。 |
| #3 | 旧M=1.6/新M=1.4、mesh IDs=[1,1]/frame IDs=[2,1]で本来1候補の最終frameが0候補（`baseline-issue3.*`）。 | sampled/exact共通の `element_materials` で各frameの確定IDを解決。空/欠落listのみlegacy fallback、None/型不正/長さ不一致/未知IDは明示拒否。非活動要素のIDも検査。MATSET前後、lambda/kappa/nu変更、同一ID、legacy、sampled→exact実CLIを照査。 |
| #2 | 構造moment/axes/force/history不一致、ZIPのhash/size/inventory不一致が-O/envで受理された。exactの変更入力も最適化で受理（`baseline-issues2-3.*`, `baseline-issue3.*`）。 | beamの照査を明示例外へ、数値比較は形状と有限性も確認。指定interim unpackのmanifest/ZIP検査とexactのhash/実装hash/shape/正値/整合検査も明示例外へ変更。通常/-O/PYTHONOPTIMIZE=1を別processで対照。構造不一致はmanifestなし、ZIP不一致は展開先/成功表示なし、exact不一致はoutputなし。指定2パスは起点に実在。 |
| #1 | 軽量の使い捨てstudyで親だけへのSIGINT後に子が残り、lockが消失。待機/status例外も後始末なし（`baseline-issues1-4.*`）。 | 所有Popenを保持。起動中の協調signalをハンドル代入まで遅延し、割込み・例外時にterminate/必要ならkillとwaitで当該子の終了を確認してからunlock。確認不能なら親/子PID・中断・手動確認事項をlockに残す。progress書込失敗とは独立して記録。親だけのSIGINT/SIGTERM、起動途中signal、wait例外、起動直後status例外、終了確認不能、正常成功/子失敗と安全な再起動を確認。同名process一括終了はしない。 |

## 最終Linux結果

Rust全414件、関係するPython87件、toolsのPython2件が成功。各コマンドの終了コードは0。Rustの失敗/無視/測定/フィルター除外、Pythonの失敗/skipは0。追加Rust2件とPython24件はこの件数に内包され、再加算していない。fmt、Clippy（警告をエラー扱い）、locked release buildも終了0。

## 試験結果の読み方

最終件数・終了コードは `validation.json` と各 `final-*.json/log`。途中の失敗・コンパイルエラー・旧CLIによる失敗を最終成功数へ加算しない。修正前の最初のSRM試験にはテスト実装上のcompile errorとValue経由sealの誤計算があり、そのログも保持した。これらは不具合の再現根拠に使わず、型付きpayloadの正しいsealに直した `baseline-issue11-sealed.*` を根拠とした。

追加回帰はPython24 test methods（subTestの条件数を別test件数に水増ししない）とRust2 test methods。Linuxでは追加回帰と関係する既存Pythonモジュール、全Rust試験、fmt、Clippy `--locked --all-targets -- -D warnings`、locked release buildを実行する。Python全スイートの実行とは記載しない。Windowsで2つのPOSIX signal試験がskipになるのはOS制約であり、Windowsの実停止成功の証拠ではない。

`tested-sources.json` は最終試験開始前のソース指紋、`preservation.json` は起点treeとの変更範囲照合。新規証拠以外の既存interimで変更するのはIssueで名指しされた2スクリプトだけ。Fs区間、Griffiths L2、柱Tet4、杭、降雨、湿潤JOINT等の保護記録・過去の失敗・既存passed=false・低プレロードNC全撤去未達は保持する。許容値、参照値、過去の合格数値、履歴許容値は追加・変更していない。

## 残件とWindows CI

P2修正・回帰が同ブランチへ反映されたIssueだけをcloseする。反映後の実際のclose結果と新先端CI状態はIssue/PRコメントにも記録する。このコミットの作成時点で新先端のWindows CIは未確認であり、Linux成功をWindows成功へ転用しない。

旧run [42699446488](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/42699446488) はGitHubのjob結果を読み直した（`prior-windows-ci.json`）。集約failure、rust-core/rust-integration-1/verificationはsuccess、python/rust-integration-0はcancelled。cancelledを成功や内部の数値失敗に数えない。今回はIssue修正を優先し、CIのtimeout/チェックを緩めず、CI分割の追加変更は行わない。

一般の塑性・水・JOINT・構造・空間履歴収束、低プレロードNC全撤去、P3b/P4各行、P6は未完。合成fixtureの診断検証は実Rust MATSET施工全体の新しい受入を意味しない。WindowsのCIM/保持handle/実停止・Native GUI・実大規模FEMは今回未実行。SIGKILL、OSクラッシュ等の協調不能な親終了ではfinally自体が動かないので、残るlockの所有を手動確認してから扱う。

## ユーザーPC手順

既存作業・studyは変更せず、別の検証checkoutで同ブランチの新コミットを取得する。`git rev-parse HEAD^` が起点SHAと一致することを確認する。PR #13はマージしない。

```powershell
git fetch origin codex/plastic-spatial-convergence
git switch codex/plastic-spatial-convergence
git pull --ff-only origin codex/plastic-spatial-convergence
git rev-parse HEAD
git rev-parse HEAD^
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
python -m pip install -r requirements-reference.lock.txt
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo test --locked --release
cargo build --locked --release
python -m unittest discover -s tests -p 'test_p2_*.py' -v
```

PowerShellでは各コマンド直後に `$LASTEXITCODE` を確認し、非ゼロを後続の成功と混ぜない。追加Python回帰は通常実行だけで内部からnormal/-O/envの別processを起動する。POSIX signalの2試験はWindowsではskipであり、実Windows終了を確認したい場合は別の隔離環境の使い捨て親子だけを使う。既存のsafe-pause候補を実studyへ `--execute` して試さない。

生存子を終了確認できず `.reference-resume.lock` が残った場合は、自動で削除しない。記録されたparent_pid/child_pidと起動時刻・image・studyを確認し、当該所有子の終了を確かめた後だけ復旧する。名前一致のprocess一括終了は行わない。

既存sampled診断は生成時の実装hashにも結び付く。過去reportのhashを書き換えて新実装へ転用しない。新しい診断が必要なら新出力先へsampled→exactを生成し、固定済みの過去reportは保持する。

生ログの末尾空白・空行も記録として保持した。生ログを含めた任意のgit diff空白検査はexit 2（`diff-with-raw-logs.*`）、生ログを除いたコード/JSON/文書の検査はexit 0（`diff-source.*`）。この任意検査をテスト成功件数へ加算せず、元の失敗ログを整形・削除していない。CRLFの既存Pythonは改行形式を保持した。
