# Windows CI：固定末尾を作る試験fixtureの改行を修正

**原因は試験用一時文書のtext書込みによるLF→CRLF変換。監査のバイト比較を維持し、fixtureを明示したUTF-8バイトで作る。崩壊分類D、NC-01／P6の未達は維持する。**

起点はPR #13の最新先端 `092f1322ffce293edb3b93b03d0b7649e55f0641`。同じ `codex/plastic-spatial-convergence` へpushし、PRをマージしない。前回の測定、固定SHA、失敗記録、入力・参照値は変更しない。

## 失敗したCIと原因の検証

[FEM run 43373828275](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/43373828275) の `python-core-0 on Windows`（job 131715567090）を実際に取得した。341件、error 1。`CollapseAppendixTests.test_both_pinned_suffixes_checked_before_old_prefix` の `NC collapse document appendix suffix differs: docs/LOW_PRELOAD_NC.ja.md` を確認した。[取得したログ](validation/retained-windows-python-core-0.log)を保持する。ログはコネクタが返した文字列をUTF-8バイトへ戻した保存表現で、未取得のZIP生バイトのSHAとは区別する。

[取得した全job状態](ci-failure.json)ではFEMのmatrix 7件success・通常0群1件failure、全群成功を要求する集約jobもfailureだった。[gate run 43373828255](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/43373828255) のpwsh／powershell両jobはsuccess。集約jobの失敗を独立した数値不具合として数えず、依存失敗を正しく伝播した結果として保持する。

CIログでは、checkout前の `git config --global core.autocrlf false` は成功し、`PYTHONUTF8=1`、`PYTHONIOENCODING=utf-8` も指定されていた。checkoutが文書をCRLFへ変えたという仮説は今回の原因にしない。

該当fixtureは一時文書を `Path.write_text('original prefix\n'+old_suffix+new_suffix)` で作っていた。Windowsではtext出力の既定newline変換により、文字列中のLFがCRLFになる。Gitのautocrlf設定とは別の処理である。監査は固定manifest中のLFを含むUTF-8末尾と生バイトで比較するため、この変換を正しく拒否した。

[再現測定](measurements.json)は、Linuxでも `write_text(..., encoding='utf-8', newline='\r\n')` を使ってWindowsと同じ変換を強制した結果である。

|一時文書の作成方法|固定監査の結果|
|---|---|
|LFを明示したtext出力|受入|
|CRLFを明示したtext出力|同じcollapse suffixエラーを再現|
|UTF-8のwrite_bytes|受入。LF版と同じ全バイト・SHA|
|旧suffixだけCRLF、新suffixは固定LFのまま|新suffix検査は通るが、旧suffixとの連結を検査するphysical suffix監査で拒否|

各文書のraw SHA-256、CRLF/LF数、診断用にCRLFをLFへ戻した場合のSHAを両方保存した。CRLF版のraw SHAはLF版と異なり、診断用LF化後のSHAは一致する。**この正規化は原因の照合にだけ使い、監査の比較には使わない。** パス区切りや文字コードを補正して合格させていない。

## 修正と検出力

変更は `tests/test_nc_collapse.py` のfixtureと、旧suffix改変試験の文書書込みを `write_bytes(text.encode('utf-8'))` にするもの。期待する固定バイトをOSに依存せず作る。旧suffix改変試験も改行の副作用に先に遮られず、意図した旧suffixの改変を検査できる。

追加2試験は、Windowsのtext変換をどのOSでも再現して新suffixのCRLF化を拒否すること、新suffixが正しくても旧suffixだけのCRLF化を拒否することを確認する。既存の追記改変・manifest再封印・manifest欠落・旧suffix改変の拒否試験も保持した。診断単体試験は19→21件で成功した。

`collapse_document_appendices`、`physical_document_appendices`、`audit` の関数本文は起点と同じ。[関数本文・固定文書・manifestのSHA](validation/preserved-bytes.tsv)を保持する。P6追補台帳を追加登録するためのSHA定数更新以外、監査器の処理を変更しない。旧manifestの固定SHA、旧文書、旧登録のprefixを保持し、比較の正規化・許容値緩和・参照値変更は行わない。

## 崩壊分類と未達の保持

分類は **D（判定不能）**。前回の[不均一場の診断](../2026-11-13-nc-lambda-collapse-classification/README.ja.md)を変更していない。

|根拠|維持する判断|
|---|---|
|fold約0.95140010、荷重方向と非直交な近傍弱モード、負の二階仕事・音響指標|局所的な限界点候補を支持する|
|Kは数値rank54、拡大系rank55、厳密零点・別枝・メッシュ収束は未確定|A物理的崩壊、B分岐、C純粋な数値失敗のいずれも確定しない|
|固定pc・GP離散モデルのλ_LB約0.95191145、弱モード上側約0.95227175|foldはこの区間より下。連続体FELA／硬化MCCの厳密上下界は未達|
|材料点99点は応力rayを6 kPa延長、63点は局所逆解法未確定|材料限界到達GP数は未確定のまま|

仕様の(a)(b)(c)は実装せず、現時点では(c)未達維持が根拠に合う。NC-01・P6を完了にしない。今回のCI修正を材料・求解の数値的改善や崩壊分類の変更へ読み替えない。

## 固定SHA登録と検証

このフォルダの3 JSONを先行する不変コミットへ固定し、既存 `register_sources.py` によりP6-01追補台帳へ追加登録する。固定commit・Git blob SHA-1・SHA-256・bool/null値を全て検査し、監査の除外や未登録資料の黙認をしない。

`validation.json` は最終全件前のsnapshotであり、成功へ書き換えない。Linuxの全Python・Rust、fmt・clippyの最終件数・実終了コード・ログSHAは `validation/final-summary.tsv` に記録する。前回092f132のWindows失敗は保持し、新しい先端のCI結果と区別する。

再現は既存の環境を用いて以下を実行する。Linuxでは既存PowerShell互換試験にもpwshを利用できるPATHを指定する。監査や試験をskipして通さない。

```text
python -m unittest discover -s tests -p test_nc_collapse.py -v
python tools/run_ci_python_tests.py --group python-core-0 --output runs/NEW-linux-core-0
python tools/run_ci_python_tests.py --group python-core-1 --output runs/NEW-linux-core-1
python tools/run_ci_python_tests.py --group python-mean-evidence --output runs/NEW-linux-mean
python tools/run_ci_python_tests.py --group python-spatial-evidence --output runs/NEW-linux-spatial
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo test --locked --release --all-targets
cargo test --locked --release --doc
```
