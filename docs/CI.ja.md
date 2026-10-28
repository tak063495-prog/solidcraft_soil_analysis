# FEM本体の新規Windows環境での検査

公開7b6f3e1650bba94db9c466c07b8e159f687b8fdfの[Windows CI40643394162](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/40643394162)は成功し、artifact12797471087のZIP SHA・600ソース・Rust265件・Python137件・固定依存とRust1.99.1を照合しました。174ファイルはバイト一致、426ファイルはcheckoutのCRLF化だけで一致します。初回ea3126の8.3短縮パス比較1件の失敗と修正を別記録に保持しています。後続607ソースの公開f7aaa1fも[CI40645436331](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/40645436331)が成功し、artifact12799530192のZIP SHA・607ソース・Rust265件・Python159件・固定依存とRust1.99.1を照合しました。178ソースはバイト一致、429ソースはcheckoutのCRLF化だけで一致します。Native総合CI6327件と実GPUの限定受入はFEM専用CIと別に記録します。

以下は過去版の記録です。


後続公開b4fb3d1312f98536e3cace51842d33fbfbf83947の[Windows実行40427004281](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/40427004281)は成功した。artifact12701391682のZIP SHA、598ソース、Rust263件・Python124件（無視0）、Rust1.99.1と固定Python依存を照合済み。170ファイルはバイト一致、428ファイルはcheckoutのCRLF化だけで一致した。記録は `runs/fem-ci-b4fb3d-success-20261027-v1/verification.json`。後続q6のMCC試行・板対応/3はローカル263/128件と新バイナリの独立200→400比較を通過しているが、新公開コミットのリモートCIとは区別する。Native総合CIと実GPU受入は別報告に保持する。

公開dfb7c29のWindows実行40304406676は、Rust249件・fmt/clippy・CLI buildが合格し、Python94件のうち履歴参照再開の1件が失敗した。元ログとartifactは `runs/fem-ci-dfb7-failure-20261026-v1` に保持する。履歴ソースのパスを文字列だけで照合する処理を、実ファイルの正規化したパス＋変更しないバイトSHAへ修正し、別表記の同一ファイルを検査する回帰試験を追加した。修正f962806の[Windows CI40309767151](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/40309767151)は成功した。ダウンロードしたartifactのRust249件・Python94件（無視0）、Rust1.99.1とPython固定依存、ソース578件を公開時の原本と照合した。157件はバイト一致、421件はcheckoutのCRLF化だけで一致する。記録は `runs/fem-ci-f962-success-20261026-v1/verification.json`。後続のローカルPython・native構造部材表示の変更はこの実行の対象外であり、SolidCraft総合CIの未達も保持する。

混合MCCと梁・板の連成を加えたローカル版では、Rust全249件・Python全94件、fmt/clippyが合格した（`runs/mixed-structures-20261026-v1-*`）。以下の公開f12d374に対するWindows成功とは版を区別する。SolidCraft総合CIの未達も継続して保持する。

`.github/workflows/fem-ci.yml` は `windows-2025` のGitHubホスト環境で、固定Rust 1.99.1、Python 3.13、既存の固定依存ファイルを使う。mainへのpush、pull request、手動実行を対象とする。

fmt、全ターゲットclippy、releaseの全Rust試験、CLI build、Pythonの全unittestを順に実行する。Pythonの曲面Tet10検査がRust CLIを呼ぶため、通常の実行ファイルも明示的にビルドする。ログとソース指紋を14日間のartifactに残す。Actions自体は確認したリリースのコミットSHAで固定し、ソース取得後のGit認証情報を保持しない。

この検査はFEMリポジトリを対象とする。SolidCraftの全体CI、実アプリでの画面操作、長時間の独立3D解析と空間収束、全P0〜P6の完了検証を含めたものではない。初回のリモート実行成功を確認するまでは、新規環境での再現に合格したとは扱わない。

2026-10-26、公開コミット `f12d374b01d108fa96950d4215279cbd31c37923` の[実行40235522373・試行2](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/40235522373)が成功した。Windowsの新規環境でfmt・clippy・CLIビルド・Rust247件・Python67件（無視0）を確認し、artifactのソース562件と固定依存を照合した。153ファイルはバイト一致、409ファイルはGitのcheckoutによるCRLF化だけで一致した。後続のローカル再開ツール・内部壁対応はこのコミットのCI対象には含まれない。記録は `runs/goal-resume-20261026-v1/clean-windows-ci-verification.json`。

通常の利用環境は [Setup.ps1](../scripts/Setup.ps1)、[README](../README.md) を参照。CIはPythonのメッシュ検査に必要な依存を導入するが、公式SolidCraftアプリのダウンロード・GUI起動は行わない。
