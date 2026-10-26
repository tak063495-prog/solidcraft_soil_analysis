# FEM本体の新規Windows環境での検査

混合MCCと梁・板の連成を加えたローカル版では、Rust全249件・Python全94件、fmt/clippyが合格した（`runs/mixed-structures-20261026-v1-*`）。以下の公開f12d374に対するWindows成功とは版を区別する。SolidCraft総合CIの未達も継続して保持する。

`.github/workflows/fem-ci.yml` は `windows-2025` のGitHubホスト環境で、固定Rust 1.99.1、Python 3.13、既存の固定依存ファイルを使う。mainへのpush、pull request、手動実行を対象とする。

fmt、全ターゲットclippy、releaseの全Rust試験、CLI build、Pythonの全unittestを順に実行する。Pythonの曲面Tet10検査がRust CLIを呼ぶため、通常の実行ファイルも明示的にビルドする。ログとソース指紋を14日間のartifactに残す。Actions自体は確認したリリースのコミットSHAで固定し、ソース取得後のGit認証情報を保持しない。

この検査はFEMリポジトリを対象とする。SolidCraftの全体CI、実アプリでの画面操作、長時間の独立3D解析と空間収束、全P0〜P6の完了検証を含めたものではない。初回のリモート実行成功を確認するまでは、新規環境での再現に合格したとは扱わない。

2026-10-26、公開コミット `f12d374b01d108fa96950d4215279cbd31c37923` の[実行40235522373・試行2](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/40235522373)が成功した。Windowsの新規環境でfmt・clippy・CLIビルド・Rust247件・Python67件（無視0）を確認し、artifactのソース562件と固定依存を照合した。153ファイルはバイト一致、409ファイルはGitのcheckoutによるCRLF化だけで一致した。後続のローカル再開ツール・内部壁対応はこのコミットのCI対象には含まれない。記録は `runs/goal-resume-20261026-v1/clean-windows-ci-verification.json`。

通常の利用環境は [Setup.ps1](../scripts/Setup.ps1)、[README](../README.md) を参照。CIはPythonのメッシュ検査に必要な依存を導入するが、公式SolidCraftアプリのダウンロード・GUI起動は行わない。
