# FEM本体の新規Windows環境での検査

`.github/workflows/fem-ci.yml` は `windows-2025` のGitHubホスト環境で、固定Rust 1.99.1、Python 3.13、既存の固定依存ファイルを使う。mainへのpush、pull request、手動実行を対象とする。

fmt、全ターゲットclippy、releaseの全Rust試験、CLI build、Pythonの全unittestを順に実行する。Pythonの曲面Tet10検査がRust CLIを呼ぶため、通常の実行ファイルも明示的にビルドする。ログとソース指紋を14日間のartifactに残す。Actions自体は確認したリリースのコミットSHAで固定し、ソース取得後のGit認証情報を保持しない。

この検査はFEMリポジトリを対象とする。SolidCraftの全体CI、実アプリでの画面操作、長時間の独立3D解析と空間収束、全P0〜P6の完了検証を含めたものではない。初回のリモート実行成功を確認するまでは、新規環境での再現に合格したとは扱わない。

通常の利用環境は [Setup.ps1](../scripts/Setup.ps1)、[README](../README.md) を参照。CIはPythonのメッシュ検査に必要な依存を導入するが、公式SolidCraftアプリのダウンロード・GUI起動は行わない。
