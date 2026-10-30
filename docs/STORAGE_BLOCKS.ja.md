# 元結果を保持する分割読込

この候補では、projected MCCの元結果を小さな履歴ブロックへ分け、既存の力学検査へ全工程・全積分点・初期種状態・施工遷移・最終確定状態を渡します。元結果と各ブロックのSHA256を確認します。表示用の平均値で材料履歴を置き換えません。

元結果のencoded読込はprojected MCCに限り32 GiBまでです。各componentのencoded上限256 MiB、payloadとキャッシュの会計上限512 MiB、表示パケットの上限は保持します。監査側のNumPy配列や保持する参照は別であり、プロセス全体のメモリーを512 MiBに保証するものではありません。MC/IDH等の従来の受付上限は変えません。

実結果1,459,184,313 bytes・8,061 Tet10・217,647 GP・全5工程で、元の混合MCC検査と既存の全有効GP/板/節点VTU比較を通過した読込実装が基礎です。既存APIの板座標tupleを保つ修正を含む、この統合候補の品質検査・実producerと大きな14,909要素ケースは別に検証します。Native、一般10万要素、空間収束、正式配布や全P6の受入はこの記録だけで認定しません。

通常の環境構築はrequirements.lock.txt、独立参照を含む試験はrequirements-reference.lock.txtを使用します。追加のijsonは3.6.0に固定しています。Windows x64・CPython3.13のwheelだけを使う検証では、requirements-verification-win-cp313.lock.txtに9パッケージのwheel SHA256を固定しています。

検証用の9 wheelをwheelhouseへ保存し、新しい仮想環境へ次のようにオフライン導入します。この入口はCPython 3.13.14・Windows x64に固定しています。既存の環境パスは拒否し、変更しません。

```powershell
.\scripts\Setup-Verification.ps1 -EnvironmentPath '.\new verification env' -Wheelhouse '.\wheelhouse' -Python python
& '.\new verification env\Scripts\python.exe' -W error::SyntaxWarning -m unittest discover -s tests -p 'test_*.py' -v
```

セットアップは、9 wheelのサイズとSHA256、ロックとの対応を導入前に検査し、`--no-index --require-hashes --no-cache-dir`で導入します。pipの導入元がすべてローカルの検証済みwheelであること、全9パッケージのバージョンと実際のimport先、新規環境のprefixとuser-site無効化、導入後の元wheel・ロック・スクリプト不変を確認します。証拠は新規環境内の`verification-environment.json`、`pip-install-report.json`、`pip-install.log`、`pip-check.log`に保存します。

一般的なRust/CAD環境の構築は従来の`Setup.ps1`を使います。この検証用入口はPythonの固定依存環境を作成するもので、Rust・CAD・GPUの配布受入を代替しません。

Rustを含む試験には固定したFEM実行ファイルをtarget/releaseへ置くか、同じソースからビルドします。別の場所でビルドしたexeのハッシュを元の受入証明書へ書き直さず、新しいジョブと監査を作成してください。同一PCの新規配置と新規PC・別GPUの受入は区別します。

pipの環境変数と設定ファイルは導入・整合検査の子プロセス内で無効化し、既存の設定による導入先の変更を防ぎます。親プロセスやユーザーのpip設定ファイルは変更しません。
