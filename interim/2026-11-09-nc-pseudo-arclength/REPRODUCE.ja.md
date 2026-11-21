# 再計測と再実行

この資料は5 kPaの既存小ケース専用。合格・全解放・安定性を認定するツールではない。Windowsでの実行と新先端のCIは未確認。

## 求解しない再計測（Linuxで実施）

リポジトリのルートで、numpyのあるPythonを使用する。出力先は未存在のフォルダを指定する。

```sh
python interim/2026-11-09-nc-pseudo-arclength/remeasure.py --output work/nc-arc-remeasure
```

`archive-sha256.json` の順番で `evidence.tar.xz.partNN` を結合し、分割・全体・全エントリーをSHA256照合する。`out/evidence/` に展開し、現在の監査ソースの一致も検査したうえで、`out/recomputed/measurements.json` を再生成する。求解・材料積分は実行しない。元JSONとの一致は `out/remeasurement.json` に記録する。元アーカイブや元保存点を上書きしない。展開には約3～4 GB、結合した圧縮アーカイブ用にも空き容量が必要。

## 小ケースの数値診断を再実行

Rustはリポジトリ指定のtoolchain／lockfile、Pythonは既存のnumpy／scipyを使う。新しい物理パラメータ、許容値緩和、別構成則を加えない。

```sh
python interim/2026-11-09-nc-pseudo-arclength/reproduce.py --output work/nc-arc-reproduced --list-only
python interim/2026-11-09-nc-pseudo-arclength/reproduce.py --output work/nc-arc-reproduced
```

一覧確認は計算しない。実行モードは隔離crateの生成・ビルド後、記録した24ドライバを順に起動する。個々のコマンドが終了コード0でも、診断JSONの `passed=false` と停止理由を必ず読む。数値ドライバごとの時間制限・刻み上限はソースに固定し、大規模計算を起動する分岐はない。元のsrcや本体保存点を変更しない。指定先端e21778cの数値ソースと現在のsrcが異なる場合は準備時に停止する。

準備だけを照合する場合：

```sh
python interim/2026-11-09-nc-pseudo-arclength/prepare_evidence.py --output work/nc-arc-source-check --source-only
```

これは69ソースの再生成までで、材料求解はしない。本セッションではこの全バイト一致を実際に検査した。24件をこの一括ラッパーで再度全実行したとは主張しない。各数値ドライバは本セッションで個別起動した。

## Windowsで行う場合（未実施）

PowerShellでも上記Python引数をそのまま使える。`work\...` の未存在出力先と、同じlockfileを使う。`material_diagnostic.exe` の拡張子は自動選択する。Windows版バイナリはLinux版と別SHAになるため、Linuxの封印済み本体jobへ置き換えて再開しない。本セッションの再開確認はPythonの診断状態64＋64点であり、封印済み本体jobの完了試験ではない。

LinuxとWindowsではLAPACK・丸め差・線形検査の停止位置が変わり得る。値が異なる場合は生の新規出力を保存し、時間切れ、未収束、既存残差検査超過を破壊や合格へ変換しない。旧先端のCI成功は使用しない。
