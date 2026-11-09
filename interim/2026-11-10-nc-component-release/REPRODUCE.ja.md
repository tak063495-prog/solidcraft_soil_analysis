# ユーザーPCでの再計測と追加探索

この手順は既存5 kPaケースの診断専用。全解放・破壊・安定限界を認定するツールではない。元のPR #13をマージせず、`codex/plastic-spatial-convergence` の本フォルダを含む先端をcheckoutする。元の封印jobへ実行版を付け替えたり、checkpointの束縛を再封印したりしない。

## 配布証拠だけの再計測

リポジトリルートで、既存 `requirements-reference.lock.txt` の環境と指定Rust toolchainを用意する。PythonのNumPyが必要。出力は未存在フォルダを指定する。

```sh
python interim/2026-11-10-nc-component-release/remeasure.py --output work/nc-component-remeasured
python -m unittest discover -s interim/2026-11-10-nc-component-release -p 'test_*.py' -v
```

`archive-sha256.json` の順で分割を読み、各分割・全体・全エントリーのSHA/サイズと診断ソースsnapshotを照合する。結合一時ファイルやhardlinkを使わず新規展開し、全保存応力・実線形系・凍結・対応場差を再計測する。`measurements_byte_identical=true` を確認する。求解や材料積分は実行しない。証拠検査8件の成功も、未達の全解放を合格にしない。

## 元入力からの小ケース再実行

```sh
python interim/2026-11-10-nc-component-release/reproduce.py --output work/nc-component-pc --list-only
python interim/2026-11-10-nc-component-release/reproduce.py --output work/nc-component-pc
```

準備は前の45 MBアーカイブ全1,038エントリーを照合し、同じ69ファイルの観測crateを再生成・照合・locked offline releaseビルドする。Cargo依存が未取得なら、通常の `cargo fetch --locked` で先に用意する。別のソースを使って観測crateの一致検査を迂回しない。

独立の400分割RK4前史を入力から新規計算し、除荷終了場が得られない場合は停止する。native材料場を独立側の代用へ渡さない。続いて2写像×2撤去手順×2刻み、水平先行の2再開の計10診断を順に起動する。本セッションはそれぞれを並行起動したため、PCの壁時計予算で到達点が変わり得る。結果を旧数値へ合わせない。

各診断は粗300秒・細600秒、最大2,048受理増分で止める。CLIの終了コード0と `passed=false` は両立する。必ず `report.json` の停止理由、h/v、残存荷重と全場を読む。

## 証拠不足を埋める追加探索（未完了・合格にしない）

独立側の半刻みは本セッションでは時間切れ。たとえば水平先行を、実際に書かれた同じcheckpointからさらに有界に再開する：

```sh
python interim/2026-11-10-nc-component-release/study.py --map independent --order horizontal-first --reference work/nc-component-pc/reference400 --step-size .00390625 --resume work/nc-component-pc/independent-horizontal256/restart.json --output work/nc-component-pc/independent-horizontal-tail --count 4096 --budget-seconds 1800
```

鉛直先行も `--order vertical-first`、`independent-vertical256/restart.json`、新規出力先として同様に行う。元ファイルを上書きしない。これは同じ下降弧長経路の延長ではなく、今回の成分別手順の未完区間である。ただし、時間を延ばしただけで成功や到達不能を推定しない。

同じsource/model/刻み/material分割のbindingとpayload SHAをそのまま検査する。刻みや材料分割を変更してcheckpointを付け替えない。材料分割4→8や前史400→800を比較する場合は、それぞれ入力から別の前史・別経路を生成し、全履歴を比較する必要がある。この収束確認は今回未実施であり、既存値を参照合格値にしない。

λ=1が得られても、全解放後ゼロ増分の材料・接触凍結、全場釣合い、両実装の全経路刻み感度、同じ荷重定義、完了までの途中再開を揃えるまで `passed=false` を維持する。未収束・最小増分・点数上限を破壊にしない。

## Linux全試験とWindows

```sh
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release
cargo test --locked --release
python -m unittest discover -s tests -p 'test_*.py' -v
python -m unittest discover -s tools -p 'test_*.py' -v
python -m unittest discover -s interim/2026-11-10-nc-component-release -p 'test_*.py' -v
```

`CARGO_BUILD_JOBS=2`、`RUST_TEST_THREADS=2`、BLAS/OMP各1スレッド、UTF-8を用いる。sparse checkoutの場合は既存試験が参照する旧証拠を省かない。参照エラーをスキップへ変更しない。

Windows PowerShellでもPythonの引数は同じ。`reproduce.py` は `material_diagnostic.exe` を自動選択する。WindowsとLinuxのバイナリSHAが違うので、Linux診断のnative checkpointをWindowsバイナリへ付け替えず、Windowsで入力から開始する。本セッションでWindows版の成分別経路を求解したとはしない。

push後は GitHub Actions の `FEM build and verification` で先端SHAに紐づく `rust-core`、`rust-integration-0`、`rust-integration-1`、`python`、`verification` の5群と集約を確認する。PRのテストmerge SHAとPR先端SHAを区別する。最新先端が queued/in_progress の間は成功とは書かない。
