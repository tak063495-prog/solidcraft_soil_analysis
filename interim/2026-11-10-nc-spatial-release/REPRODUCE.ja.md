# ユーザーPCでの再現

既存5 kPaケースの診断専用。全解放・安定限界・破壊を認定するコマンドではない。PR #13をマージせず、同じブランチの本フォルダを含む先端で実行する。

```sh
git fetch origin codex/plastic-spatial-convergence
git switch codex/plastic-spatial-convergence
git pull --ff-only origin codex/plastic-spatial-convergence
git merge-base --is-ancestor fdd3018e59bfa04a7808f59ce0d19db3ca5d72aa HEAD
```

リポジトリ指定のRust 1.99.1、Python 3.13と `requirements-reference.lock.txt` の環境を用意する。Cargo依存が未取得なら `cargo fetch --locked` を先に実行する。以下の `python` はその仮想環境の実行ファイルを使う。PowerShellなら `.venv/Scripts/python.exe`、Linuxなら `.venv/bin/python` へ置き換えられる。

## 求解せず配布証拠を再監査

未存在の出力先を指定する。

```sh
python interim/2026-11-10-nc-spatial-release/remeasure.py --output work/NEW-nc-spatial-remeasured
```

分割アーカイブ・結合全体・全エントリー・ソースsnapshotのSHA/サイズを照合して展開し、全保存応力から力を再組立する。全場差、凍結、実線形系、別プロセス再開も再計測する。出力の `verification.json` で `measurements_byte_identical=true`、`numerical_acceptance=false` を確認する。SHAやJSONの不一致を、許容値緩和・スキップ・値の書換えで回避しない。

結合アーカイブの複製は作らず、分割を逐次読み出す。JSONLは元バイトのSHAを検査してからgzip形式で展開し、再監査時に読み戻す。再計測用の追加ディスクを抑えるための処理で、保存場・試行を省かない。

Windows等の別プラットフォームでは、同じ保存データに対する派生ノルムの最下位桁が異なり得る。必要な場合は `--allow-platform-float-differences` を明示する。アーカイブの全バイト、離散情報、全場の力基準、凍結・再開の厳密一致は引き続き要求し、派生floatの全差分と `measurements_byte_identical` の真偽を保存する。新しい誤差許容値・合格条件は設けない。CIのWindows証拠試験はこのモードで全場を再監査し、差分をログへ出す。ソースcheckoutのCRLF/LF差だけも両SHAを記録し、その他のソース差は拒否する。

## 入力から数値経路を再実行

BLAS/OMP各1スレッド、UTF-8を使う。PowerShellでは `$env:OPENBLAS_NUM_THREADS='1'` 等、Linuxでは `export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1` を設定する。

```sh
python interim/2026-11-10-nc-spatial-release/reproduce.py --output work/NEW-nc-spatial --list-only
python interim/2026-11-10-nc-spatial-release/reproduce.py --output work/NEW-nc-spatial
```

準備では従来の `2026-11-09-nc-turning-path-diagnosis` のアーカイブを全エントリー照合し、同じ69ソースの観測crateを生成・照合・locked offline releaseビルドする。独立400分割前史も入力から新規計算する。除荷終了場が得られなければ停止し、本体場で代用しない。

続いて2写像×2空間順序×2撤去刻み、独立左側先行の材料8分割・イベント刻み下限対照、12受理点の再開2件、時間切れ2経路からの継続2件を順に起動する。中心側先行の独立経路は、実行中に生成した途中checkpointをチェックサム照合してそのままコピーする。終了後の場からcheckpointを捏造・再封印しない。途中checkpointが生成されなければ、その再開を完了と書かない。

各主経路は粗600秒・細900秒、最大2,048受理増分。対照は600秒、短い再開は300秒/12増分、継続は600秒/1,024増分の上限である。本セッションでは複数経路を並行実行した。PCの性能・並行負荷によって時間切れの到達点は変わる。逐次実行の時間予算合計は最大約150分に準備時間が加わるが、数値停止した経路は先に終了する。旧到達値に合わせるために場や設定を変更しない。

全 `report.json` の `passed`、停止理由、s/a/b、残存保持荷重を確認する。終了コード0は観測処理の正常終了であり、数値合格ではない。独立左側先行のRK4停止を物理的限界に変換しない。sは元の比例撤去λではない。

本番CLIの全4工程途中再開は、別の新規フォルダで実行する。

```sh
cargo build --locked --release
python tools/check_low_preload_nc_restarts.py --solver target/release/soil-fem-3d --output work/NEW-nc-phase-restarts
```

Windowsではソルバーを `target/release/soil-fem-3d.exe` とする。数値経路の再現ドライバは `material_diagnostic.exe` を自動選択する。各OSで入力から開始し、Linuxのnative診断checkpointをWindowsバイナリへ付け替えない。同じソース・入力・数値設定・実行版で作られたcheckpointだけを使い、変更は別入力・別実行として始める。

## 試験とCI

```sh
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo test --locked --release
cargo build --locked --release
python -m unittest discover -s tests -p 'test_*.py' -v
python -m unittest discover -s tools -p 'test_*.py' -v
```

Rustのビルド・試験スレッドは各2を使う。既存試験が参照する過去の証拠をsparse checkoutから省かない。今回の新規Rust2/Python14件は既存5群CIの自動探索へ含まれる。新規保存証拠4件は材料再求解をせず、配布証拠を展開して再監査する。

`FEM build and verification` の実先端SHAに紐づく `rust-core`、`rust-integration-0`、`rust-integration-1`、`python`、`verification` と集約 `Rust and Python on Windows` を確認する。テストmerge SHAとPR先端SHAを区別し、queued/in_progressを成功と書かない。

全解放、全解放後の凍結、独立塑性3Dの全撤去経路と両側の刻み感度、完了までの再開が揃わない限り `passed=false` を維持する。
