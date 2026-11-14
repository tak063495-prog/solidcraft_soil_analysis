# ユーザーPCでの再現

PR #13をマージせず、`bb30e468445493c965693f7590c873d3f1904626` を祖先に持つ同じブランチで実行する。既存の保護記録を変更しない。今回の対象は5 kPa追加プレロードの小さい既存NCケースである。

```sh
git fetch origin codex/plastic-spatial-convergence
git switch codex/plastic-spatial-convergence
git pull --ff-only origin codex/plastic-spatial-convergence
git merge-base --is-ancestor bb30e468445493c965693f7590c873d3f1904626 HEAD
```

Windowsでも数値観測crateの既存SHAを照合するため、LFのソースcheckoutを用いる。別の新規checkoutを作る場合は以下を使える。既存の未保存変更をresetしてはならない。

```sh
git -c core.autocrlf=false clone --branch codex/plastic-spatial-convergence https://github.com/tak063495-prog/solidcraft_soil_analysis.git soil-nc
```

Rust 1.99.1、Python 3.13、`requirements-reference.lock.txt` の環境を用意する。Cargo依存が未取得なら `cargo fetch --locked` を先に実行する。以降の `python` はその仮想環境の実行ファイルを使う。既存試験が参照する過去の証拠も省かない。

## 保存証拠を材料再求解なしで再監査

未存在の出力先を指定する。

```sh
python interim/2026-11-10-nc-mean-stress-release/remeasure.py --output work/NEW-nc-mean-remeasured
```

全アーカイブ分割、結合全体、全エントリー、ソースsnapshotのSHAとサイズを照合する。実線形系はfloat64の生ビットを前行との差分XORとgzipで無損失保存しており、元のJSONLへ復元したSHAも照合する。丸め・間引き・低精度化は行わない。全保存場から力を再組立し、材料域、凍結、荷重分割の上限、刻み差、全解放後のゼロ増分、実checkpointからの再開を検査する。

`verification.json` の `archive_verified=true`、`measurements_byte_identical=true`、`numerical_acceptance=false` を確認する。再監査成功は全上部撤去の合格ではない。特に補助荷重を含む探索経路のλ=1を、成分上限内の撤去の合格に転用しない。

Windows等で派生ノルムの最下位桁だけが異なる場合は、以下を明示する。

```sh
python interim/2026-11-10-nc-mean-stress-release/remeasure.py --output work/NEW-nc-mean-remeasured-platform --allow-platform-float-differences
```

これはアーカイブ、float64/元JSONLの全バイト、離散情報、凍結、再開の厳密一致を緩めない。力分割の再組立差を含む全派生float差と `measurements_byte_identical` の真偽を記録する。ソースのCRLF/LF差も両SHAを記録し、それ以外の差は拒否する。

## 入力から経路を新規求解

BLAS/OMP各1スレッドとUTF-8を使う。Linuxでは以下、PowerShellでは対応する `$env:変数名='値'` を設定する。

```sh
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
export PYTHONUTF8=1 PYTHONIOENCODING=utf-8
python interim/2026-11-10-nc-mean-stress-release/reproduce.py --output work/NEW-nc-mean --list-only
python interim/2026-11-10-nc-mean-stress-release/reproduce.py --output work/NEW-nc-mean --jobs 4
```

準備では既存アーカイブを全エントリー照合し、既存69ソースの観測crateを照合・locked offline releaseビルドする。独立側の400/800分割前史は入力のみから新規計算する。独立除荷終了場が得られなければ停止し、本体場で代用しない。本体の初期・載荷・除荷も新規CLIで再計算し、診断起点の除荷終了全状態と比較する。

実行するのは、補助荷重を含む平均応力分割、その局所材料RK4の刻み回復対照、独立全体荷重速度RK4、元の保持荷重の各成分上限内に制限した分割である。両順序と1/128・1/256、独立4→8積分分割、前史400→800を保存する。新しい形状や構成則は作らない。力成分の分割をclipする制御と、禁止している材料状態の降伏面への射影を混同しない。

各主経路は600または900秒、最大4,096受理増分の有界診断である。全体の逐次実行の制御予算合計は約7時間50分で、準備・ビルド時間が別に加わる。既定は4並行で、数値停止すれば早く終了する。PC速度と同時負荷により時間切れの到達点は変わる。旧到達値に合わせるために設定や場を編集しない。一時JSONLと配布証拠のための空きディスクを確保する。

途中checkpointは、その実行が実際に書いた `early.json` / `restart.json` の生バイトを渡す。モデル・ソース・数値設定・実行バイナリとペイロードのSHA束縛を保つ。LinuxのcheckpointをWindowsのバイナリへ付け替えたり、ソース変更後に再封印したりしない。Windowsでは `.exe` を自動選択する。

補助荷重経路は途中から終端まで再開する。成分上限内の経路は24受理増分の再開を比較し、未到達の完了までの一致を主張しない。粗側の時間停止からの継続は、終端が既にs=1なら不要として省く。内部材料呼出しが予算を超えた場合は外側watchdogの停止も別記する。時間切れ・点数上限・最小増分・未収束を破壊や合格に変更しない。

`report.json` / `measurements.json` の `passed`、停止理由、s、両成分解放率、残存保持荷重を確認する。sは非比例経路の座標であり、途中のsを元の比例撤去λと呼ばない。s=1かつ保持荷重ゼロの場合だけ全解放λ=1である。診断プロセスの終了コード0は数値合格を意味しない。

## Linux試験とWindows CI

```sh
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo test --locked --release
cargo build --locked --release
python -m unittest discover -s tests -p 'test_*.py' -v
python -m unittest discover -s tools -p 'test_*.py' -v
```

Rustビルド・試験スレッドは各2を使う。新規保存証拠試験も通常のPython自動探索に含まれる。Windows CIは既存workflowの `rust-core`、`rust-integration-0`、`rust-integration-1`、`python`、`verification` と集約 `Rust and Python on Windows` を、実先端と紐づけて確認する。テストmerge SHAとPR先端SHAを区別する。未起動、queued、in_progress、GitHub内部起動障害をWindows成功と呼ばない。

CIの保存証拠再監査はLinuxで保存した全場の照合であり、この新しい全数値経路をWindowsで入力から再求解した証明ではない。新規求解には上の再現ドライバを使う。
