# ユーザーPCでの再現・途中再開

PR #13 はマージせず、指定起点 `5bee3b47f61878e614cae1e7268a57449d2f12d1` を祖先に持つ同じブランチを使う。既存の未保存変更を reset しない。

```sh
git fetch origin refs/heads/codex/plastic-spatial-convergence:refs/remotes/origin/codex/plastic-spatial-convergence
git switch codex/plastic-spatial-convergence
git pull --ff-only origin codex/plastic-spatial-convergence
git merge-base --is-ancestor 5bee3b47f61878e614cae1e7268a57449d2f12d1 HEAD
```

ローカルにこのブランチがまだない場合は、`git switch --no-track -c codex/plastic-spatial-convergence origin/codex/plastic-spatial-convergence` を使う。別 checkout を作る場合は LF を保持する。

```sh
git -c core.autocrlf=false clone --branch codex/plastic-spatial-convergence https://github.com/tak063495-prog/solidcraft_soil_analysis.git soil-nc
```

Rust 1.99.1、Python 3.13 と既存の `requirements.lock.txt` / `requirements-reference.lock.txt` の環境を使う。以下の `python` はその仮想環境を指す。Cargo 依存がなければ先に `cargo fetch --locked` を実行する。大きい既存リポジトリに加え、一時展開・生成場・ビルドの空き領域を確保する。

## 保存証拠の再監査

未存在の出力先を指定する。新しい MCC 求解は行わない。

```sh
python interim/2026-11-10-nc-bounded-refinement/remeasure.py --output work/NEW-nc-bounded-replay
python interim/2026-11-10-nc-bounded-refinement/audit_linear_roundoff.py --evidence work/NEW-nc-bounded-replay/evidence --output work/NEW-nc-bounded-linear-roundoff.json
```

全分割・結合全体・全エントリーの SHA とサイズを照合し、封印したソースを隔離した別 Python プロセスで実行する。全受理材料場・接触履歴・変位、未 commit の内部速度場、元の成分荷重境界、残存域釣合い、棄却時の状態保持、終盤再開、刻み差を再計測する。全線形求解の右辺・解・行列 SHA・残差と、各プロセス末尾128件の実行列を保存している。末尾より前の実行列全量の再構築は、この新アーカイブの検証範囲に含めない。

`verification.json` の `archive_verified=true`、`measurements_byte_identical=true`、`numerical_acceptance=false` を確認する。厳密な派生 JSON のバイト比較には Linux / WSL と記録した NumPy 環境を使う。別プラットフォームのノルム最下位桁の違いは、アーカイブ破損や全解放合格へ読み替えず、その違いを別記する。

追加の線形残差監査器は、アーカイブ作成後に追加した別添ソースである。全求解要約を受理/棄却列へ対応させ、元の配列と連続配列での残差を区別する。末尾の実行列では Fortran 配列順序と拡張精度の計算も比較する。全実行列がない範囲へ拡張精度監査を一般化しない。

旧 f19ce64 証拠の全再監査は次の既存コマンドで行う。修正後は旧アーカイブ内のソースを使い、現行 checkout との差を記録する。これは旧証拠の再現であり、現行生産コードの数値認定ではない。

```sh
python interim/2026-11-10-nc-mean-stress-release/remeasure.py --output work/NEW-historical-mean-replay
```

## 入力から新規に求解

BLAS/OMP を各1スレッド、Rust ビルド・試験を各2スレッドにする。PowerShell では対応する `$env:変数名='値'` を設定する。

```sh
export CARGO_BUILD_JOBS=2 RUST_TEST_THREADS=2
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
export PYTHONUTF8=1 PYTHONIOENCODING=utf-8
cargo build --locked --release
python interim/2026-11-10-nc-bounded-refinement/reproduce.py --output work/NEW-nc-bounded --list-only
python interim/2026-11-10-nc-bounded-refinement/reproduce.py --output work/NEW-nc-bounded --jobs 4
python interim/2026-11-10-nc-bounded-refinement/measure.py --evidence work/NEW-nc-bounded --output work/NEW-nc-bounded/measurements.json
```

対象は既存 `make_free_contact.py --mcc --pc 100 --preload-kpa 5 --force 1 .5 0` だけである。準備は既存1,038エントリーのアーカイブを照合する。観測器は f19ce64 の Rust ソースと既存追加観測コードを使い、既存69ソースの全 SHA と照合して locked/offline release ビルドする。現行 CLI は同じ5 kPaの初期・載荷・除荷を新規計算し、除荷終了全状態を照合する。独立400/800分割の前史は入力のみから新規計算する。独立除荷終了場が得られなければ、本体の場で代用しない。

新経路は各900秒、最大4,096受理増分、最小刻み `2^-28` の有界診断である。10主経路と6再開経路の時間予算合計は最大約4時間、4並行では約1時間＋準備・ビルドが目安となる。数値停止なら短縮し、PC速度・同時負荷で時間停止点は変わる。プロセス終了コード0は診断処理の終了であり、`passed=true` を意味しない。最小刻み・未収束・時間切れ・外側 watchdog をそれぞれ記録する。

準備だけを先行した場合に限り、出力先に完了した `prepared/` だけがある状態で `--reuse-prepared` を使える。他の結果を上書きするオプションではない。

## 実際の途中 checkpoint から再開

`tail.json` は差分先行では s≥.97、平均先行では s≥.515 の最初の実保存点の生バイトコピーである。`early.json`、`restart.json` も実保存点である。`--list-only` が示す対象経路と同一の数値設定行を使い、出力先だけを未存在の兄弟ディレクトリに変更して `--resume <実checkpoint>` を追加する。

入力・ソース・数値設定・実行バイナリの束縛とペイロード SHA を保持する。ソース変更後の再封印、最小刻みや反復改良設定の付け替え、Linux バイナリの封印を Windows `.exe` へ付け替える操作はしない。Windows では新規の `.exe` と自身の checkpoint を作る。NumPy と `longdouble` の仮数ビット数も束縛に含まれる。

自動ドライバは本体/独立・粗/細・両順序の6保存点を別プロセスで再開し、実際の停止点まで照合する。完了に到達していない経路については、完了までの再開一致を主張しない。

## 既存試験と Windows CI

```sh
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo test --locked --release
cargo build --locked --release
python -m unittest discover -s tests -p 'test_*.py' -v
python -m unittest discover -s tools -p 'test_*.py' -v
```

新しい重い撤去ドライバは既定試験へ追加していない。Windows CI の workflow・ジョブ群・45分制限も維持する。f19ce64 では `python` と `rust-integration-0` が45分制限で取り消された過去がある。先端の CI は実際の先端 SHA と対応させて確認し、queued / in_progress / 未取得なら **未確認** とする。Linux 成功、旧先端の成功、数値未達の診断プロセス終了を、今回の Windows 成功や全解放合格へ転用しない。
