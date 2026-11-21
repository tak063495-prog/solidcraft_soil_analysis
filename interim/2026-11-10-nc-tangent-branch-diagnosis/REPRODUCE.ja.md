# ユーザーPCでの再監査・再計算・途中再開

PR #13 はマージしない。既存の未保存変更をresetせず、指定起点を祖先に持つ同じブランチを使う。

```sh
git fetch origin refs/heads/codex/plastic-spatial-convergence:refs/remotes/origin/codex/plastic-spatial-convergence
git switch codex/plastic-spatial-convergence
git pull --ff-only origin codex/plastic-spatial-convergence
git merge-base --is-ancestor 113d83b8a37b3b68ce46eb00b6664f042334eae0 HEAD
```

初回cloneではLFを保つ。

```sh
git -c core.autocrlf=false clone --branch codex/plastic-spatial-convergence https://github.com/tak063495-prog/solidcraft_soil_analysis.git soil-nc
cd soil-nc
```

大きい既存証拠を含むリポジトリである。今回クラウドではソースを指定先端で先にチェックアウトし、必要な証拠を順に取得した。疎checkoutを使う場合も、再監査するアーカイブを省略したまま試験を実行しない。欠落は失敗であり、スキップや合格にしない。

Rust 1.99.1、Python 3.13、既存 `requirements-reference.lock.txt` を使用する。今回の計測環境はPython 3.13.17、NumPy 2.5.4、SciPy 1.18.1。WindowsでLinuxの実行版SHAを `.exe` に付け替えない。各PCで自身の実行版・新規checkpointを作る。

```sh
export CARGO_BUILD_JOBS=2 RUST_TEST_THREADS=2
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
export PYTHONUTF8=1 PYTHONIOENCODING=utf-8
```

PowerShellでは各変数を `$env:変数名='値'` で設定する。以下の `python` は固定環境のものを指す。必要なら先に `cargo fetch --locked` を行う。

## 保存証拠だけから再監査

新しいMCC釣合い求解は行わない。出力先は未存在にする。

```sh
python interim/2026-11-10-nc-tangent-branch-diagnosis/remeasure.py --output work/NEW-tangent-replay
```

新アーカイブの全分割・結合全体・全エントリーをSHA照合する。次に前回の2,034,108,934バイト・1,349エントリーも全SHA照合し、今回診断が使った入力証明と一致させる。新アーカイブ内の封印ソースを隔離プロセスで使い、弧長全場の釣合い・成分境界・凍結・刻み差・途中再開・全拡大線形系を再計測する。JOINTの28/56受理場の履歴伝播も再生する。

`verification.json` で `archive_verified=true`、`previous_archive_verified=true`、`measurements_byte_identical=true`、`joint_propagation_checks_exact=true`、`numerical_acceptance=false` を確認する。厳密な派生JSONのバイト比較は記録環境と一致するLinux/WSLを推奨する。WindowsのBLAS・固有値・ノルムの最下位桁が異なれば、この監査器は差を隠さず失敗させる。生証拠のSHA照合と、派生浮動小数のバイト一致を区別し、差を全解放合格に読み替えない。

主アーカイブ封印後の別添監査は、上でSHA照合した同じ生応答を使う。

```sh
python interim/2026-11-10-nc-tangent-branch-diagnosis/audit_partitions.py --tangent work/NEW-tangent-replay/evidence/diagnosis/tangent/measurements.json --controls work/NEW-tangent-replay/previous-evidence/native-complement-refined128/controls.json --output work/NEW-partition-global.json
```

この追加器は両入力と既存の保持荷重定義ソースのSHAを照合し、固定分割・中央/片側差分が弱い全体モードへ与える影響だけを再計測する。主アーカイブ内の封印ソースへ後付けしたとはしない。

生証拠はアーカイブ内 `diagnosis/`、JOINTは `joint-history/`。詳細な全GP・左右特異ベクトル・実行列固有値は `diagnosis/tangent/measurements.json` にある。本体の固定履歴要求・応答、各弧長の全状態・全線形系、実checkpointも保存する。旧保存場自体は前回のアーカイブを参照し、今回の配布物へ重複格納しない。

## 指定ケースを新規に診断する

```sh
cargo build --locked --release
python interim/2026-11-10-nc-tangent-branch-diagnosis/extract.py --output work/NEW-tangent-prior
python interim/2026-11-10-nc-tangent-branch-diagnosis/prepare_observer.py --output work/NEW-tangent-observer
python interim/2026-11-10-nc-tangent-branch-diagnosis/reproduce.py --evidence work/NEW-tangent-prior --observer work/NEW-tangent-observer --output work/NEW-tangent-study --solver target/release/soil-fem-3d
python interim/2026-11-10-nc-tangent-branch-diagnosis/measure.py --source work/NEW-tangent-study --evidence work/NEW-tangent-prior --output work/NEW-tangent-study/measurements.json
python interim/2026-11-10-nc-tangent-branch-diagnosis/joint_study.py --output work/NEW-joint-history
```

Windowsでは `--solver target/release/soil-fem-3d.exe` を使う。観測器は起点から変えていない現在のRustソースの隔離コピーへ、既存の読み取り専用追加観測コードだけを足す。本番srcへ追加しない。物理入力は5 kPaプレロードの既存ケースであり、100 kPa追加プレロードへ変えない。

診断は実残存域λ=0の非零初期推定Newton、本体の初期/載荷/除荷の新規求解、保存実行列/材料の再評価、5弧長主経路、2途中再開を順に行う。各主経路は300秒・128または256点、2並行。数値停止なら短くなる。遅いPCでは時間予算停止になり得る。停止理由をそのまま保持し、時間切れや点数上限を破壊・静的限界・合格にしない。プロセス終了0は診断処理の終了であって `passed=true` ではない。

全解放候補が将来得られても、次を全て満たすまで合格にしない。

|必要な証拠|今回の状態|
|---|---|
|元保持荷重の成分境界|今回の全弧長受理場・試行で検査|
|残存域λ=0の実Newton釣合い|本体・独立で確認|
|撤去中の全材料・接触履歴と棄却時状態保持|今回到達範囲のみ。旧全前史は出所とSHAを保持|
|全解放λ=1・実ゼロ保持荷重|未達|
|全解放後の実ゼロ増分による全材料・接触履歴凍結|未達|
|入力のみからの独立塑性3D完了全場|未達|
|両実装の完了までの刻み感度|未達。今回の同じ弧長での差だけでは代用不可|
|実checkpointから完了までの別プロセス再開|未達。今回の2再開は部分経路の終端まで|

## 実弧長checkpointから途中再開

新規の本体・細刻み経路の例。出力先だけを新しい兄弟ディレクトリにする。

```sh
python interim/2026-11-10-nc-tangent-branch-diagnosis/arc.py --evidence work/NEW-tangent-prior --binary work/NEW-tangent-observer/target/release/examples/material_diagnostic --output work/NEW-tangent-native-resumed --map native --step-size 0.00005 --count 256 --material-steps 4 --budget-seconds 300 --resume work/NEW-tangent-study/arc-native-5e-05-m4/middle.json
```

独立側は `--map independent` と対応する `arc-independent-5e-05-m4/middle.json` を使う。Windowsの観測器には `.exe` を付ける。保存したdriver SHA・実行版・起点・数値設定の束縛とペイロードSHAを維持する。途中で設定やソースを変えて再封印しない。この途中再開は、生産CLIの全撤去完了の再開認定ではない。

## Linux・Windowsの全試験

```sh
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo test --locked --release
cargo build --locked --release
python -m unittest discover -s tests -p 'test_*.py' -v
python -m unittest discover -s tools -p 'test_*.py' -v
```

Windows CIと同じ分割をPCで実行する場合は、各出力先を新しくし、次の全4群を省かず実行する。

```sh
python tools/run_ci_python_tests.py --group python-core-0 --output work/NEW-python-core-0
python tools/run_ci_python_tests.py --group python-core-1 --output work/NEW-python-core-1
python tools/run_ci_python_tests.py --group python-mean-evidence --output work/NEW-python-mean
python tools/run_ci_python_tests.py --group python-spatial-evidence --output work/NEW-python-spatial
```

各 `python-test-inventory.json` の全発見IDと全4群の和集合を照合する。どれか1群の成功だけで全件成功にしない。toolsの2試験も割り当てに含まれる。既定試験へ重い新規弧長求解や新アーカイブの再監査を追加していない。新規の拒否試験は小さい入力だけを使用する。

起点のWindows Pythonでは、隔離子プロセスのcp1252によるUTF-8読取失敗と、45分上限のキャンセルを確認した。旧再監査器の子プロセスには隔離を保持した `-I -X utf8` を明示した。今回の分割は全試験の維持を目的とし、許容値・assert・スキップ条件を緩めない。最終先端のCIは、実際のSHAに対応するrunと全ジョブ・集約の状態を確認する。queued/in_progress/未取得は成功未確認である。

旧証拠の修正後再監査を単独で実行する場合は次を使う。これは旧封印証拠の再現で、現在の生産コードの全撤去認定ではない。

```sh
python interim/2026-11-10-nc-mean-stress-release/remeasure.py --output work/NEW-historical-utf8-replay
```

Windowsでは既存の `--allow-platform-float-differences` 選択も使える。元の全SHA・離散構造の一致検査は保持し、派生浮動小数の差を列挙するもので、差を全解放合格へ変換しない。Linux/WSLでの厳密バイト比較とは区別する。

## JOINTとその他の残件

新監査は既存サイクルの保存履歴を既存写像で再生する整合性検査である。一般的な全履歴収束、独立した本体全場比較、一般塑性・部材・水の同時経路の完成証拠にはならない。

部分亀裂先端は利用者の形状・先端位置・共有節点の方針が必要。接触BIRTHは既存履歴と新しい設置基準・水理初期値の明示仕様が必要。板/地盤補間整合も高次板または弱結合の仕様・仮想仕事と独立3D証拠が必要。これらを既存の全面分離柱から捏造しない。

Richards湿潤・実2D互換JOINT、一般地質・実CAD・湛水/流出、一般塑性変位、杭列/群杭、梁塑性/座屈等には、それぞれ利用者入力、独立3Dの力/仕事/水収支、時空間分割、棄却と途中再開が必要。既存小ケースの限定結果を転用せず、証拠がないものは拒否または未検証として残す。P3b/P4の各行、構造・水・JOINT全体、P6全体は未完である。
