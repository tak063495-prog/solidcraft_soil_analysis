# 再現・保存点継続・ユーザーPC手順

PR #13はマージしない。起点は `1c0d57b6943efbe0df85659f99ba47a7ffa15205`、親は `1b5fcd86542be615dafcc8dc9b27efca474006b9`。同じブランチを使い、既存の未保存変更を消さない。

```sh
git fetch origin codex/plastic-spatial-convergence
git switch codex/plastic-spatial-convergence
git pull --ff-only origin codex/plastic-spatial-convergence
git merge-base --is-ancestor 1c0d57b6943efbe0df85659f99ba47a7ffa15205 HEAD
```

今回のLinux環境はRust 1.99.1、Python 3.13.17、NumPy 2.5.4、SciPy 1.18.1。既存の `requirements-reference.lock.txt` を使う。旧証拠・Git pack・一時展開にも容量が必要である。今回の全Python試験は、疎取得から省いた旧証拠を戻してから成功を確認した。必要なファイルを欠落させたまま試験をスキップしない。

```sh
export CARGO_BUILD_JOBS=2 RUST_TEST_THREADS=2
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
export PYTHONUTF8=1 PYTHONIOENCODING=utf-8
```

以下の出力先は未存在のものを選ぶ。封印したLinuxの証拠の厳密なバイト比較には、同じ依存関係のLinux/WSLを使う。Windows/別BLASで浮動小数の最下位桁が変われば、厳密再計測は失敗として残す。しきい値を緩めたり、保存点のSHA・実行版封印を付け替えたりしない。

## 新アーカイブの全件監査（新しい釣合い求解なし）

```sh
python interim/2026-11-10-nc-oblique-and-joint-history/remeasure.py --output work/NEW-oblique-replay
```

全分割・結合全体・全エントリーのSHA/サイズ/安全なパスを検査する。封印したソースを `-I -X utf8` の別プロセスで使い、NCの保存場再組立JSONの完全一致、JOINTの全受理点の材料/接触/水圧/反力監査と粗密比較の一致を要求する。`verification.json` の各監査フラグを確認する。`numerical_acceptance=false` は維持される。これは新しい釣合い解・破壊判定・全撤去合格ではない。

## NCの新規比較

対象は既存 `make_free_contact.py --mcc --pc 100 --preload-kpa 5 --force 1 .5 0` だけ。pcの100 kPaと追加プレロード5 kPaを区別する。

まず旧bounded-refinementの全122分割を配置する。疎取得の場合も、SHA検査を省かない。

```sh
git sparse-checkout add '/interim/2026-11-10-nc-bounded-refinement/*.part*'
cargo build --locked --release
python interim/2026-11-10-nc-tangent-branch-diagnosis/extract.py --output work/NEW-oblique-prior
python interim/2026-11-10-nc-tangent-branch-diagnosis/prepare_observer.py --output work/NEW-oblique-observer
python interim/2026-11-10-nc-mean-stress-release/fresh_prefix.py --solver target/release/soil-fem-3d --evidence work/NEW-oblique-prior/prepared/origin --output work/NEW-oblique-prefix
```

extractは旧アーカイブ2,034,108,934バイト・全1,349エントリーを照合し、必要な保存場だけを残す。観測器は現在のRustソースの隔離コピーへ、既存の読み取り専用材料観測関数を加える。旧ヘルパーの固定コミットラベルは今回の呼出しコミットではない。元RustソースSHA・観測器SHA・新しいprefixの全状態一致を別に確認する。生産Rustをパッチしない。

```sh
python interim/2026-11-10-nc-oblique-and-joint-history/reproduce.py --evidence work/NEW-oblique-prior --binary work/NEW-oblique-observer/target/release/examples/material_diagnostic --output work/NEW-oblique-study --invocation-commit "$(git rev-parse HEAD)" --budget-seconds 180
python interim/2026-11-10-nc-oblique-and-joint-history/extend.py --study work/NEW-oblique-study --evidence work/NEW-oblique-prior --binary work/NEW-oblique-observer/target/release/examples/material_diagnostic --budget-seconds 120
python interim/2026-11-10-nc-oblique-and-joint-history/measure_nc.py --source work/NEW-oblique-study/paths --output work/NEW-oblique-measurements.json
```

reproduceは最大2子プロセスで、本体/独立、2折線経路、初期刻み1/128・1/256の8比較を行い、実midpointからの部分再開を照合する。extendは各未完経路の実terminal checkpointをそのまま別プロセスへ渡す。時間予算とPC速度で終端は変わる。処理の終了コード0と数値受入を区別する。

独立側は入力のみから計算済みの自身の400分割前史をSHA照合して使い、新しい撤去履歴を自身のRK4材料写像で積む。本体の応力・塑性履歴・Hを独立解へ代入しない。ただし400分割前史を今回の新規求解と呼ばず、全撤去の独立全場/刻み収束も未達のままとする。

同じ既存の符号保存mean/complement分割H=A+Bに対し、s=0.5で解放率(a,b)=(0.8,0.2)または(0.2,0.8)、両端(0,0)→(1,1)を結ぶ2本の折線を使う。第1群Aはcomplement、第2群Bはclip済みmeanである。荷重は `Ffinal+(1-a)A+(1-b)B`。各元H成分の符号・上限を全評価で検査し、追加拘束圧・逆符号荷重・元ゼロ成分への追加を拒否する。sは元の比例撤去λではなく、s=1かつH実ゼロだけが全解放に対応する。生産施工手順として採用していない。

力1e-8、真の線形残差1e-10、根補正1e-11 m、材料増分ひずみL1上限0.001、Newton20回、既存線探索・反復改良を維持する。最小刻みは2^-22。時間・点数・最小刻みの停止を静的限界や破壊へ読み替えない。

## JOINTの既存小ケース

```sh
python interim/2026-11-10-nc-oblique-and-joint-history/joint_study.py --output work/NEW-joint-bulk-history
python tools/audit_wet_joint_bulk_history.py --history work/NEW-joint-bulk-history/g32-l8/accepted --result work/NEW-joint-bulk-history/g32-l8/result.json --output work/NEW-joint-bulk-audit.json
```

既存54節点・12 Tet10・324材料GP・32接触点の同じCONTINUOUS湿潤サイクルだけ。16/32全体分割、接触内8分割、既存解析水圧体積積分・GL7、全7工程を用いる。新しい形状・亀裂先端・構成則・接触BIRTH・高次板/弱結合を作らない。入力から得た独立弾性MC枝＋開離/すべり参照の保存場整合監査であり、生産側の独立全履歴比較、一般塑性、部材/水の全体完了ではない。

## WindowsユーザーPCでの全件試験と新しい生産job

PowerShellでは上記各環境変数を `$env:変数名='値'` で設定する。Pythonは `.venv/Scripts/python.exe`、実行版は `.exe` に読み替える。

```powershell
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release
cargo test --locked --release
.venv/Scripts/python.exe tools/run_ci_python_tests.py --group python-core-0 --output work/NEW-py-core0
.venv/Scripts/python.exe tools/run_ci_python_tests.py --group python-core-1 --output work/NEW-py-core1
.venv/Scripts/python.exe tools/run_ci_python_tests.py --group python-mean-evidence --output work/NEW-py-mean
.venv/Scripts/python.exe tools/run_ci_python_tests.py --group python-spatial-evidence --output work/NEW-py-spatial
.venv/Scripts/python.exe tools/probe_nc_midpoint.py --solver target/release/soil-fem-3d.exe --output work/NEW-windows-nc --scheme midpoint --budget-seconds 300
# 同一実行版・同一入力の、自分のWindows jobだけを再開する。
.venv/Scripts/python.exe tools/probe_nc_midpoint.py --solver target/release/soil-fem-3d.exe --output work/NEW-windows-nc-resumed --scheme midpoint --budget-seconds 300 --resume-from work/NEW-windows-nc/job
```

これは元の生産側比例撤去の別検証であり、今回の折線診断の合格ではない。Linux保存点をWindows実行版へ再封印しない。Windowsの新規求解・JOINT監査は今回の先端CIやユーザーPCで実結果を確認するまで成功と記載しない。

## 未完了の受入項目

NCは撤去中の全材料/接触履歴、非零推定からの残存域λ=0実Newton釣合い、全解放H=0とその後の実ゼロ増分凍結、入力だけからの独立塑性3D完了全場、両実装の完了までの刻み感度、実checkpointから完了までの再開を揃える。部分場・部分再開・100 kPa追加プレロードの別合格で代用しない。

部分亀裂先端には利用者の形状/先端/共有節点方針、接触BIRTHには履歴・設置参照・水理初期値の仕様、高次板/弱結合には補間/仮想仕事の仕様と独立3D証拠が必要。今回それらを推測せず、証拠不足は拒否またはユーザーPCの残件とする。P3b/P4各行、構造・水・JOINT全体、P6全体は未完である。
