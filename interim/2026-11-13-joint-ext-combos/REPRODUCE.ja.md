# ユーザーPCでの再現

PR #13の `codex/plastic-spatial-convergence` を取得する。マージは不要である。Rust 1.99.1、Python 3.13と固定依存を使う。Windowsでは短い作業パスを推奨する。以下はPowerShellの例で、各コマンドの終了コードが0であることを確認する。

```powershell
git switch codex/plastic-spatial-convergence
git pull --ff-only
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements-reference.lock.txt
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:CARGO_BUILD_JOBS='2'
$env:RUST_TEST_THREADS='2'
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo test --locked --release
cargo build --locked --release --examples
python -m unittest discover -s tests -p 'test_*.py' -v
python -m unittest discover -s tools -p 'test_*.py' -v
python tools/verify_joint_ext_combos.py --solver target/release/soil-fem-3d.exe --rehash target/release/examples/joint_ext_checkpoint_rehash.exe --output runs/joint-ext-user-001
```

続けて途中状態の独立再組立も実行する。

```powershell
python interim/2026-11-13-joint-ext-combos/audit_checkpoints.py --evidence runs/joint-ext-user-001 --output runs/joint-ext-user-001/checkpoint-reassembly.json
```

出力先は毎回未作成のパスを使う。Linuxでは実行ファイルの `.exe` を除く。ドライバが実行版を出力先へ固定コピーし、別CLI再開・同一受理刻み列・再ハッシュ改変拒否を検査する。`measurements.json` の `passed` を確認し、途中停止・失敗は保存したまま別名でやり直す。最終Linux全件の件数は `validation.json` を参照し、Windows固有試験の件数をLinuxと混同しない。

## 小型入力を個別に調べる

同梱 `evidence.zip` の各ケースに `mesh.json`、`ground.json`、`case.json`、`run.json` がある。`restart-*.zip` には再開・改変拒否の元データがある。別の新規出力先へ `solve-stages` で再計算する。`checkpoint-stages MESH GROUND CASE 受理増分数 NEW_CHECKPOINT`、`resume-stages MESH GROUND CASE CHECKPOINT NEW_RUN` で再開する。初期平衡の受理増分も数える。測定JSONの保存地点一覧を使う。

入力の `hydraulic_initialization` と `states` を省略して湿潤BIRTHを試した場合は拒否が正しい。静水位、初期有効応力、材料履歴、配置を利用者が指定する。既存の定常流を初期水圧の代用にせず、今回未受入のSOIL／CONTINUOUS／PRESCRIBEDの流れ付きBIRTHとRichards湿潤JOINTを有効にするためガードを外さない。

## 大規模・実CAD・GUI/GPU（この受入では未実行）

実CADでは、元ボディ・材料・面集合・先端辺3節点・両所有要素を保存し、案1の明示先端分離手順でTet10メッシュを作る。未宣言の節点スターを切らず、共有先端の係数相殺と体積・面の向き・対応の一意性を検査する。先端の空間細分化、材料非線形、部分水没、構造部材との複合を小型結果から合格としない。

先端Aは2000GP／6000物理節点成分の厳密投影上限を維持する。大きい場合はまず小型部分モデルで入力と独立釣合いを確認し、上限を外す変更は別の検証に分ける。時間刻み・空間刻みを個別に細分化し、非収束を含む全試行を保存する。却下済みNCチューニングを再実行しない。

VTU/PVD/CSVは `tools/export_joint_birth.py --help` と案2の手順で書き出し、実際のSolidCraft/ParaView上で、地盤・先端・接合世代・初期力・誕生基準・水圧・凍結履歴の表示値を生JSONと照合する。GPU描画・大規模実行の成功を、データ往復だけで認定しない。使用機器・版・実行時間・全入力・ログ・スクリーンショットを別の証拠パケットへ保存する。

Richards結合を将来検証する場合はONE_WAYか双方向かを最初に仕様化し、変位節点と水頭節点の写像、開離／再接触と接合誕生／撤去の時刻、独立の水量収支と接合面水圧の荷重再組立、浸潤・排水の刻み細分化を実装・試験する必要がある。現行の拒否は、この手順を書いただけでは解除しない。
