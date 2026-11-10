# 再現手順（ユーザーPC）

PR #13の `codex/plastic-spatial-convergence` を取得する。マージ操作は不要。Rust 1.99.1とPython 3.13、既存の `requirements-reference.lock.txt` を使う。以下はリポジトリルートで実行し、出力先には未作成のディレクトリを指定する。PowerShellは失敗したコマンドの `$LASTEXITCODE` が0であることを各段で確認する。

```powershell
git fetch origin codex/plastic-spatial-convergence
git switch codex/plastic-spatial-convergence
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
python -m venv .venv
.venv\Scripts\python -m pip install -r requirements-reference.lock.txt
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo test --locked --release
cargo build --locked --release
.venv\Scripts\python -m unittest discover -s tests -p 'test_*.py' -v
.venv\Scripts\python -m unittest discover -s tools -p 'test_*.py' -v
.venv\Scripts\python tools/plate6_study.py runs/plate6-user-pc
```

LinuxではPythonを `.venv/bin/python` に読み替える。BLASの並列数を1にする場合は `OPENBLAS_NUM_THREADS=1`、`MKL_NUM_THREADS=1`、`OMP_NUM_THREADS=1` を設定する。小規模studyは同じ実行バイナリをコピーして、5地点の別CLI再開を実施する。`acceptance.json` の `passed: true` を確認し、数値残差・独立比較・出力・改変拒否の全ゲートを検査する。時間切れ、CLI非ゼロ終了、未完了、未収束を合格にしない。

## 入力

コンパクトな面宣言の各surfaceに `"formulation":"mindlin_tri6"` を加える。対象地盤面は直線・平面の6節点面である必要がある。直接部材入力の場合は `nodes:[頂点0,頂点1,頂点2,中間01,中間12,中間20]` と同じformulationを指定する。省略時と `cst_dkt` は従来方式。`thickness_offset_m` は非ゼロを指定しない。荷重は面tractionで指定し、6節点均等配分を外部で作らない。

## 大規模・実CAD・GUI/GPU（クラウドでは未実行）

1. 既存の有限幅壁/掘削モデルを別ジョブへ複製し、旧入力・ソース指紋・4.740/3.466 μmを含む旧結果を保存する。面宣言の方式だけを明示的に切り替え、断面・厚さ・荷重・参照面・土のメッシュ/構成則・全許容値を保持する。JOINT、水、土BIRTHの既存の拒否を回避しない。
2. `tools/prepare_native_plates.py SOURCE NEW_OUTPUT` で展開する。`plate-surface-binding.json` の節点順、面の所有者、面積、板種別、元地盤の座標/接続不変を確認する。曲面等を検査のために平坦化しない。
3. 従来のジョブ管理手順でLinux/Windows CLIを実行する。構造容量・記憶域は既存の `storage_budget` で明示し、超過を切り捨てない。元のNewton/線形/再組立許容値を保持する。粗密・厚さ・drilling感度は別系列として記録する。
4. 受理点の保存・同じバイナリでの別プロセス再開、棄却時の確定状態保持、設置基準、せん断を含む全エネルギー、全反力・モーメント、VTU再読込を検査する。`plate6_study.py` の小規模直接解を大規模モデルの独立解の代わりにしない。
5. SolidCraftの既存の表示/出力経路へ渡し、triangle6を6節点の二次補間として表示できることを確認する。3頂点だけで描画する旧表示では受入にしない。全工程の節点・辺内・面内の変位、生の曲率/せん断、局所軸、厚さ、断面、撤去時消去、CAD失効/Undoを公式GPUで確認する。現時点でこの新方式の実Native表示は未検証であり、表示側が未対応なら残件として止める。

大規模系列、実CADでの界面伝達/空間収束、別PC・GPU、Windows CIの成功は今回のLinux小規模結果から推定しない。P6一覧の未達をそのまま保持する。
