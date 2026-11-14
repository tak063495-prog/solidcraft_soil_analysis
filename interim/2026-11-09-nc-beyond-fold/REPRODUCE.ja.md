# ユーザーPCでの再計測・延長

対象は既存の小さい5 kPa低プレロードNCだけ。100 kPa追加プレロードの結果は使用しない。
このフォルダのプログラムは診断用で、本番jobに増分を確定せず、合格を出さない。

## 準備

PR #13の `codex/plastic-spatial-convergence` を取得する。起点は `046fd22ddced2e5b66b80be48c552508affb7053`。
本追補を含む先端を使う。mainだけの配置で実行しない。Rustは `rust-toolchain.toml`、依存は既存lockfileを使う。
Python 3.13と `requirements-reference.lock.txt` のnumpy/scipy等が必要。出力先は毎回未存在のフォルダにする。

PowerShellの例（リポジトリのルート）：

```powershell
git fetch origin codex/plastic-spatial-convergence
git switch codex/plastic-spatial-convergence
git merge-base --is-ancestor 046fd22ddced2e5b66b80be48c552508affb7053 HEAD
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements-reference.lock.txt
rustup show
cargo fetch --locked
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
```

既存venvがあれば再作成しない。Linuxでは既存venvを有効化して同じPython引数を使い、環境変数は `export` で設定する。
WindowsバイナリのSHAはLinuxと異なる。封印済みLinux本番jobのバイナリや封印を付け替えない。

## 求解しない再計測

```powershell
python interim/2026-11-09-nc-beyond-fold/remeasure.py --output work/nc-beyond-remeasured
python interim/2026-11-09-nc-beyond-fold/compare_history.py --evidence work/nc-beyond-remeasured/evidence --output work/nc-beyond-remeasured/full-path-sensitivity.json
```

分割アーカイブ・全エントリー・診断ソースのSHA256を検証し、全受理点の応力から力を再組立する。
保存された拡大系の真の残差、固有値、仕事、全場の刻み差と再開一致も再計測する。材料更新や釣合い求解は行わない。
`remeasurement.json` の `exact_saved_measurements=true` を確認する。異なる場合は両JSONを保存し、値を手で合わせない。
補足の全場刻み監査は `full-path-sensitivity.json` の両写像 `matched_points=1024` と保存版との一致を確認する。
この補足JSONは別の後処理で、アーカイブ内の主計測再生成とは分けている。新しい材料求解はしない。
ディスクは新証拠の展開に数GB、前回証拠も展開する再実行では約15 GB以上の空きを見込む。
再計測は同程度のCPUで数分を想定するが、ディスク速度にも依存する。

## 今回の有界な数値診断を再実行

```powershell
python interim/2026-11-09-nc-beyond-fold/reproduce.py --output work/nc-beyond-pc --list-only
python interim/2026-11-09-nc-beyond-fold/reproduce.py --output work/nc-beyond-pc
```

前回証拠4,513エントリーを照合・再計測し、同じ69ファイルの観測crateを再生成してビルドする。
Windowsでは `.exe` を選ぶ。前回 `prepare_evidence.py` のLinux固定バイナリ名の後処理は使わず、`--source-only` と別ビルドを使う。
native・独立の各写像について、粗刻み512点、細刻み1,024点、連続粗刻み1,024点、粗刻み途中再開512点、細刻み延長1,024点を実行する。
力1e-8、真の線形残差1e-10、平面算術残差1e-12 m、Newton20回、弧長単位換算1 mは前回から維持する。
各ドライバの予算は600秒。終了コード0でも停止理由を必ず読む。スクリプト内の許容値を変更して通さない。

同程度のCPUでは、ビルドを含め十数分～数十分が目安。PCでは10～60分程度を見込み、遅いPCでは各600秒予算で途中停止し得る。
今回の10数値起動の実測時間の合計は **1249.239765567秒**（並列実行を含むため壁時計時間ではない）。初回の隔離ビルドはログで2分10秒だった。
これは時間の保証ではない。実測時間は各 `report.json.elapsed_seconds` と `processes.json` に残る。
この一括ラッパー全体はクラウドで再実行していない。クラウドでは同じ10数値ドライバを個別に実行した。
前回の独立前史を再生成する場合は、先に前回フォルダの `reproduce.py` を別出力先で実行する。その所要時間は上記に含めない。

## さらに先へ延長する場合

前回終了点を改変せず、同じ刻み・写像・材料分割のまま再開する。以下は独立写像の例：

```powershell
python interim/2026-11-09-nc-beyond-fold/continue_path.py --map independent --reference work/nc-beyond-pc/prior-verified/evidence/independent-input-rk4-fine --arc work/nc-beyond-pc/prior-verified/evidence/independent-arc-history400 --resume work/nc-beyond-pc/extended-independent-whole/restart.json --output work/nc-beyond-pc/independent-extension --step-size 0.0001 --count 1024 --budget-seconds 600
```

nativeでは `--map native --evidence work/nc-beyond-pc/prior-verified/evidence/prior --binary work/nc-beyond-pc/prepared/target/release/examples/material_diagnostic.exe` とし、native側の保存点を指定する。
Linuxでは拡張子を除く。λは[0,1]に限定し、正方向に1を跨ぐ予測の場合のみλ=1の端点補正を試す。
単一実行は最大4,096点・1,800秒まで。境界外試行、時間切れ、未収束、線形検査超過は判定不能で止まる。
無期限ループにせず、各区間を保存・点検する。刻みを半減する比較は同じ起点から別経路として行う。
途中で刻み・材料分割・束縛対象ソース・バイナリを変えたcheckpointは拒否する。束縛するファイルは `controls.json.binding.source_sha256` に列挙される。

独立材料RK4の4→8分割を調べる場合は `--material-steps 8` で同じ独立起点から新規実行し、4分割のcheckpointを流用しない。
前史400→800の差も別に確認する。今回の弧長刻み感度だけで全履歴収束を認定しない。

## 確認するJSON

- 各 `report.json`：`stop_reason`、`equilibrium_samples`、`last_progress`、`full_removal_reached`、`elapsed_seconds`。
- `analysis/measurements.json`：`passed`、`verdict`、`static_limit_confirmed`、`history_converged`。
- `paths`：`maximum_reassembled_relative_force_residual`、`maximum_recomputed_true_linear_residual`、`linear_systems_above_1e_10`、`removed_material_frozen`、`contact_history_frozen`。
- `remaining_region_lambda_zero_*`、`step_sensitivity` の全節点変位・全GP応力/塑性ひずみ/pc/間隙比、`restarts.*.all_saved_fields_exact`。
- `positive_to_negative_eigen_crossings`、仕事の符号と両写像の差。対称化接線の負固有値だけで物理的破壊を認定しない。
- λ=1へ到達した場合も `full-release.json` の釣合い・ゼロ増分凍結、全解放後の履歴、独立全場、刻み感度、完了までの再開を確認する。診断コードは自動で `passed=true` にしない。

時間切れ・未収束・カットバック枯渇は破壊にも合格にもしない。下降後の再上昇が今回の有界区間にないことは、全域で枝が存在しない証明ではない。
静的限界を確定するには、両写像と前史を含む独立3D証拠・刻み感度・妥当な安定性の説明が必要。未探索区間を残したまま確定しない。

## 試験

```powershell
$env:NC_EXTENDED_MEASUREMENTS='interim/2026-11-09-nc-beyond-fold/measurements.json'
python -m unittest discover -s interim/2026-11-09-nc-beyond-fold -p 'test_*.py' -v
$env:NC_ARC_MEASUREMENTS='interim/2026-11-09-nc-pseudo-arclength/measurements.json'
python -m unittest discover -s interim/2026-11-09-nc-pseudo-arclength -p test_diagnostic.py -v
cargo test --locked --offline --release --test mcc_midpoint --test free_contact
```

新しいスカラー多項式の試験は弧長アルゴリズム・端点処理・checkpoint拒否のソフトウェア試験で、地盤の新形状・材料仕様・参照解ではない。
Windowsでこれらを実行した結果は今回未取得。既存Windows CIは `tests/` と `tools/` を探索し、このinterimの新規12 Python試験は自動対象に含まない。
旧先端CIやLinux成功を今回のWindows成功へ転用しない。
