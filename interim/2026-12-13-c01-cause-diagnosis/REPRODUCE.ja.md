# C-01 原因診断の再現と PC 手順

PR #13、`codex/plastic-spatial-convergence`、起点 `e8133f045a852c6ca374a552070526067f11dbe9`。
マージしない。以下はリポジトリルートの PowerShell。NEW を未存在の出力名へ変える。
PC 計算は今回実行していない。時刻の実単位は定義されない境界係数サイクルなので、履歴比較は工程ごとの進行度に対応させる。

## 準備

```powershell
git fetch origin codex/plastic-spatial-convergence
git switch codex/plastic-spatial-convergence
git pull --ff-only origin codex/plastic-spatial-convergence
git merge-base --is-ancestor e8133f045a852c6ca374a552070526067f11dbe9 HEAD
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
py -3.13 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements-reference.lock.txt
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:CARGO_BUILD_JOBS='2'
$env:RUST_TEST_THREADS='2'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
cargo build --locked --release --example wet_joint_diagnosis_probe
```

Linux/WSL は `.venv/bin/python`、実行版は拡張子なし。数値停止・時間切れ・上限到達は失敗のまま残す。
証拠の SHA は実バイトで確認し、JSON/JSONL/ログ/part を改行変換しない。文書比較だけは `read_text` の改行正規化で CRLF も扱う。

## 保存証拠だけの再計算（新規平衡求解はしない）

```powershell
python interim/2026-12-13-c01-cause-diagnosis/c01_cause_evidence.py --output work/NEW-c01-cause-replay
```

新旧 part・ZIP と選択エントリーの SHA を検査し、各 GP のエネルギー・散逸と全履歴 L2、GP 変更比較、参照 256/512 比較を再計算する。
出力は `diagnostic/`、`saved-switch/paths/`、`saved-path/paths/`、`comparison/comparison.json`。
数分〜15分、RAM 8〜16 GB、ディスク 6〜12 GB の目安。Windows 実測ではない。
元の非重み付き L2、面積重み付き L2、進行度台形積分 L2、符号付きの端点/台形散逸診断を混同しない。
新しい定義が単調でも合格へ変更しない。

## 本体の内部接触分割数を独立に採取

```powershell
python interim/2026-12-13-c01-cause-diagnosis/c01_cause_diagnostics.py --extract switch --output work/NEW-c01-switch-saved
python interim/2026-12-13-c01-cause-diagnosis/c01_cause_diagnostics.py --extract path --output work/NEW-c01-path-saved
foreach ($n in 16,32,64,128,256) {
  python interim/2026-12-13-c01-cause-diagnosis/c01_cause_capture_native.py --saved work/NEW-c01-switch-saved/paths --output work/NEW-c01-counts-switch/n$n --solver target/release/examples/wet_joint_diagnosis_probe.exe --control switch --steps $n
  python interim/2026-12-13-c01-cause-diagnosis/c01_cause_capture_native.py --saved work/NEW-c01-path-saved/paths --output work/NEW-c01-counts-path/n$n --solver target/release/examples/wet_joint_diagnosis_probe.exe --control path --steps $n
}
python interim/2026-12-13-c01-cause-diagnosis/c01_cause_diagnostics.py --input work/NEW-c01-switch-saved/paths --native-counts work/NEW-c01-counts-switch --output work/NEW-c01-post-switch
python interim/2026-12-13-c01-cause-diagnosis/c01_cause_diagnostics.py --input work/NEW-c01-path-saved/paths --native-counts work/NEW-c01-counts-path --output work/NEW-c01-post-path
```

全7工程、54節点・12 Tet10・32接触GPの固定ケース。観測専用 example は既存 Observer のコピーだけを受け取り、停止・分割・Newton・解を変更しない。
`raw/trials.jsonl` の収束試算を、保存済み受理状態の旧/新変位の全バイト・工程・進行度で一意に照合する。
一致しなければ推定分割数で代用せず停止する。`counts.json` と `measurement.json` が出力。
10経路で想定5〜30分、RAM16 GB、追加ディスク2〜6 GB。Windows固有差を許容値で吸収しない。
`--native-counts` を省いた再計算は受理端点からの分割数下限という別診断で、実際の内部履歴の再現とはしない。

## 小規模の GP 比較（明示した独立参照だけを変更）

```powershell
foreach ($n in 32,64) {
  foreach ($q in 4,5,6) {
    python interim/2026-12-13-c01-cause-diagnosis/c01_cause_experiment.py --control switch --steps $n --contact-order $q --output work/NEW-c01-experiments/switch-n$n-q$q
  }
}
python interim/2026-12-13-c01-cause-diagnosis/c01_cause_experiment.py --control switch --steps 512 --contact-order 4 --output work/NEW-c01-experiments/switch-n512-q4
python interim/2026-12-13-c01-cause-diagnosis/c01_cause_compare.py --saved-switch work/NEW-c01-switch-saved/paths --saved-path work/NEW-c01-path-saved/paths --experiments work/NEW-c01-experiments --output work/NEW-c01-comparison
```

Jacobi 4×4/5×5/6×6 ×2面、32/50/72接触GP。`--contact-order` の既定4は既存独立参照をそのまま使う。
本体入力・水荷重・体積積分・制御は変更しない。GPが異なる履歴を補間せず、共通節点変位と面積積分した量だけ比較する。
各経路は `independent/*.npz`、`result.json`、`integrals.json`、`geometry.json`、`summary.json`。
上記7経路で想定2〜15分、RAM16 GB、追加ディスク1〜3 GB。見る値は全7工程の完了、元の力残差1e-10、積分面積・二次patch、同じ4×4と旧保存解の差0、面積積分の差、共通節点変位差。
これらの積分量は元の全GP全履歴 L2 の代替受入ではない。

## PCのみ：参照自体の長い系列

```powershell
foreach ($n in 512,1024) {
  python interim/2026-12-13-c01-cause-diagnosis/c01_cause_experiment.py --control switch --steps $n --contact-order 4 --pc-series --output work/NEW-PC-c01-reference/switch-n$n-q4
}
python interim/2026-12-13-c01-cause-diagnosis/c01_cause_compare.py --saved-switch work/NEW-c01-switch-saved/paths --saved-path work/NEW-c01-path-saved/paths --experiments work/NEW-PC-c01-reference --reference-levels 256 512 1024 --output work/NEW-PC-c01-reference-comparison1
python interim/2026-12-13-c01-cause-diagnosis/c01_cause_experiment.py --control path --steps 512 --contact-order 4 --output work/NEW-PC-c01-reference/path-n512-q4
```

切替制御の参照512/1024は想定5〜60分以上、RAM16 GB、ディスク5〜15 GB。実測の保証ではない。
誤差制御512はクラウドで既存2000試行/工程の上限に達した同じ設定を再現する失敗手順。
上限・許容値・誤差scaleは緩和しない。`failure.json`、途中 `independent/manifest.json` と全履歴を保存し、完了扱いしない。
長系列は各 `result.json` の `path` から同一工程・実在する `i/256`、`i/512` の座標だけを取り、
256/512/1024の同一座標で元の非重み付き全GP L2・絶対差・方向余弦を比較する。
1024も既存2000試行/工程以内で完了した場合だけ比較する。2048は基準受理数だけでこの上限を超えるため候補にしない。
元の全場比較が必要な場合は既存の汎用 PC study/measure を用いる（次項）。

## PCのみ：本体・独立の全場、接触内上限、空間系列

```powershell
cargo build --locked --release --example wet_joint_path_probe --example wet_joint_history_probe
python interim/2026-12-12-c01-switch-event-localization/study.py --mode init --solver target/release/examples/wet_joint_path_probe.exe --output work/NEW-PC-c01-full
python interim/2026-12-12-c01-switch-event-localization/study.py --mode native --pc-series --steps 256 512 1024 --timeout 14400 --output work/NEW-PC-c01-full
python interim/2026-12-12-c01-switch-event-localization/study.py --mode independent --pc-series --steps 256 512 1024 --output work/NEW-PC-c01-full
python interim/2026-12-12-c01-switch-event-localization/measure.py work/NEW-PC-c01-full work/NEW-PC-c01-full-measure --steps 256 512 1024
python interim/2026-11-15-wet-joint-convergence/user_pc.py --solver target/release/examples/wet_joint_history_probe.exe --output work/NEW-PC-c01-local --levels 0 --steps 256 512 1024 --maximum 0.000025 --native-timeout-seconds 14400
python interim/2026-11-18-wet-joint-step-convergence/pc_summary.py --input work/NEW-PC-c01-local --output work/NEW-PC-c01-local-summary --steps 256 512 1024 --h 0 --maximum 0.000025
python interim/2026-11-15-wet-joint-convergence/user_pc.py --solver target/release/examples/wet_joint_history_probe.exe --output work/NEW-PC-c02-space --levels 0 1 2 --steps 64 128 256 --maximum 0.000025 --native-timeout-seconds 14400
```

前4行の全場系列は想定1〜6時間以上・30〜60 GB・RAM16 GB。停止は停止のまま残す。
接触内上限の交差対照は想定30分〜2時間以上・30 GB・RAM16 GB。
12/96/768 Tet10の空間系列は想定数時間〜1日・100 GB・RAM64 GB以上。
出力は `gN/` または `hH-gN-m2.5e-05-cut/`、別集計先の `measurements.json`。
見る値は全18場の元の全履歴 L2・方向・次数、元の力/線形残差、GPカテゴリ、全工程完了。
共通物理点・接触方式自体の空間検証は別途必要で、C-02〜C-04、NC-01、P6、A/Bを完了にしない。

## 全件検査と固定台帳

```powershell
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release
cargo test --locked --release
python tools/run_ci_python_tests.py --group python-core-0 --output work/NEW-c01-check/core0
python tools/run_ci_python_tests.py --group python-core-1 --output work/NEW-c01-check/core1
python tools/run_ci_python_tests.py --group python-mean-evidence --output work/NEW-c01-check/mean
python tools/run_ci_python_tests.py --group python-spatial-evidence --output work/NEW-c01-check/spatial
```

PowerShell実物が必要。4群の発見IDと実行IDの重複なし和集合を照合し、成功/失敗/未実行/skipを `validation.json` に記録する。
Linuxでの成功はpush後のWindows CI成功ではない。新先端Windows CIは未確定。
新しい状態証拠JSONを先に不変コミットへ保存し、
`interim/2026-12-10-ci-p601-inventory-fix/register_sources.py --commit FULL_SHA --id UNIQUE_ID interim/2026-12-13-c01-cause-diagnosis/FILE.json`
の出力を既存 `inventory.toml` に追記する。既存登録は再生成しない。`SUPPLEMENTAL_SHA256` を更新し、P6-01/P6-03を再検査する。
固定コミット・SHA-256・Git blob・全boolean/null pointerの照合を保持し、SOURCE_NOT_INVENTORIEDの例外やskipは追加しない。

C-01は未達、passed=false。
