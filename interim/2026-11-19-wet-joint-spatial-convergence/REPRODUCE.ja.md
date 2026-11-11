# C-02 の再現と Windows PC 手順

保存先ラベルは 2026-11-19、実施開始日は 2026-11-11（日本時間）。起点は `2d2d54c48070581d088959b52984ceb348e0e85e`、PR #13 の `codex/plastic-spatial-convergence`。PR はマージしない。各 `NEW-...` は未存在の出力先に置き換える。旧証拠・checkpoint の hash/snapshot・失敗は上書きしない。

## 準備

PowerShell、Python 3.13、Rust 1.99.1。終了コードが 0 以外なら停止してログ・途中場を保存する。以下の時間と資源は Windows での実測ではなく計画値であり、時間切れを成功にしない。

```powershell
git fetch origin codex/plastic-spatial-convergence
git switch codex/plastic-spatial-convergence
git pull --ff-only origin codex/plastic-spatial-convergence
git merge-base --is-ancestor 2d2d54c48070581d088959b52984ceb348e0e85e HEAD
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
cargo build --locked --release
cargo build --locked --release --example wet_joint_spatial_probe
$Study='interim/2026-11-19-wet-joint-spatial-convergence'
$Run='work/NEW-c02-solve'
python "$Study/study_c02.py" --mode init --solver target/release/examples/wet_joint_spatial_probe.exe --output $Run
```

Linux/WSL は `source /workspace/soil-environment/activate.sh`（このクラウド固有）または通常の仮想環境を使い、`.exe` を外す。ソルバー実行版は開始時にコピーして SHA を固定する。旧 checkpoint の実行版 SHA を新ビルドへ付け替えない。

## 小～中規模の全履歴新規求解

```powershell
foreach ($h in 0,1) {
  foreach ($g in 16,32,64,128,256) {
    python "$Study/study_c02.py" --mode native --h $h --steps $g --timeout 14400 --output $Run
    if ($LASTEXITCODE -ne 0) { throw "native failed: h=$h g=$g" }
    python "$Study/study_c02.py" --mode independent --h $h --steps $g --output $Run
    if ($LASTEXITCODE -ne 0) { throw "independent failed: h=$h g=$g" }
  }
}
foreach ($g in 32,64) {
  python "$Study/study_c02.py" --mode native --h 2 --steps $g --timeout 14400 --output $Run
  if ($LASTEXITCODE -ne 0) { throw "native failed: h=2 g=$g" }
}
python "$Study/study_c02.py" --mode restart --h 2 --steps 64 --timeout 14400 --output $Run
foreach ($q in 6,8) {
  python "$Study/ablation_c02.py" --mode solve --order $q --steps 256 --output $Run
  if ($LASTEXITCODE -ne 0) { throw "quadrature diagnostic failed: q=$q" }
}
$Measure='work/NEW-c02-measure'
python "$Study/measure_c02.py" $Run $Measure
python "$Study/ablation_c02.py" --mode measure --steps 256 --output $Run --target $Measure
python "$Study/ablation_c02.py" --mode water --output $Run --target $Measure
python "$Study/audit_raw_c02.py" $Run $Measure
```

`h=0/1/2` は同じ既存メッシュを一様 8 分割した 12/96/768 Tet10、54/250/1458 節点、324/2592/20736 材料 GP、32/128/512 接触 GP。各経路は初期平衡・shear・open・recontact・reverse・reload・return の全前史を入力から解き直す。既存の力 `1e-8`、真の線形残差 `1e-10`、独立系 `1e-10`、接触増分上限 `5e-5 m`、affine-cut 体積水圧を保持する。

主出力は `$Run/mesh-hH.json`、`hH-gG/native-raw/{accepted.jsonl,progress.jsonl,run.json,checkpoint-*.json}`、`hH-gG/{native,independent}/{manifest.json,*.npz}`、元の GP enum、全受理場監査・所要時間・入力と実行版 SHA。再開は最細水準 return 工程途中の**実ファイル**を別 PID が読み、全後続確定状態と Run 全体を時間だけ除いて比較する。全工程・一般塑性の再開受入を意味しない。

計測先は `measurements.json`、`common-coordinates.json`、全比較点寄与の `*-entities.npz` と SHA、`raw-material-audit.json`、`switch-brackets.json`、`contact-quadrature-diagnostic.json`、`water-affine-cut-diagnostic.json`。異なる受理時刻は補間しない。空間の比較点・補間方法・補間誤差の限界は README の定義に従う。

この一式は RAM 16～32 GB、ディスク 20～40 GB、1～3 時間を計画する。旧証拠を含む全件試験には別途 20～40 GB を確保する。同時多重実行は CPU を共有するため所要時間を増やす。クラウド実測は `timings.json` に記録し、Windows の保証時間にはしない。

## PC 専用：細かい固定時間・独立最細水準

クラウドで未測定の 768 要素・128/256/512 全体刻みと、最細の密行列独立系を追加する。最初は 128/256 のみ、資源を確認して 512 へ延長する。

```powershell
foreach ($g in 128,256,512) {
  python "$Study/study_c02.py" --mode native --h 2 --steps $g --timeout 28800 --output $Run --user-pc
  if ($LASTEXITCODE -ne 0) { throw "native incomplete: g=$g" }
  python "$Study/study_c02.py" --mode independent --h 2 --steps $g --output $Run --user-pc
  if ($LASTEXITCODE -ne 0) { throw "reference incomplete: g=$g" }
}
python "$Study/study_c02.py" --mode restart --h 2 --steps 512 --timeout 28800 --output $Run --user-pc
python "$Study/measure_c02.py" $Run work/NEW-c02-fine-time-measure
```

RAM 64 GB、ディスク 100 GB 以上、数時間～数日を計画する。独立系は既存の dense 3D 連立系なので最細水準では特に高価である。入力の `i/N` 全受理列・全 GP/全節点の全場・元の残差基準・実 checkpoint の別プロセス一致が必須。128/256 だけで「十分細かい」と認定せず、空間差の時間依存と負次数・方向反転を報告する。

## PC 専用：さらに細かい空間水準

6144 要素はクラウドで実行しない。新しい出力先に 4 水準を生成し、まず固定 64、次に固定 256、必要なら 512 と延長する。中間水準を欠かした 2 水準の次数は作らない。

```powershell
$Fine='work/NEW-PC-c02-h0123'
python "$Study/study_c02.py" --mode init --max-level 3 --user-pc --solver target/release/examples/wet_joint_spatial_probe.exe --output $Fine
foreach ($h in 0,1,2,3) {
  foreach ($g in 64,256) {
    python "$Study/study_c02.py" --mode native --h $h --steps $g --timeout 86400 --user-pc --output $Fine
    if ($LASTEXITCODE -ne 0) { throw "fine spatial route incomplete" }
  }
}
# 同時細分化 96/32 → 768/64 → 6144/128（不足経路だけ追加）
python "$Study/study_c02.py" --mode native --h 1 --steps 32 --timeout 86400 --user-pc --output $Fine
python "$Study/study_c02.py" --mode native --h 3 --steps 128 --timeout 86400 --user-pc --output $Fine
python "$Study/study_c02.py" --mode restart --h 3 --steps 128 --timeout 86400 --user-pc --output $Fine
python "$Study/measure_c02.py" $Fine work/NEW-PC-c02-h0123-measure --max-level 3
```

RAM 64～128 GB、空きディスク 500 GB 以上、1～数日を計画する（未実測）。大きい独立 dense 系を無条件に起動せず、768 要素の所要時間から予算を見積もる。GUI/GPU・実 CAD はこの数値ドライバの検証範囲外。

## 合否と全件ソフトウェア検証

整合確認は、全 7 工程・全受理列と全場の欠落なし、有限値、元の力・線形基準、接触摩擦円・仕事恒等式・履歴保持、弾性材料場の適合、独立参照、実 checkpoint の完全一致を要求する。不足・未収束・時間切れは失敗として保存する。

C-02 の数値判定は別である。同一物理座標・同じ時刻集合の絶対差ノルムから次数と方向を出す。非単調・負次数・補間不確かさ・開閉境界の不一致が残れば「未達」とする。正次数・面平均・7 端点ゲート・整合監査・CI 成功から新しい許容値を作らず、外挿場や Richardson 補正を出さない。**今回 C-02 は未達**。C-01・C-03・C-04、NC-01、P6 も未達のまま。

```powershell
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release
cargo test --locked --release
python tools/run_ci_python_tests.py --group python-core-0 --output work/NEW-c02-check/core0
python tools/run_ci_python_tests.py --group python-core-1 --output work/NEW-c02-check/core1
python tools/run_ci_python_tests.py --group python-mean-evidence --output work/NEW-c02-check/mean
python tools/run_ci_python_tests.py --group python-spatial-evidence --output work/NEW-c02-check/spatial
```

全発見 ID と 4 群の和集合が重複・欠落なく一致し、fail/error/skipped/ignored/filtered-out/expected-failure/unexpected-success が 0 であることを数える。sparse checkout は既存 CI と同じ `git sparse-checkout set --no-cone '/*' '!/interim/2026-11-10-nc-postfold-branches/evidence.tar.gz.part*'` にする。これは未使用保存系列の取得除外であり、試験の skip ではない。全件試験に必要な証拠の取得が完了するまで待つ。push 後の Windows 成否は今回の新先端の run で別途確認し、旧成功・キャンセル・時間切れを転用しない。

## 途中停止からの明示的な続行

予算終了で `failure-native-hH-gG-*.json` が残り、未完了の `native-raw/` に封印済み checkpoint がある場合だけ使う。失敗した初回経路は保持する。

```powershell
python "$Study/study_c02.py" --mode continue-native --h 2 --steps 32 --timeout 14400 --output $Run
```

初回を `native-interrupted-raw/` に保持し、最後の checkpoint を別プロセスが `native-resumed-raw/` へ続行する。重複する全状態の完全一致を要求し、比較用 prefix/suffix の由来を `native-raw/continuation-provenance.json` に記録する。完了済み経路を上書きするコマンドではない。初回の時間切れを成功へ変更しない。

## 封印証拠の再計測

```powershell
python "$Study/evidence_c02.py" --mode replay --output work/NEW-c02-sealed-replay
```

全 part・ZIP 全体・全エントリーの SHA を照合し、新規展開した封印ソースを `-I -X utf8` の別プロセスで実行する。7 測定 JSON の完全一致を要求する。全数値の再計測であり、新規求解・checkpoint の実再開の代用ではない。新規求解は前節のコマンドを使う。展開は数十 GB・数十分を計画し、既存の出力先には書かない。

厳密な封印再計測の JSON 完全一致は今回 Linux で確認する。Windows の BLAS/LAPACK 等の丸めで完全一致しない場合は、その差と失敗を保持し、許容値を追加して成功へ変更しない。厳密再計測には同じ依存版の Linux/WSL 環境を用いる。新規 Windows 求解の系列評価と、旧 Linux 封印値のビット再現は区別する。

## PC：細かい時間側での同時 3 水準

初期 16/32/64 の同時系列は粗い全体刻みの影響が大きい。前節の追加計算後、12/64 → 96/128 → 768/256 を次の小さな集計スクリプトで別に比較する。外挿・補正値は作らない。全点集合は同じ 672 接触点・23652 材料点・1458 節点である。

```powershell
@'
import sys
from pathlib import Path
sys.path.insert(0,str(Path('interim/2026-11-19-wet-joint-spatial-convergence').resolve()))
import mapping as m
from measure_c02 import triple
c=m.c
source=Path(sys.argv[1]).resolve()
target=Path(sys.argv[2]).resolve()
target.mkdir(parents=True,exist_ok=False)
maps=m.Maps([c.read(source/f'mesh-h{h}.json') for h in range(3)],c.read(source/'stages.json'))
result=triple(source,target,maps,[(0,64),(1,128),(2,256)],'native','simultaneous-fine')
c.write(target/'common-coordinates.json',maps.coordinates())
c.write(target/'measurements.json',dict(passed=False,c02_judgement='未達',series=result,extrapolated_field_created=False))
'@ | Set-Content -Encoding utf8 work/c02-fine-sim.py
python work/c02-fine-sim.py $Run work/NEW-PC-c02-fine-sim-measure
```

同じ関数へ 12/128 → 96/256 → 768/512 を指定すれば次の同時系列になる。`passed=false` は自動で受入へ変えず、補間不確かさ・固定時間の空間系列・開閉境界・両実装を併せて再評価する。
