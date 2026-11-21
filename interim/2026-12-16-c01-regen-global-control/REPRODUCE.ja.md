# C-01 再生参照・2接線・全体平衡対照の再現

PR #13 / `codex/plastic-spatial-convergence`。マージしない。フォルダ名は依頼された記録識別子で、実施日は2026-11-14（日本時間）。PC実行はしていない。
リポジトリルートのPowerShell 7で実行する。`NEW`は存在しない出力名へ置換する。

## 準備と旧証拠の保全

既存の作業を上書きせず、新しいcloneを推奨する。Git履歴と作業treeは50 GB規模になり得る。
`core.autocrlf=false`はこのcloneの証拠・ソース生バイトのため。文書比較試験はCRLFにも対応するが、封印SHA自体は正規化しない。

```powershell
git clone --filter=blob:none --single-branch --branch codex/plastic-spatial-convergence --config core.autocrlf=false https://github.com/tak063495-prog/solidcraft_soil_analysis.git NEW-C01
Set-Location NEW-C01
git merge-base --is-ancestor e8aebaba9f08655b626f8b6ac9f5a9a1a7b50287 HEAD
# 完全再現には旧証拠も必要。空き容量を確認してから全証拠を取得する。
git sparse-checkout disable
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
cargo build --locked --release --example wet_joint_local_path_probe
```

Linuxでは `.venv/bin/python`、実行版の拡張子なしを使う。
全件試験用のアーカイブ集合は既存 `.github/workflows/fem-ci.yml` の各ジョブの和集合。
`reference400-fixture.zip`も含め、ファイル名が`evidence.zip`でない証拠を落とさない。
新しい全体平衡対照は診断用Python参照だけに存在し、本体CLIの入力やJSON schemaにはオプションを加えていない。
`c01_regen_reference._reference(..., diagnostic_schedule=None)`が既定。指定時だけ別の平衡問題を解く。

## 同じ21保存jumpの×64/×128再生

```powershell
python interim/2026-12-16-c01-regen-global-control/c01_regen_finite.py --prepare --inputs work/NEW-c01-saved --output work/NEW-c01-finite
```

12-11/12-12/12-13のpart・索引・ZIP SHAを検証し、12-15と同じ保存jump・実分割数を取り出す。
×64/×128の全21経路を再生し、旧12-15の封印traceから×16/×32を検証して抽出する。
出力は`x64/`、`x128/`の`decomposition.json`・全GP`point-increments.jsonl`・`finite-reference.json`、旧traceの`old32/`、全4水準の`finite-summary.json`。
元の51比較の定義を保持する。各(i)〜(iii)ベクトルの変化、全差への符号付き投影、交差項、64→128の差を確認する。
(iv)の真の端点間経路・連成効果は独立ベクトルとして同定不能であり、0や厳密な誤差上界とは読まない。
想定10〜60分、RAM4〜12 GB、出力・展開4〜8 GB（Git履歴は別）。機種に依存しWindows実測ではない。

## 解析解のある回転経路

```powershell
python interim/2026-12-16-c01-regen-global-control/c01_regen_rotating.py --solver target/release/examples/wet_joint_local_path_probe.exe --output work/NEW-c01-rotating
```

3経路、外側2/4/8/16/32/64、内部1/2/4/8の72経路と、その固定3D回転コピー72経路。
入力と解析解は`tau=c(t)n(t)`、`ds/dt=gamma*n(t)`、`dx=dtau/kt+ds`から構成する。
`request.json`・`response.json`・回転コピーと全増分`rotating.json`、軽い`summary.json`を出す。
見る値はすべり・traction・増分D・累積Dの誤差、全系列の次数／単調性、散逸増分の非負性、すべり増分方向・traction方向、固定座標回転の整合性。
円弧を外側端点の弦で置換するため、内部数だけを増やした時に残る幾何学的経路欠損も比較する。
想定数秒〜1分、RAM1 GB未満、出力20〜80 MB。動く法線に対する有限回転輸送の検証ではない。

## 同じ入力からの新規全体平衡と開閉固定対照

```powershell
python interim/2026-12-16-c01-regen-global-control/c01_regen_global.py --output work/NEW-c01-global --levels 64 128 256
```

54節点・12 Tet10・324材料GP・32接触GP、固定水位0.75 m、全7工程、入力だけの3D別組立を使う。
通常64/128/256を新規求解し、256の全GP保存gapポリラインの根から開区間を作る。
その同一状態列を固定する対照64/128/256を全体Newtonで新規求解する。
既存の接触内上限5e-5 m、切替二分割4、Newton60、力1e-10、真の線形1e-10を保持する。新しい刻み制御はない。
本体の通常力基準1e-8より厳しい既存独立参照の力基準1e-10を変えていない。

固定開区間では力0、固定閉区間では法線ばねを両側へ延長し、摩擦上限は`mu*max(-qn,0)`。
閉鎖を固定した点が開く向きのgapを持つ時には引張反力も出る。これは指定状態列を保つための**変更問題**であり、元の片側接触問題の解ではない。
開閉時刻で局所積分を切ることは、処方した状態履歴の積分であり新しい全体刻み制御ではない。
通常／対照の全場・入力・力と真の線形残差・途中失敗を各フォルダへ保存する。合否には使用しない。
`fine-open-schedule.json`、`comparison.json`の元のL2・差方向余弦、介入差の交差内積を見る。
想定3〜30分、RAM1〜4 GB、出力200〜500 MB。Windows実測ではない。
局所再生を全体平衡の代用にしていない一方、この有限参照の介入を(ii)/(iii)の普遍的で一意な因果比率とはしない。

## PCのみ：元の真の段階間差

以下はクラウドで新規実行せず、ユーザーPC向けに残す。既存の上限・基準のまま、失敗も保存する。

```powershell
python interim/2026-12-13-c01-cause-diagnosis/c01_cause_diagnostics.py --extract switch --output work/NEW-PC-switch-saved
python interim/2026-12-13-c01-cause-diagnosis/c01_cause_diagnostics.py --extract path --output work/NEW-PC-path-saved
foreach ($n in 512,1024) {
  python interim/2026-12-13-c01-cause-diagnosis/c01_cause_experiment.py --control switch --steps $n --contact-order 4 --pc-series --output work/NEW-PC-reference/switch-n$n-q4
  # 非ゼロ終了の場合はfailure.jsonと途中fieldを保存し、完了扱いにしない。
}
python interim/2026-12-13-c01-cause-diagnosis/c01_cause_compare.py --saved-switch work/NEW-PC-switch-saved/paths --saved-path work/NEW-PC-path-saved/paths --experiments work/NEW-PC-reference --reference-levels 256 512 1024 --output work/NEW-PC-comparison
python interim/2026-12-13-c01-cause-diagnosis/c01_cause_experiment.py --control path --steps 512 --contact-order 4 --output work/NEW-PC-reference/path-n512-q4
cargo build --locked --release --example wet_joint_path_probe
python interim/2026-12-12-c01-switch-event-localization/study.py --mode init --solver target/release/examples/wet_joint_path_probe.exe --output work/NEW-PC-full
python interim/2026-12-12-c01-switch-event-localization/study.py --mode native --pc-series --steps 256 512 1024 --timeout 14400 --output work/NEW-PC-full
python interim/2026-12-12-c01-switch-event-localization/study.py --mode independent --pc-series --steps 256 512 1024 --output work/NEW-PC-full
python interim/2026-12-12-c01-switch-event-localization/measure.py work/NEW-PC-full work/NEW-PC-full-measure --steps 256 512 1024
```

独立512/1024は想定5〜60分以上・RAM16 GB・5〜15 GB、元の全場系列は1〜6時間以上・RAM16 GB・30〜60 GB。
見る値は元の全18場の全履歴L2、絶対差方向、次数、力／線形残差、全7工程完了。
誤差制御512が既存2000試行／工程上限に達した失敗を保持し、上限を変えない。1024も完了した場合だけ段階間差を比較する。
外挿、Richardson補正、許容値変更、定義変更で合格を作らない。保存済み512も1024も厳密解とは呼ばない。

## 全件検証と封印・台帳

```powershell
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release
cargo test --locked --release --all-targets
cargo test --locked --release --doc
foreach ($group in 'python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence','python-c01-switch','python-c01-path','python-c01-cause') {
  python tools/run_ci_python_evidence_shards.py --group $group --output work/NEW-validation/$group
  if ($LASTEXITCODE -ne 0) { throw "failed: $group" }
}
python interim/2026-12-16-c01-regen-global-control/c01_regen_evidence.py
python interim/2026-12-06-p6-01-status-reconciliation/audit_status.py
```

Python全発見IDと7群実行IDの和集合・重複0を確認し、成功／失敗／未実行／skipを`validation.json`へ正確に記録する。
新しい試験は軽い通常群。重い再生・全体対照はCIで再計算せず、48 MiB分割アーカイブのSHAを検査する。
C-01の既存3シャード・75分予算とNC関連コード・設定は変更しない。
`acceptance.json`と`validation.json`は先に不変commitへ保存し、既存`register_sources.py`で固定commit・SHA-256・Git blob・全boolean/null pointerを台帳へ追記する。
旧台帳を再生成せず、SOURCE_NOT_INVENTORIEDの例外・skipを追加しない。push後のWindows CIは未確定。

C-01は未達、passed=false。
