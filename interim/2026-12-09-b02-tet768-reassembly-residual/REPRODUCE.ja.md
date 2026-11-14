# B-02 768 Tet10 残差・共通点差の診断を再現する

保存先の `2026-12-09` は依頼されたラベルであり、実行日ではない。
PR #13 と参照 draft PR をマージしない。元の B-02 記録、NC、許容値、参照値、hash/snapshot は変更しない。
以下はリポジトリルート、Rust 1.99.1 と requirements-reference.lock.txt の固定依存を前提とする。
出力名 `NEW-*` は毎回未存在の名前へ変える。時間切れ、未収束、拒否、非ゼロ終了を合格にしない。

## Windows PC：全件検証と新規求解

```powershell
git fetch origin codex/plastic-spatial-convergence
git worktree add --detach D:\Work\b02-residual-repro origin/codex/plastic-spatial-convergence
Set-Location D:\Work\b02-residual-repro
git rev-parse HEAD
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
python -m venv .venv
.venv\Scripts\python -m pip install -r requirements-reference.lock.txt
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:CARGO_BUILD_JOBS='2'
$env:RUST_TEST_THREADS='2'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
$B='interim/2026-12-01-tri6-offset-composite-b02'
$D='interim/2026-12-09-b02-tet768-reassembly-residual'
.venv\Scripts\python "$B/validate.py" --part rust --output runs/NEW-b02-residual-rust
.venv\Scripts\python "$B/validate.py" --part python --output runs/NEW-b02-residual-python
.venv\Scripts\python "$D/build.py" work/NEW-b02-residual-build
.venv\Scripts\python "$D/run.py" work/NEW-b02-residual-core --binary target/release/b02-residual-observer.exe --levels 0 1 --steps 4 8 16
.venv\Scripts\python "$D/run.py" work/NEW-b02-residual-primary --binary target/release/b02-residual-observer.exe --levels 2 --steps 4
.venv\Scripts\python "$D/build_frames.py" work/NEW-b02-frame-build
.venv\Scripts\python "$D/run_frames.py" work/NEW-b02-residual-fine --cli target/release/soil-fem-3d.exe --observer target/release/b02-frame-observer.exe
.venv\Scripts\python "$D/join.py" work/NEW-b02-residual-formal --sources work/NEW-b02-residual-core work/NEW-b02-residual-primary work/NEW-b02-residual-fine
.venv\Scripts\python "$D/analyze_matrix.py" work/NEW-b02-residual-formal runs/NEW-b02-residual-analysis
.venv\Scripts\python "$D/localization.py" work/NEW-b02-residual-formal runs/NEW-b02-residual-analysis runs/NEW-b02-localization.json
```

完全 checkout の旧証拠が必要。必要なら `git sparse-checkout disable` を行う。
全件検証は fmt、locked Clippy 全ターゲット・warnings-as-errors、locked release build、
Rust 全ターゲットと doc、Python 全発見集合4群を実行する。各 process JSON・ログ・ID inventory を保持する。
旧 NC mean/spatial の Linux 準備エラーによる4+4件未実行を成功やskipにしない。
Windows の元の OS 分岐を変更せず、今回の実結果を旧記録と分ける。
初回ビルドを含む全件検証は20〜120分をPC未実測の目安とする。

上記は12/96/768 Tet10 × 荷重刻み4/8/16の9条件。
最初の7条件は本体観測付き新規求解、入力だけの独立新規求解、観測なし本体の対照求解。
本体と対照の Run 全体は `elapsed_seconds` だけを除き一致する必要がある。
768/8,16は本体CLIで元の40/80増分を最後まで新規求解し、完成Runの5工程端点を
別の公開API観測器で再評価する。独立側は全40/80受理点を保存する。
この端点観測を本体全中間場の観測に読み替えない。
生産の manifest/lock/src は変えず、別 manifest から公開 API の診断実行版をビルドする。
`native/accepted.jsonl`、`independent.jsonl`、`native/operators.json`、`native/run.json`、
`plain/run.json`、入力と process JSON が生出力である。空き20 GB、RAM 8〜16 GB、
60〜300分をPC未実測の目安とする。完全checkoutと履歴の空きは別に80 GB以上を見込む。768/16の受理点ごとの再開観測は通常の一括求解より重い。
小さい動作確認だけなら `--levels 0 --steps 4` を使用できるが、9条件の代わりにはしない。

`analyze_matrix.py` は生出力を読む観測であり、新規求解に数えない。
768/4の全受理点と、他の各条件の元残差最大点・載荷終点を詳細再組立する。
元残差と本体対独立の全場比較は、最初の7条件では全受理点、768/8,16では
本体の5保存端点を対象にする。比較ごとの `common_state_points` を確認し、
20点と5点の非重み付きL2をそのまま同じ空間次数の計算に混ぜない。
`*-residual.json` で次を確認する。

- `residuals.original` と既存1e-10との関係。本体・独立を区別する。
- `soil_forward/reverse/kahan/extended`、`structure_split/relative/compensated/extended`、`both_extended`。
- `stored_stress_minus_elastic_from_total_u`、土・板・梁の寄与、節点/方向/座標付き `maximum_contributors`。
- 同じ最大自由度の `decimal80_residual_kn`。80桁は入力binary64値の算術再評価であり、連続体の80桁解ではない。
- `native_solver_attempts` の非線形力残差・各真の線形残差と設定。
- 独立 `all_direct_linear_observations` の通常/拡張精度残差、後退誤差、捨てた1回改良の差。

`longdouble` の仮数ビット数は出力に保存する。Windowsで53ビットなら拡張精度には数えない。
最大自由度の Decimal80 は同じ環境でも利用可能。全ベクトルの64ビット仮数の拡張演算は
Linux/WSL2上で同じ手順を実行し、`extended_mantissa_bits=64` を確認する。
Linuxでは `.venv/bin/python` と拡張子なしの実行版、`export NAME=value` に読み替える。

`field-refinement.json` の `space` と `load` を分け、せん断の回復点補間/解からの直接評価/未投影値、
瞬間仕事と4荷重成分の差を見る。未投影値は別の観測であり、本体の生せん断の代用ではない。
`*-fields.json` の回復点 fit 誤差、元生値と直接評価の差、要素 owner 数、仕事再合成誤差も確認する。
平滑化、外挿、Richardson補正、新しい合格許容値は導入しない。

## 6144 Tet10：ユーザーPCのみの大規模手順

この系列はクラウドでは未実行。本体単独の診断とする。
既存の独立直接参照は6500土自由度上限を超えるため、上限を外して実行しない。

```powershell
.venv\Scripts\python "$B/prepare_pc.py" runs/NEW-b02-residual-6144 --level 3 --steps 16 --offset 0.05
$L='runs/NEW-b02-residual-6144'
& target/release/soil-fem-3d.exe solve-stages "$L/mesh.json" "$L/ground.json" "$L/stages.json" "$L/run.json"
& target/release/soil-fem-3d.exe checkpoint-stages "$L/mesh.json" "$L/ground.json" "$L/stages.json" 40 "$L/checkpoint-40.json"
& target/release/soil-fem-3d.exe resume-stages "$L/mesh.json" "$L/ground.json" "$L/stages.json" "$L/checkpoint-40.json" "$L/resumed.json"
.venv\Scripts\python -c "import sys;sys.path.insert(0,'$B');from study import exact;from cases import read;assert exact(read('$L/run.json'),read('$L/resumed.json'))"
```

RAM 16〜32 GB、空き30 GB以上、各求解1〜8時間、再開照合まで3〜24時間をPC未実測の目安とする。
`scope.json` の節点・要素・参照上限、`run.json` の `complete/failure_reason`、
各 `attempts[].relative_force_residual` と `linear_solves[].true_relative_residual`、
`frames` の支持反力・荷重仕事、板生回復点のせん断・エネルギー、梁先端変位/接続回転を見る。
点荷重・接続位置・支持・物性・偏心・既存ゲートは同一に保つ。
完全な共通物理点系列を追加する場合は `convergence.py` の既存定義を保ち、
生回復値と直接評価を両方保存する。未到達の解や未対応な独立解を埋めない。
40は初期平衡16＋設置16＋載荷途中8の受理数で、物理秒ではない。

同じPCで求解した同じ固定実行版の checkpoint だけを再開する。
クラウド/Linuxの封印 checkpoint をWindowsの別実行版へ渡さない。
容量・幾何・保存整合ガードに掛かったら拒否結果をそのまま保存する。
GUI/GPU、実CAD、塑性土、水、JOINT、アンカーとの複合の受入には転用しない。

## 生証拠の全SHA再検証

```powershell
.venv\Scripts\python "$D/package.py" --verify "$D/evidence-manifest.json" --scratch work/NEW-b02-residual-verify
```

数分〜20分、結合ZIPの2倍以上の空きを見込む。
全part/memberのサイズとSHAが一致することを確認する。これは再求解ではない。
B-02全体は `passed=false`、NC-01、P6、A-01〜A-04、B-01/B-03、C-01〜C-04も未達のまま。
push後のWindows CIは、該当SHAの実結果を別途確認するまで未確定とする。

## 停止条件だけの対照（正式系列の代用にしない）

```powershell
.venv\Scripts\python "$D/stopping_controls.py" work/NEW-b02-residual-formal/mesh-2 runs/NEW-b02-stopping-controls --cli target/release/soil-fem-3d.exe --budget 180
```

同じ768/4入力のコピーに、(1) 線形停止1e-12だけ、(2) 線形1e-12＋力停止1e-10を指定する。
元の入力・許容値・既定値・再組立1e-10ゲート・参照・NC設定は書き換えない。
両方とも厳しくする入力実験で、B-02の合格化・限定受入の拡張には使わない。
`stopping-controls.json` の全設定差、complete、実測時間、受理/棄却、力残差と
真の線形残差を確認する。時間切れはresource_timeout、未完了として保存する。
出力Runの有無とfailure_reasonを確認し、予算終了を数値的停止・成功にしない。
想定各1〜10分。180秒で未完了なら、PCで別出力名・十分な予算にして追試する。

## 768の全中間場を8/16刻みでも観測する長時間系列（PCのみ）

```powershell
.venv\Scripts\python "$D/run.py" work/NEW-b02-residual-all-intermediate --binary target/release/b02-residual-observer.exe --levels 2 --steps 8 16
```

受理点ごとの厳密な保存再開・経路認証のため、一括求解より大幅に重い。
合計1〜6時間、RAM 8〜16 GB、追加空き20 GBをPC未実測の目安とする。
`native/accepted.jsonl` がそれぞれ40/80行、`comparison.json` の受理点数が40/80、
`native/run.json` のcomplete=true、`plain/run.json` と計測時間だけを除く厳密一致、
`source-binding-before/after.json` の一致を確認する。
今回のクラウドの5端点結果から未保存の全中間場の最大残差を推定しない。

## 全件試験の今回件数と環境回復

今回の最終固有IDはRust461成功、Python665件中663成功・2失敗、未実行/skip 0。
旧NC mean/spatial各4件は今回成功したが、過去の8件未実行・準備エラー記録は保持する。
P6-01/P6-03の固定台帳検査2件は、指定起点でも後続PC3 validation未登録を拒否し、
本追記のmeasurements/validationも未登録と判定する。古い台帳を変えて通さない。
したがって全件validationのpassed=falseを保持する。

LinuxではPowerShell 7の実pwshもPATH上に配置する。未配置によるskipで終わらせず、
準備後に元の同じIDを再実行する。今回の23 ID再実行は次で再現する。

```powershell
.venv\Scripts\python "$D/retry_tests.py" runs/NEW-b02-recovered-ids
.venv\Scripts\python "$D/final_validation.py" --rust runs/NEW-b02-residual-rust --python runs/NEW-b02-residual-python --retry runs/NEW-b02-recovered-ids --output runs/NEW-b02-final-validation.json --scratch work/NEW-b02-original-validation.json
```

全件グループのログを消さず、最終の元665 IDに対して最新の実結果を対応付ける。
再実行23件を発見件数へ重複加算しない。PowerShellの10項目は既存Pythonメソッド内の検査である。
疎配置を使う場合、新しいP6-01のledger.json等も元blobから復元する。
既存のJSON・許容値・Linux厳密比較・snapshot/hashを変更しない。

生証拠ZIPは約1.09 GB、展開ファイル計約3.66 GB。検証用の結合/展開先も別に確保する。
補足のvalidation/preparation/handoff-interrupted.zipは意図的に中断した重複観測であり、
完成Run・未到達工程・正式な全履歴観測の証拠として使わない。
