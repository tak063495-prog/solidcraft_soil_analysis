# A-03クラウド追補の再現

PR #13、`codex/plastic-spatial-convergence` の89e9dfc以降。マージしない。
2026-12-30は保存先ラベル。基準は[SPEC](SPEC.ja.md)と旧11-26/12-22のまま。
毎回未存在の出力先を使い、過去の成功・失敗・生証拠・封印を上書きしない。

## Windows PCの準備

GPD WIN Max 2、使用可能約48GBを前提とする。以下は**Windows未実行の手順**。
解析は54節点・12 Tet10だけ。見込み40GB超・1経路1日超の計画、768/6144要素の新手順はない。
一般空間収束でこの上限を超える検証が必要な部分は**規模上限のため確認不可・未達**とする。

```powershell
git status --short
git config --local core.autocrlf false
git fetch origin codex/plastic-spatial-convergence
git merge-base --is-ancestor 89e9dfc26cbe444199179d7faa1edb1437849309 HEAD
if ($LASTEXITCODE -ne 0) { throw '指定起点の子孫ではありません' }
$Py = '.venv/Scripts/python.exe'
if (-not (Test-Path $Py)) { py -3.13 -m venv .venv }
& $Py -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw '固定依存の配置失敗' }
$env:CARGO_BUILD_JOBS='2'
$env:RUST_TEST_THREADS='2'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
if ($LASTEXITCODE -ne 0) { throw 'toolchainの配置失敗' }
$D = 'interim/2026-12-30-a03-wet-soil-birth-cloud-followup'
cargo build --locked --release --all-targets
if ($LASTEXITCODE -ne 0) { throw '本体ビルド失敗' }
& $Py "$D/a03follow_build.py" --output work/NEW-a03follow-build
if ($LASTEXITCODE -ne 0) { throw '隔離診断ビルド失敗。ログを保持' }
```

生産Rustはコピーして読むだけで、切替は専用ビルドにしか入らない。
現在のhelperは補助関数をcfg(test)モジュール前へ配置するv2。保存した主求解はv1で、
元版の封印を保持したまま4基準Runの厳密一致を別検査した。指紋を付け替えない。
隔離コピーのall-targets clippyも行う場合は、validation/isolated-fixture-restoration.jsonの14個の
JSONを元のGit blob/SHAどおりsource/examplesへコピーする。fixture欠落をskipで回避しない。
ビルドはメモリ2〜8GB、合計10〜90分の保守的な予算（Windows未実測）。
通常checkoutには旧証拠を含む十分なディスク空きが必要。今回の取得・SHA照合と
A-02のgit-recovery記録を参照し、欠落アーカイブをskipや合成値で代用しない。
partial cloneでは固定出典のGit blobも必要。P6監査の前に取得し、元のGit blob/SHAで照合する。
新partは通常python-core群の実sparse規則に含まれることをvalidationの契約で確認する。

## 3つの小型全体平衡対照

```powershell
$env:A03FOLLOW_PROCESS_SECONDS='3600'
try {
  & $Py "$D/a03follow_run.py" --build work/NEW-a03follow-build --output runs/NEW-a03follow-PC
  if ($LASTEXITCODE -ne 0) { throw '数値未完了。途中状態とログを保持' }
  & $Py "$D/a03follow_report.py" --run runs/NEW-a03follow-PC --output runs/NEW-a03follow-PC-primary.json --vectors runs/NEW-a03follow-PC-vectors.json
  if ($LASTEXITCODE -ne 0) { throw '集計未完了。求解成功に読み替えない' }
  & $Py "$D/a03follow_surface.py" --formal runs/NEW-a03follow-PC --output runs/NEW-a03follow-PC-surface
  if ($LASTEXITCODE -ne 0) { throw '静水圧専用面則の順序対照が未完了' }
  & $Py "$D/a03follow_surface_report.py" --formal runs/NEW-a03follow-PC --surface runs/NEW-a03follow-PC-surface --primary runs/NEW-a03follow-PC-primary.json --output runs/NEW-a03follow-PC-report.json --vectors runs/NEW-a03follow-PC-surface-vectors.json
  if ($LASTEXITCODE -ne 0) { throw '全面則の集計未完了' }
} finally { Remove-Item Env:A03FOLLOW_PROCESS_SECONDS -ErrorAction SilentlyContinue }
```

2形状×2初期応力のN4/8/16、各N8で独立GL5/GL7。
同じwide/k05/N8/M1で要素順E・面宣言順F・Tri6面積分点順Qの全8組合せ。
同じwide/k05のN4/8/16/32×M1/2/4を両実装で求解する。
Mは接触局所仕事の最小分割数で、既存の必要分割上限・候補時の再分割を無効化しない。
全体格子はi/Nで入力から固定し、独立側へ本体の受理列・解・接線を渡さない。
通常CLIは4つのN8基準入力を別に解き、専用ビルドの既定と比較する。
観測なしwholeも全本体入力で解く。再集計やローカルjump再生は新規全体求解に数えない。
主行列の独立Qは接触/GL25面則までで、追加4経路は専用静水圧ray則の内外点も逆順にする。
部分Q対照を保持したまま、全8組合せを完成させる。正式64経路・2,816状態、
観測なしと通常CLIも含め96経路・4,224状態を別計数する。

|出力|確認する量|
|---|---|
|plan.json、source-binding.json、source-snapshot/、isolated-build-binding.json|実行前の行列・ソース・版SHA。実行中の差替えなし|
|各入力のmesh/ground/case.json、spec.json、*-process.json、*.log|入力・PID・終了コード・時間切れ。数値未完了と範囲拒否を区別|
|native/run.json、accepted.jsonl、constitutive.jsonl、preparations.json|全受理状態・材料履歴・全体Run・準備保持力と公開再組立|
|native-fields.jsonl、independent-fields.jsonl、*-layout.json|全54場。物理要素の接続と面内形状関数で対応付け、最近接点を使わない|
|報告JSONのraw/work_differences|接触仕事、数値散逸、弾性エネルギー、荷重仕事、反力仕事、自由残差、支持込み残差の生値|
|geometry/ordering/cross、差ベクトルJSON|GL感度、順序主効果、N×Mの主効果/交互作用、2つの3格子組の観測次数|

PC全行列の予算は10〜90分、ピーク概算1〜6GB、各プロセス上限1時間（Windows未実測）。
再開の内容再演算を含むためN比例の速度保証はしない。予算内で完了しなければ数値未完了・未達を保持する。
一般塑性/自由水面/大規模空間収束を追加するためのガード解除・メモリ上限解除は行わない。

元の変位/応力/traction相対1e-7、節点力1e-10 kN、拘束1e-10 m、仕事恒等式1e-10を保持する。
入力の力1e-11/線形1e-12、分母下限1e-12も維持。全54場の包括認定はfalse。
ノルム/仕事の集計は物理ID→受理点→成分の固定順、math.fsum。元のnumpy集計は別欄に残す。
Windowsの集計差を許容値で吸収しない。文書の末尾検査はCRLFをuniversal newlineで正規化する。
水量・流量・貯留・誕生移送はnullのままで、0や水収支合格にしない。

## 封印と全件検証

```powershell
& $Py "$D/a03follow_evidence.py" --verify $D --extract runs/NEW-a03follow-archive
if ($LASTEXITCODE -ne 0) { throw 'part/結合ZIP/entryのSHA不一致' }
& $Py "$D/a03follow_report.py" --run runs/NEW-a03follow-archive/formal --output runs/NEW-a03follow-reopened-primary.json --vectors runs/NEW-a03follow-reopened-vectors.json
if ($LASTEXITCODE -ne 0) { throw '保存場の再集計失敗' }
& $Py "$D/a03follow_surface_report.py" --formal runs/NEW-a03follow-archive/formal --surface runs/NEW-a03follow-archive/retained/a03follow-surface-full --primary runs/NEW-a03follow-reopened-primary.json --output runs/NEW-a03follow-reopened.json --vectors runs/NEW-a03follow-reopened-surface-vectors.json
if ($LASTEXITCODE -ne 0) { throw '全面則の保存場再集計失敗' }
& $Py "$D/a03follow_validation.py" rust --output runs/NEW-a03follow-validation/rust
& $Py "$D/a03follow_validation.py" python --output runs/NEW-a03follow-validation/python
& $Py "$D/a03follow_validation.py" collect --output runs/NEW-a03follow-validation
& $Py interim/2026-12-06-p6-01-status-reconciliation/audit_status.py
if ($LASTEXITCODE -ne 0) { throw 'P6固定出典/台帳の不一致' }
```

全件は既存のRust全ターゲット＋doc、Pythonの7群（testsとtools）をそのまま実行する。
新CI試験は小型契約・CRLF・順序・封印SHAだけで、rust-integration-1への重い追加は0。
全件の予算は30〜180分、メモリ2〜8GB（Windows未実測）。実測と全発見/実行IDを保存する。
単独スモーク・登録後の再確認を全件の固有件数へ二重加算しない。
封印の再読込・保存場再集計は新規求解ではない。外挿や最細値で合格を作らない。

Linuxでは `source /workspace/soil-environment/activate.sh` 後に同じPythonコマンドを使える。
当該push後のWindows CIは未確定。PCで必要な追加確認は同じ小型入力のOS間再現であり、
S1〜S5の未定義仕様をPCで解く手順ではない。
A-01/02/04、B-01〜03、C-01〜04、NC-01、P6も全体未達を維持する。

A-03は未達、passed=false。
