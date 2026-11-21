# P6-01 再照合とユーザーPC手順

保存先の2026-12-06は依頼ラベル、実施日は2026-11-12（日本時間）。起点は
`8342762bd714e33523ce5fef772572699fc21881`、PR #13は未マージ。
以下の時間/RAMは既存手順の計画値であり、保証値ではない。
未存在の出力先を使い、時間切れ・キャンセル・未収束・初回失敗を保存する。

## 台帳と監査だけの再実行

監査はPython標準ライブラリだけを使い、ソルバーを呼ばない。Gitの固定blobも照合するので
`.git`と根拠ファイルを含む配置で実行する。指定起点のcommit/treeオブジェクトが必要で、shallow cloneはその起点まで取得する。数秒〜1分、RAM 1 GB程度。

```powershell
git fetch origin refs/heads/codex/plastic-spatial-convergence:refs/remotes/origin/codex/plastic-spatial-convergence
$rev = git rev-parse origin/codex/plastic-spatial-convergence
git -c core.autocrlf=false worktree add --detach D:\Work\p6-status-reconciliation $rev
Set-Location D:\Work\p6-status-reconciliation
git merge-base --is-ancestor 8342762bd714e33523ce5fef772572699fc21881 HEAD
if ($LASTEXITCODE -ne 0) { throw 'wrong branch ancestry' }
$P='interim/2026-12-06-p6-01-status-reconciliation'
New-Item -ItemType Directory work/NEW-p601-audit -ErrorAction Stop | Out-Null
py -3.13 "$P/audit_status.py" > work/NEW-p601-audit/audit-result.json
$auditExit=$LASTEXITCODE
py -3.13 "$P/audit_status.py" --table > work/NEW-p601-audit/ledger.ja.md
py -3.13 -m unittest discover -s $P -p 'test_*.py' -v 2>&1 | Tee-Object work/NEW-p601-audit/tests.log
$testExit=$LASTEXITCODE
@{audit_exit=$auditExit;test_exit=$testExit} | ConvertTo-Json | Out-File work/NEW-p601-audit/process.json
```

終了0は根拠整合の判定であり、P6数値受入ではない。`integrity_ok`、`row_count`、
`state_counts`、全`findings`を見る。整合errorは終了1、入力読取エラーは2。
warningの時点/範囲差も残す。同時点・同範囲の受入/未達はerrorで報告する。
日本語の自由文全体の意味判定は自動化しきれないため、登録した引用・scope・claimsを人も読む。
WindowsのCRLF化でSHAが変わる場合は合格にせず、LFを保った別の新規配置で実行する。
比較器に改行・丸め許可を追加しない。2文書の追記は台帳に保存した**完全一致の末尾**だけ許可し、元本文は起点blobと照合する。

## Linux / Windowsの全件ソフトウェア検証

本記録のLinuxはRust 1.99.1、Python 3.13.17、固定依存の環境を使用。
数値計算の新しい系列ではなく、既存の全件試験を実行する。
WindowsはPython 3.13と`requirements-reference.lock.txt`で仮想環境を作る。
全ビルド/全件には20〜120分以上、空き30〜60 GBを計画する。
旧mean/spatialアーカイブだけで約6.4 GBある。fixture取得時間は別枠。

```powershell
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
py -3.13 -m venv .venv
.venv\Scripts\python -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw 'dependency setup failed' }
$env:CARGO_BUILD_JOBS='2'; $env:RUST_TEST_THREADS='2'
$env:OPENBLAS_NUM_THREADS='1'; $env:MKL_NUM_THREADS='1'; $env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'; $env:PYTHONIOENCODING='utf-8'
New-Item -ItemType Directory work/NEW-p601-check -ErrorAction Stop | Out-Null
cargo fmt --all -- --check 2>&1 | Tee-Object work/NEW-p601-check/fmt.log
$fmt=$LASTEXITCODE
cargo clippy --locked --all-targets -- -D warnings 2>&1 | Tee-Object work/NEW-p601-check/clippy.log
$clippy=$LASTEXITCODE
cargo build --locked --release --bins --examples 2>&1 | Tee-Object work/NEW-p601-check/build.log
$build=$LASTEXITCODE
cargo test --locked --release --all-targets --no-fail-fast 2>&1 | Tee-Object work/NEW-p601-check/rust-tests.log
$rust=$LASTEXITCODE
cargo test --locked --release --doc 2>&1 | Tee-Object work/NEW-p601-check/rust-doc.log
$doc=$LASTEXITCODE
@{fmt=$fmt;clippy=$clippy;build=$build;rust=$rust;doc=$doc} | ConvertTo-Json | Out-File work/NEW-p601-check/rust-exits.json
foreach ($g in @('python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence')) {
  .venv\Scripts\python tools/run_ci_python_tests.py --group $g --output "work/NEW-p601-check/$g" 2>&1 | Tee-Object "work/NEW-p601-check/$g-driver.log"
  # 終了1でもログとinventoryを保持し、次の群を実行する。
}
```

Linuxでは同じcargoコマンドと`python tools/run_ci_python_tests.py`を使う。
このクラウド固有の環境は`source /workspace/soil-environment/activate.sh`。
全IDの和集合・4群への一意割当を確認し、最終成功数を固有IDで数える。
初回と再試験を加算しない。fail/error/skip/ignored/filtered/xfailを全て別記する。

旧NC2群が厳密なsetUpClass検査で止まったら、各4件は未実行、各準備エラー1と記録する。
skipや参照・許容値変更で通さず、過去の8件未実行も保持する。
疎配置の不足なら、まず初回ログを保存し、`git ls-tree HEAD <path>`の元blobを
`git show HEAD:<path>`で復元する。SHA-256とGit blob SHAを記録して該当群を全件再実行する。
これらの全件成功とNC/P6の数値受入は別である。

## 一般3D空間収束：PCで行う最細系列

今回は新規の独立数値確認を追加していない。小型・平滑解や2水準差では一般3Dの残件を
覆えず、既存基準のまま最細時間/空間と独立系を揃えるにはPCの資源が必要。
次は既存C-02の固定水圧・MC弾性bulkという**限定した問題**の追加であり、これだけで一般3D完了とはしない。

```powershell
cargo build --locked --release --example wet_joint_spatial_probe
if ($LASTEXITCODE -ne 0) { throw 'probe build failed' }
$C='interim/2026-11-19-wet-joint-spatial-convergence'
$Run='work/NEW-PC-p601-c02'
.venv\Scripts\python "$C/study_c02.py" --mode init --solver target/release/examples/wet_joint_spatial_probe.exe --output $Run
if ($LASTEXITCODE -ne 0) { throw 'input generation failed' }
foreach ($h in 0,1,2) {
  foreach ($n in 128,256) {
    .venv\Scripts\python "$C/study_c02.py" --mode native --h $h --steps $n --timeout 28800 --user-pc --output $Run
    if ($LASTEXITCODE -ne 0) { throw 'native incomplete: preserve all files' }
    .venv\Scripts\python "$C/study_c02.py" --mode independent --h $h --steps $n --user-pc --output $Run
    if ($LASTEXITCODE -ne 0) { throw 'independent incomplete: preserve all files' }
  }
}
.venv\Scripts\python "$C/study_c02.py" --mode restart --h 2 --steps 256 --timeout 28800 --user-pc --output $Run
if ($LASTEXITCODE -ne 0) { throw 'real checkpoint restart failed' }
.venv\Scripts\python "$C/measure_c02.py" $Run work/NEW-PC-p601-c02-measure
```

12/96/768 Tet10、32/128/512接触GP、共通672接触点・23652材料点・1458節点。
力1e-8、真の線形/独立系1e-10、接触増分上限5e-5 m、既存の水圧積分と物理入力を保持する。
RAM 64 GB、空き100 GB以上、数時間〜数日。まずh=2,n=128の所要時間を見て256の予算を決める。
独立系はdenseなので時間とRAMが増える。勝手に近似系・許容値を変更しない。

見る値は`measurements.json`の`spatial.fixed128/fixed256`、`temporal_contamination`、
`implementation_comparison`、`restarts`と全差分/方向/補間不確かさ。
全7工程のi/N受理列、全場・残差・摩擦円/仕事/凍結を検査し、別PID再開は時間だけ除いた全Run一致が必要。
負次数・方向反転・時間経路/補間の混入が残れば未達。正次数だけで合格にせず、未定義の全履歴閾値を新設しない。
6144要素・512以上は[既存PC手順](../2026-11-19-wet-joint-spatial-convergence/REPRODUCE.ja.md)の
「さらに細かい空間水準」に従う。64〜128 GB/空き500 GB/1〜数日は未実測の計画値。

## Griffiths元5%・正確なFs・実断面SRM

e278dbfの48/300要素F=1は、元661点・F=1.3500000000000003の比較ではない。
[既存の詳細手順](../2026-12-04-griffiths-ilut-p6-02/REPRODUCE.ja.md)の入力整合・solver固定を先に行う。
元参照と同じF/材料/全y拘束・奥行き1 m・ny=1を保持する。

```powershell
$S=(Resolve-Path target/release/soil-fem-3d.exe).Path
$G='interim/2026-12-04-griffiths-ilut-p6-02'
$H='runs/NEW-PC-p601-griffiths-h'
$Ref='tests/data/griffiths_2d_reference.json'
.venv\Scripts\python "$G/pc_series.py" prepare $H --reference $Ref
if ($LASTEXITCODE -ne 0) { throw 'prepare failed' }
.venv\Scripts\python "$G/pc_series.py" preflight $H --reference $Ref
if ($LASTEXITCODE -ne 0) { throw 'fixed input mismatch' }
foreach ($level in @('xy1-ny1','xy2-ny1','xy4-ny1')) {
  $q="$H/$level"
  & $S solve-stages "$q/mesh.json" "$q/ground.json" "$q/case.json" "$q/result.json" 2>&1 | Tee-Object "$q/solve.log"
  if ($LASTEXITCODE -ne 0) { throw 'incomplete: no 5%/Fs verdict' }
}
.venv\Scripts\python "$G/pc_series.py" compare $H --reference $Ref --solver $S --output runs/NEW-PC-p601-griffiths-comparison.json
```

xy1/xy2/xy4=1200/4800/19200 Tet10。xy1は10〜60分、xy2は0.5〜6時間、xy4は数時間〜1日超、
RAM 8〜32 GB/空き30 GB以上。ILUTは反復数が減っても遅くなり得る。
全受理力残差<=1e-8、元の真の線形<=1e-10、全工程completeが必要。
元の中央面L2 15.582/15.886/15.593/15.672%不合格を保持し、比較は元5%で行う。
連続側は同条件の面内h系列と比較器の条件を満たすまで判定不能。

正確なFsは同手順の「Fsと実断面」を使い、`run-srm-search-job`と`verify-exact-fs`で
各率の元入力からの履歴・独立釣合い排除を照合する。数時間〜数日。
不収束・上限・時間切れは独立排除ではない。元の2つの一様Fs区間とその中点を正確な斜面Fsにしない。
実断面は実CAD/メッシュ/材料/荷重/出典/資源予算が揃うまで未実行。入力のない実断面コマンドを仮造しない。

## NC・他のP6残件

- NCは[PC-3](../2026-12-05-nc-independent-arc-residual/REPRODUCE.ja.md)と既存`rerun-pc.ps1`を使用。
  400/800前史、材料内4/8、元のfloat64残差<=1e-10、全場/モード/荷重、診断restartを確認する。
  前史各5〜120分、局所数十秒〜数分が目安。PC-2の超過/0点は残し、Linux8点をPC成功へ転用しない。
  e940badは未確定・強く支持、NC-01全撤去は未達。案2は未承認なので実装しない。
- A/B/Cの長系列は各最新interimのREPRODUCE（台帳sourceパス参照）を使う。
  Richards＋JOINT、真曲面、一般折れ目等の現行拒否は、仕様/承認/独立検証を揃えるまで保持する。
- MCC空間/奥行き、杭/旧降雨、一般10万要素、Native/実GUI/GPU、別PC/配布は
  [P6全体の整理](README.ja.md)のPC欄へ。未達や判定不能を「完了」に置き換えない。
