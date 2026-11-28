# P6-02 ILUT：再現と Windows PC の続行手順

2026-12-04は依頼ラベル、実施日は2026-11-12（日本時間）。PR #13はマージしない。
すべて未存在の出力名を使い、失敗/不収束/時間切れ/キャンセルも保存する。
以下のPC時間・RAMは未実測の予算目安。クラウド実測はREADMEの48/300要素だけ。

## 専用ワークツリー・ビルド・全件試験

```powershell
git fetch origin refs/heads/codex/plastic-spatial-convergence:refs/remotes/origin/codex/plastic-spatial-convergence
$rev = git rev-parse origin/codex/plastic-spatial-convergence
git worktree add --detach D:\Work\griffiths-ilut-p602 $rev
Set-Location D:\Work\griffiths-ilut-p602
git rev-parse HEAD | Out-File revision.txt
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
py -3.13 -m venv .venv
.venv\Scripts\python -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw 'dependency setup failed' }
$env:OPENBLAS_NUM_THREADS='1'; $env:MKL_NUM_THREADS='1'; $env:OMP_NUM_THREADS='1'
$env:CARGO_BUILD_JOBS='2'; $env:RUST_TEST_THREADS='2'
$env:PYTHONUTF8='1'; $env:PYTHONIOENCODING='utf-8'
$B='interim/2026-12-04-griffiths-ilut-p6-02'
$V='interim/2026-12-01-tri6-offset-composite-b02'
.venv\Scripts\python "$V/validate.py" --part rust --output runs/NEW-p602-rust
# Rustの終了コードを記録してからPythonも全件実行する。
$rustExit=$LASTEXITCODE
.venv\Scripts\python "$V/validate.py" --part python --output runs/NEW-p602-python
$pythonExit=$LASTEXITCODE
@{rust_exit=$rustExit;python_exit=$pythonExit} | ConvertTo-Json | Out-File runs/NEW-p602-tests.json
```

ビルド数分、全件20〜120分、RAM 4〜16 GB・空き30 GB以上を目安にする。
元Git blobの旧分割証拠をすべて配置する。Linux今回の疎取得では旧NC封印archive合計約6.4 GBの
取得がHTTP 503で止まった。欠損をskip・検査削除・参照書換えで成功へ変えない。
`processes.json`、全ログ、4群の`python-test-inventory.json`を確認する。
旧NC mean/spatialの各4件がsetUpClassで止まれば未実行8・準備エラー2を明記し、理由も記録する。
ファイル不足とLinux末尾浮動小数差は別の失敗理由。過去の未実行記録を保持する。

以下ではビルド済みCLIを複写して固定する。実行中のcargo再ビルドでSHAを入れ替えない。

```powershell
New-Item -ItemType Directory runs/NEW-p602-bin -ErrorAction Stop | Out-Null
.venv\Scripts\python tools/freeze_sources.py runs/NEW-p602-bin
if ($LASTEXITCODE -ne 0) { throw 'source freeze failed' }
Copy-Item target/release/soil-fem-3d.exe runs/NEW-p602-bin/soil-fem-3d.exe
$S=(Resolve-Path runs/NEW-p602-bin/soil-fem-3d.exe).Path
Get-FileHash $S | Format-List | Out-File runs/NEW-p602-bin/sha256.txt
```

## オプションとクラウドの小型対照

入力は既存CLIのCASE.jsonに指定する。新しい位置引数や環境変数はない。

```json
"linear_solver": {
  "method": "sparse_gmres_ilu0",
  "options": {
    "restart": 60, "maximum_iterations": 3000, "pivot_tolerance": 1e-14,
    "ilut": {"drop_tolerance": 0.0001, "fill_limit": 200}
  }
}
```

`ilut`省略/nullは従来ILU(0)。指定時は自然順ILUT、report.methodは`csr_gmres_ilut`。
dropは有限[0,1]、fillは整数1..4096で、**各strict L/U行それぞれの保持数**。
SciPyのfill倍率やCOLAMDではない。対角シフト・ピボット置換は行わない。
力1e-8、真の線形1e-10、Newton35、GMRES3000を維持する。

起点版とのビット比較には、別ワークツリーで指定先端をビルドした実行版が必要。
WindowsでLinuxの保存点を使い回したり、solver SHAを再封印して移植しない。

```powershell
git worktree add --detach D:\Work\p602-base 2a716084dc66754b681f29e062e5d9155d505647
cargo build --locked --release --manifest-path D:\Work\p602-base\Cargo.toml
if ($LASTEXITCODE -ne 0) { throw 'baseline build failed' }
.venv\Scripts\python "$B/study.py" runs/NEW-p602-small --cli $S --baseline D:\Work\p602-base\target\release\soil-fem-3d.exe --strides 5 2 --budget 600
.venv\Scripts\python "$B/build_observer.py" work/NEW-p602-observer
if ($LASTEXITCODE -ne 0) { throw 'observer build failed' }
.venv\Scripts\python "$B/restarts.py" runs/NEW-p602-small runs/NEW-p602-restarts --observer target/release/p602-observer.exe --cli $S
```

5〜30分目安。実測JSONの両complete、既存受理残差、全節点変位/全GP応力差を確認する。
前処理間の包括一致ゲートは「判定基準なし」のまま。数値差から閾値を作らない。
既定の全Runは認証した実行版SHAとelapsed_secondsだけ区別してビット一致すること。
同じ実行版の別PID再開はelapsed_secondsだけ除いて後続全状態/全Run一致が必要。
入力範囲/型、前処理/drop/fill食違い、checksum未更新、再hash改変をすべて拒否する。

## xy2：失敗時と同じ入力からの ILUT 新規実行

元失敗入力は実際の保存場所へ読み替える。下記は同じmesh/groundをバイト保持し、caseはILUT項目だけ追加する。
元失敗のF、ny、材料/積分/刻み、初期種、拘束が同じか必ず保存入力と照合する。
F=1と固定参照F=1.3500000000000003を混ぜない。

```powershell
$src='C:\soil-inputs\griffiths-xy2'
$out='runs/NEW-xy2-ilut'
New-Item -ItemType Directory $out -ErrorAction Stop | Out-Null
.venv\Scripts\python "$B/prepare.py" add-ilut "$src/case.json" "$out/case-ilut.json" --drop 0.0001 --fill 200
if ($LASTEXITCODE -ne 0) { throw 'input conversion failed' }
Get-FileHash "$src/mesh.json","$src/ground.json","$src/case.json","$out/case-ilut.json",$S | Format-List | Out-File "$out/sha256.txt"
$job="$out/job"
& $S prepare-durable-stages-job "$src/mesh.json" "$src/ground.json" "$out/case-ilut.json" 8 $job
if ($LASTEXITCODE -ne 0) { throw 'job preparation failed' }
$J="$job/solver/soil-fem-3d.exe"
$sw=[Diagnostics.Stopwatch]::StartNew()
& $J run-stages-job $job 2>&1 | Tee-Object "$out/run.log"
$code=$LASTEXITCODE; $sw.Stop()
@{exit_code=$code;elapsed_seconds=$sw.Elapsed.TotalSeconds} | ConvertTo-Json | Out-File "$out/process.json"
& $J stages-job-status $job | Out-File "$out/status-read.json"
```

4800 Tet10・129600 GP、数十分〜数時間、RAM 4〜16 GB・空き10 GB以上の目安。
今回の自然順ILUTは反復数が減っても因子分解で遅くなり得る。完了時間を保証しない。
`job/progress.jsonl`のNewton・LinearSolved・Accepted/Rejected、`status.json`、
`result.json`/保存checkpointを保全する。0.44付近のGMRES3000、元停止0.4717のNewton35、
相対力残差1.109e-7超過が改善したかと、**進捗1.0の完成**を区別する。
全成功線形solveの元残差<=1e-10、全受理力残差<=1e-8、geostatic.completeとrun.completeが必要。
0.4717を越えただけでは合格、Fs、破壊としない。

別端末で協調停止する場合：`& $J stop-stages-job $job`。
終了2は受理境界での休止。元ジョブを保持して新しい復旧先を作る。

```powershell
$new='runs/NEW-xy2-recovery'
& $J prepare-stages-recovery $job $new
if ($LASTEXITCODE -ne 0) { throw 'recovery preparation failed' }
& "$new/solver/soil-fem-3d.exe" run-stages-job $new 2>&1 | Tee-Object "$out/recovery-run.log"
```

再開時のdrop/fill変更やソルバ差替えは拒否する。変更実験は元入力から別の新規ジョブにする。
固定ケースのILUT checkpoint照合には受理経路の再計算も含まれるため、再開自体にも時間がかかる。
手動終了・時間切れ・取消は未完了。復旧ファイルの存在だけで全状態一致を主張しない。

## 固定元参照661節点と同条件の面内h系列

元参照に一致するF=1.3500000000000003の別系列を準備する。xy1/xy2/xy4は1200/4800/19200 Tet10。
奥行き1m・ny=1固定、全y変位拘束、元材料・ゼロ応力種・荷重・661節点写像を保持する。
既定用の元準備器/比較器は変更しない。ILUTラッパーは新規生成case/stagesのILUT項目だけ変え、
元manifestをoriginal-preparation.jsonに保持する。再生成による全入力・manifest照合も必須。

```powershell
$ref='tests/data/griffiths_2d_reference.json'
$D='runs/NEW-h-default-F135'; $I='runs/NEW-h-ilut-F135'
.venv\Scripts\python tools/prepare_slope_h_series.py $D --reference-2d $ref
if ($LASTEXITCODE -ne 0) { throw 'default preparation failed' }
.venv\Scripts\python "$B/pc_series.py" prepare $I --reference $ref
if ($LASTEXITCODE -ne 0) { throw 'ILUT preparation failed' }
.venv\Scripts\python "$B/pc_series.py" preflight $I --reference $ref
if ($LASTEXITCODE -ne 0) { throw 'ILUT input integrity failed' }
foreach ($root in @($D,$I)) {
  foreach ($level in @('xy1-ny1','xy2-ny1','xy4-ny1')) {
    $p="$root/$level"
    & $S validate "$p/mesh.json" 2>&1 | Tee-Object "$p/validate.log"
    if ($LASTEXITCODE -ne 0) { throw 'mesh validation failed' }
    $sw=[Diagnostics.Stopwatch]::StartNew()
    & $S solve-stages "$p/mesh.json" "$p/ground.json" "$p/case.json" "$p/result.json" 2>&1 | Tee-Object "$p/solve.log"
    $code=$LASTEXITCODE; $sw.Stop()
    @{exit_code=$code;elapsed_seconds=$sw.Elapsed.TotalSeconds} | ConvertTo-Json | Out-File "$p/process.json"
    if ($code -ne 0) { throw 'incomplete solve; preserve files, no convergence/Fs verdict' }
  }
}
.venv\Scripts\python tools/griffiths_2d.py series runs/NEW-h-default-comparison.json --reference-2d $ref --folder $D --solver $S
$defaultComparisonExit=$LASTEXITCODE
.venv\Scripts\python "$B/pc_series.py" compare $I --reference $ref --solver $S --output runs/NEW-h-ilut-comparison.json
$ilutComparisonExit=$LASTEXITCODE
```

xy1は各10〜60分、xy2各0.5〜6時間、xy4各数時間〜1日超、RAM 8〜32 GB・空き30 GB以上の未実測目安。
余力がない段階は開始せず「未実行」と記録する。系列途中の結果だけで全系列を合格にしない。
両比較器は完成Runと既存の独立力学/材料/VTU監査を要求する。監査省略オプションはない。
全661点xzの相対L2を奥行き0/.5/1で元5%に照合し、最大誤差10%・奥行き差1%を別に確認する。
比較0は全離散段階のL2合格、1は測定した5%不合格、例外/2は欠損・改変・未完了等の拒否。
連続側は判定不能、Fsは未評価。旧15.582/15.886/15.593/15.672%の不合格を上書きしない。

## Fs：未低減F=1の別入力で強度低減を探索

上のF=1.35入力をさらに低減して元の絶対Fsと報告しない。未低減F=1を別に準備する。
例はxy2、要求率1.1、最大率2、宣言絶対幅0.002、最大試行20、工程1の既存SPECを固定する。

```powershell
$F='runs/NEW-srm-F1'
.venv\Scripts\python tools/prepare_slope_h_series.py $F
if ($LASTEXITCODE -ne 0) { throw 'F1 preparation failed' }
$p="$F/xy2-ny1"
.venv\Scripts\python "$B/prepare.py" add-ilut "$p/case.json" "$p/case-ilut.json"
Copy-Item examples/srm-search.json "$F/specification.json"
$srmJob='runs/NEW-xy2-srm-job'
& $S prepare-srm-search-job "$p/mesh.json" "$p/ground.json" "$p/case-ilut.json" "$F/specification.json" 8 $srmJob
if ($LASTEXITCODE -ne 0) { throw 'SRM preparation failed' }
$Q="$srmJob/solver/soil-fem-3d.exe"
Get-FileHash "$F/specification.json",$Q | Format-List | Out-File "$F/fixed-spec-sha256.txt"
& $Q run-srm-search-job $srmJob 2>&1 | Tee-Object "$F/srm-run.log"
$srmExit=$LASTEXITCODE
& $Q srm-search-job-status $srmJob | Out-File "$F/status-read.json"
$st=Get-Content "$srmJob/status.json" -Raw | ConvertFrom-Json
& $Q verify-exact-fs "$srmJob/inputs/mesh.json" "$srmJob/inputs/ground.json" "$srmJob/inputs/case.json" "$srmJob/inputs/specification.json" "$srmJob/report.json" $st.report_sha256 "$F/exact-fs-assessment.json"
$exactExit=$LASTEXITCODE
```

各率を元入力から新規履歴で解き直すため、数時間〜数日・RAM 8〜32 GB・空き30 GB以上の目安。
独立排除のない既定SPECでは、探索完了でも正確なFsを得られない場合がある。
既存の1200要素F=1.45/F=2の仮想仕事上端をxy2や条件の異なるRunへ流用しない。
新たな上端候補はそのmesh/拘束/元外力と低減材に結び付け、仮想仕事/平均応力の独立監査が必要。
上端なし、trial_limit、最大率到達、非収束、時間切れは判定不能。
`verify-exact-fs`終了0かつexact_fs_for_this_case=true、完成下端と独立排除上端と宣言幅を
確認した場合だけ、その離散モデルの**区間**を報告する。中点を点値にしない。
探索終了0やrequirement_metは正確なFs、連続体崩壊、設計安全認定ではない。
停止は`stop-srm-search-job`、復旧は`prepare-srm-search-recovery OLD NEW`で元を保持する。

実断面SRMには実測形状/材料/施工/出典が別に必要で、今回未実行。
既存準備器を使う場合も、Griffithsを実断面と呼ばない：

```powershell
$real='D:/soil-inputs/actual-section'
$bundle='runs/NEW-real-section-srm'
.venv\Scripts\python tools/prepare_real_section_srm.py prepare --solver $S --mesh "$real/mesh.json" --ground "$real/ground.json" --case "$real/case-ilut.json" --spec "$real/specification.json" --source "$real/断面出典.pdf" --section-id '出典の断面ID' --output $bundle
```

出典・実入力が未用意なら準備しない。詳細はdocs/EXACT_FS_REAL_SECTION_SRM.ja.md。
P6-02/P6、NC-01、A/B/C各残件は未達のまま。push後Windows CIは当該先端の結果だけを確認し、
Linux成功・旧CI・取消/時間切れから成功を推定しない。

## 分割証拠の再読込（新規求解ではない）

```powershell
.venv\Scripts\python "$B/package_evidence.py" --verify "$B/evidence-manifest.json" --scratch work/NEW-p602-evidence-verify
.venv\Scripts\python "$B/collect.py" work/NEW-p602-evidence-verify/raw work/NEW-p602-remeasured
.venv\Scripts\python "$B/validation_summary.py" work/NEW-p602-evidence-verify/raw work/NEW-p602-validation.json
```

全part・連結ZIP・全memberのbytes/SHAが一致し、再集計JSONが公開measurements/restarts/validationと
完全一致すること。raw内のLinux実行版はバイト認証用で、Windowsで実行しない。
Windowsでは上記の専用ビルドから新規求解・新規checkpointを発行する。
生証拠は元ファイルを保持するが、source snapshotの開発時刻/配置先は再ビルドのビット一致を保証しない。
過去の封印・参照値・失敗記録を新しい実行版SHAへ付け替えない。
