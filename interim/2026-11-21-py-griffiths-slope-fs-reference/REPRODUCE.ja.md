# Windows PCでの再現・未達の検証手順

この保存先の日付は依頼のラベル。実行確認日は2026-11-12（日本時間）。下のPC向け長系列、元661点の新規h実求解、実断面SRMはクラウドで未実行。計算時間はPCで未実測の目安であり、終了保証ではない。新しい空ディレクトリを使い、途中停止・失敗・古い出力を上書きしない。

## 1. 独立参照の環境と入力照合

Python 3.13、NumPy 2.5.4、SciPy 1.18.1は既存lockに従う。PowerShellをリポジトリのルートで開く。

```powershell
git fetch origin codex/py-griffiths-slope-fs-reference
git switch codex/py-griffiths-slope-fs-reference
git merge-base --is-ancestor 7ae0a2a669e0e1df93f3dd9791bbeda056334c3e HEAD
if ($LASTEXITCODE -ne 0) { throw '指定起点を確認できない' }
py -3.13 -m venv .venv-fs
$pyFs = (Resolve-Path .venv-fs/Scripts/python.exe).Path
& $pyFs -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw '依存の準備が未完了' }
$env:OPENBLAS_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$env:PYTHONUTF8 = '1'
$taskFs = 'interim/2026-11-21-py-griffiths-slope-fs-reference'
$runFs = 'runs/py-griffiths-fs-pc-01'
if (Test-Path $runFs) { throw '既存結果を上書きしない。別の出力先を指定する' }
New-Item -ItemType Directory $runFs | Out-Null
$sourceFs = 'https://raw.githubusercontent.com/tak063495-prog/excel-vba-soil-fem/440d246d9e32f5892c79539d4a7b797ad29442ee/docs/inconsistent-evidence24/fixtures.json'
Invoke-WebRequest $sourceFs -OutFile "$runFs/fixtures.json"
& $pyFs "$taskFs/audit_fixture.py" "$runFs/fixtures.json" --output "$runFs/fixture-audit.json"
if ($LASTEXITCODE -ne 0) { throw '出典/材料/幾何/拘束の照合が不一致' }
```

準備目安5〜15分、fixture照合1秒程度。`fixture-audit.json` は入力の整合だけを意味し、斜面の受入ではない。case.jsonの固定内容を変えると比較器は拒否する。本実装はこの一つの乾燥均質再構成だけで、実断面や湿潤条件をJSONだけ変えて解く汎用ソルバではない。

## 2. 新規試験と小規模系列

```powershell
& $pyFs -m unittest discover -s $taskFs -p checks_reference.py -v 2>&1 | Tee-Object "$runFs/new-tests.log"
if ($LASTEXITCODE -ne 0) { throw '新規試験に失敗' }
& $pyFs "$taskFs/run_reference.py" --output "$runFs/small" --grids 3 5 9 --slices 16 32 64 --h-levels 4 6 8 --trial-seconds 2 --max-trials 8 --gravity-steps 8 --width 0.025 2>&1 | Tee-Object "$runFs/small.log"
if ($LASTEXITCODE -ne 0) { throw '参照処理が未完了。残った出力を保持する' }
```

新規試験は数秒、小規模系列は概ね2〜10分。Linuxでは既定域18条件に加え、局所域 `--bounds 6 11.5 29 32 2 4` の同じ3/5/9×16/32/64の18条件、CST 96/216/384要素だけを実行した。出力は `environment.json`（版・固定入力・コードSHA）、`lem.json`（全探索条件の最小・面パラメータ・root/引張/推力線診断）、`srm-series.json`、`srm-n*/raw.json`（全試行・全Newton残差・確定した最後の重力段階・変位/応力/塑性ひずみ）、`interchange.json`、`nodes.csv`、`plastic.csv`、`legacy-diagnosis.json`。プロセス終了0は出力処理の完了であり、物理の合格や正確なFsの認定ではない。

`computational_transition=[equilibrium_lower,unresolved_upper]` は計算境界。`certified_fs_interval=null` と `accepted=null` が必須。上端が得られなければnullのまま。widthへ達しても、trial_limitやtimeoutでも、破壊側の証明は得られない。Fs中点は生成しない。タイムアウト予算はNewton境界で確認するため厳密なOS打切り時間ではない。外部キャンセルの場合は完了したように出力を補修しない。

## 3. PCだけで格子・スライス・h・載荷を分ける

まず探索格子とスライス数を交差させる。17³までで目安30分〜4時間、33³以上は半日〜数日になり得る。標準出力をログへ保存し、各条件の最小面と未解決root数を確認する。

```powershell
& $pyFs "$taskFs/run_reference.py" --lem-only --output "$runFs/lem-long" --grids 5 9 17 --slices 64 128 256 2>&1 | Tee-Object "$runFs/lem-long.log"
if ($LASTEXITCODE -ne 0) { throw '探索未完了。合格として数えない' }
# 既定探索域と分けた局所域の診断。大域探索の置換や校正には使わない。
& $pyFs "$taskFs/run_reference.py" --lem-only --output "$runFs/lem-local" --grids 5 9 17 --slices 64 128 256 --bounds 6 11.5 29 32 2 4 2>&1 | Tee-Object "$runFs/lem-local.log"
if ($LASTEXITCODE -ne 0) { throw '局所域の探索未完了' }
```

`--bounds A0 A1 B0 B1 S0 S1` は進入x/退出x/矢高の領域。既定と局所域の結果を別々に残す。領域外・より一般の非円弧・未検出rootを除外したことにはならない。同一格子のスライス差、同一スライスの格子差、最小面の変化・境界到達・力/モーメント残差・引張・推力線を観測する。新しい最低次数や許容差は未定義、`accepted=null`。小さな観測差をそのまま厳密誤差・Fs区間にしない。

独立SRMはhだけを変える系列と、自重増分を変える系列を分ける。8/12/16は384/864/1536 CST、全体1,200〜19,200 Tet10の既存モデルとは別である。目安10〜60分（各試行上限60秒×最大10試行×3水準、反復境界の超過を除く計算予算は30分）。必ず実測時間も記録する。

```powershell
& $pyFs "$taskFs/run_reference.py" --srm-only --output "$runFs/srm-h" --h-levels 8 12 16 --trial-seconds 60 --max-trials 10 --gravity-steps 16 --width 0.01 2>&1 | Tee-Object "$runFs/srm-h.log"
if ($LASTEXITCODE -ne 0) { throw 'h系列未完了' }
foreach ($gFs in @(8, 16, 32)) {
    & $pyFs "$taskFs/run_reference.py" --srm-only --output "$runFs/srm-gravity-$gFs" --h-levels 8 --trial-seconds 60 --max-trials 10 --gravity-steps $gFs --width 0.01 2>&1 | Tee-Object "$runFs/srm-gravity-$gFs.log"
    if ($LASTEXITCODE -ne 0) { throw '自重増分系列未完了' }
}
```

平衡端と未確定端、各試行状態、CSTの過剛性/局所化、載荷経路、材料活動集合・Newton停止を分けて評価する。連続極限のFsを得たと宣言しない。より一般の要素や独立破壊証明の実装は別検討であり、未確認の式や不収束から上端を生成しない。

## 4. 将来の本体出力を比較

最初は自分自身の交換形式検査だけを行える。これは本体実比較ではない。

```powershell
$refFs = "$runFs/small/srm-n8/interchange.json"
& $pyFs "$taskFs/compare.py" $refFs --reference $refFs --lem-results "$runFs/small/lem.json" --output "$runFs/interchange-self-diagnostic.json"
if ($LASTEXITCODE -ne 0) { throw '交換形式の配線が不整合' }
# 将来の本体exportから明示的に契約へ変換・監査したmanifestを指定する。
$candidateFs = 'D:/計算記録/Griffiths本体/independent-comparison-manifest.json'
& $pyFs "$taskFs/compare.py" $candidateFs --reference $refFs --lem-results "$runFs/small/lem.json" --output "$runFs/production-diagnostic.json"
if ($LASTEXITCODE -ne 0) { throw '比較拒否。単位/符号/状態/全試行/区間/同一点を確認する' }
```

交換JSONは例の `schema_version/model_id/case_sha256/units/conventions/reduction/accepted/fs_tolerance/trials/trial_count/factor_interval/endpoint_statuses/interval_kind/field_state` を全て用意する。全自重の平衡試行は `gravity_fraction=1,relative_residual,residual_tolerance<=1e-8` が必要。未収束端はそのまま `unresolved/iteration_limit/trial_limit/timeout/cancelled/constitutive_unresolved` とする。`factor_interval` は2端点が実試行に存在する正順の計算境界だけで、`interval_kind=computational_unresolved`。物理Fsの誤差許容は未定義、`fs_tolerance=null, accepted=null`。中点・点値は渡さない。

`field_state` は同じF、全自重完了、`reference=stress_free_gravity_turn_on`、`plastic_measure=engineering_tensor_components` と、以下の表ファイル名が必要。CSVヘッダーの順序・綴り・単位は固定。JSONの場合は `{"columns":[...],"rows":[[...],...]}`。

```text
node_id,x_m,y_m,ux_m,uy_m
element_id,x_m,y_m,plastic_exx,plastic_eyy,plastic_ezz,plastic_gamma_xy
```

本体3Dのx/zを独立2Dのx/yへ写すなら、その変換と全ID/物理点・拘束・全試行を別監査する。本体の独自JSONを拡張子だけ変えて渡さない。異なるF、未完了場、異なるメッシュの自動補間、圧縮正応力の黙示変換、単位タグの不一致、Fs中点化は拒否する。LEＭの包含/外れや二つの計算境界の交差は数値診断だけで、物理的合否は判定不能。新しい場の比較基準は未定義なので0終了でも `accepted=null`。

## 5. 元661節点・同条件のh系列はPCで別実求解

この節だけは**既存本体/既存比較器をPCで手動実行する別作業**。独立参照コードから呼ばない。入力・固定661点参照・許容値の説明は [既存Griffiths文書](../../docs/GRIFFITHS_2D_COMPARISON.ja.md)。F=1系列をF=1.3500000000000003の元変位へ比較しない。

```powershell
cargo build --locked --release
if ($LASTEXITCODE -ne 0) { throw '本体build未完了' }
$ref661 = 'tests/data/griffiths_2d_reference.json'
$match661 = "$runFs/matched-661"
& $pyFs tools/prepare_slope_h_series.py $match661 --reference-2d $ref661
if ($LASTEXITCODE -ne 0) { throw '元参照条件の準備未完了' }
foreach ($levelFs in @('xy1-ny1', 'xy2-ny1', 'xy4-ny1')) {
    $pathFs = "$match661/$levelFs"
    & target/release/soil-fem-3d.exe validate "$pathFs/mesh.json"
    if ($LASTEXITCODE -ne 0) { throw 'メッシュ検査拒否' }
    & target/release/soil-fem-3d.exe solve-stages "$pathFs/mesh.json" "$pathFs/ground.json" "$pathFs/case.json" "$pathFs/result.json"
    if ($LASTEXITCODE -ne 0) { throw '実求解未完了。Fs上端/破壊/合格としない' }
}
& $pyFs tools/griffiths_2d.py series "$runFs/matched-661-diagnostic.json" --reference-2d $ref661 --folder $match661 --solver target/release/soil-fem-3d.exe
```

1200/4800/19200 Tet10、2583/9963/39123節点、同じ奥行きと全奥行き拘束。想定は数時間〜数日以上、PCで未実測。出力は各入力・全Run・既存独立監査/VTU・全661点比較。元5%は全奥行きの非加重変位L2≤0.05、最大差≤0.10、奥行き差≤0.01を別々に確認する。既存比較器の終了0は全離散段階の5%だけ、1は実測不合格、2はデータ/監査拒否。連続側やP6の完了ではない。元15.582/15.886/15.593/15.672%不合格と連続側判定不能は別記録として保持する。

## 6. 実断面SRMは入力確認後のPC作業

実断面の形状・c/φ/γ・水・境界・施工履歴は本件では未確認。Griffithsや一様柱を実断面と呼ばない。既存の [実断面SRM手順](../../docs/EXACT_FS_REAL_SECTION_SRM.ja.md) を使い、出典付き入力が揃ったときだけ新規保存先へ準備する。

```powershell
$srmInputFs = 'D:/ユーザーが保有する実断面入力'
$srmBundleFs = 'D:/計算記録/実断面SRM-Fs参照01'
$srmSolverFs = (Resolve-Path target/release/soil-fem-3d.exe).Path
& $pyFs tools/prepare_real_section_srm.py prepare --solver $srmSolverFs --mesh "$srmInputFs/mesh.json" --ground "$srmInputFs/ground.json" --case "$srmInputFs/case.json" --spec "$srmInputFs/specification.json" --source "$srmInputFs/断面出典.pdf" --section-id '出典で特定した断面ID' --output $srmBundleFs
if ($LASTEXITCODE -ne 0) { throw '実断面出典または入力が未確認' }
$jobFs = "$srmBundleFs/job"
$jobSolverFs = "$jobFs/solver/soil-fem-3d.exe"
& $jobSolverFs verify-srm-search-job $jobFs
if ($LASTEXITCODE -ne 0) { throw '保存ジョブ拒否' }
& $jobSolverFs run-srm-search-job $jobFs
& $jobSolverFs srm-search-job-status $jobFs
& $pyFs tools/prepare_real_section_srm.py verify $srmBundleFs "$srmBundleFs/exact-fs-assessment.json"
```

規模・試行数が未確認なので想定時間は未定義（時間〜日単位になり得る）。準備だけは実行・合格ではない。既存仕様の宣言幅と下端平衡/独立排除上端が必要で、停止・timeout・試行上限では未達を保持する。中点を点値にしない。出典・保存時SHA・同じ保存solverを保全する。クラウドではこの節を実行していない。

## 7. 既存Python全件（CI収集を変更しない）

既存試験はRust CLIを必要とするため先にlocked release buildを行う。既存テスト証拠の6 GB超の分割ファイルも元Git blobどおり配置する。疎チェックアウトで欠けていれば取得し、skip/xfailで通さない。既存の収集はunittest、4群で全IDを一度ずつ割り当てる。

```powershell
foreach ($groupFs in @('python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence')) {
    & $pyFs tools/run_ci_python_tests.py --group $groupFs --output "$runFs/existing-$groupFs"
    if ($LASTEXITCODE -ne 0) { throw "既存全件未合格: $groupFs（ログを保持）" }
}
```

目安10〜60分、証拠取得時間は別。全群の発見集合・割当集合・実行件数・失敗/エラー/skip/expected failure/unexpected successを確認する。新規43件は別の明示起動で、既存件数へ紛れ込ませない。新先端Windows CIは、このPC手順やLinuxの成功から成功扱いにしない。

## 今回のLinux結果に残った制限

新規43件は全件成功。既存は543 IDを全てsuiteへ渡し、証拠不足の修復と該当5/6/4件の再実行後に、重複を除いて535件成功。NC平均応力/空間の封印証拠2群・計8件は `setUpClass` で厳密再計測不一致となり、試験本体は未実行。skip/expected failure/unexpected successは0だが、全件合格ではない。詳細・各ID・初回5エラー・再実行ログは [validation.json](results/validation/validation.json) と同じフォルダーに保存した。

PCで同じ準備エラーが出たら、そのまま失敗・未実行を記録する。元データや期待値を再ハッシュして差し替えたり、Linuxの厳密判定へplatform差許容を追加したり、失敗をskipへ変えて合格にしない。既存NCの原因解消は本件で実装しない。新規参照の成功と、既存全件未合格を別々に報告する。
