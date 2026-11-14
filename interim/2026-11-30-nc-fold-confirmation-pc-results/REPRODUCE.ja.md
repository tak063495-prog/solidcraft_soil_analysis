# 案1：fold 確定に不足する照合の追加PC手順

対象は Windows の既存 worktree `D:\Work\scsa_b069ec8`。PC から push しない。既存の PC 履歴・失敗出力・ZIP・checkpoint は保持する。保存先の日付は依頼ラベルで、今回の転記対象の実施日は 2026-11-12 JST。以下は**今後の実行手順**であり、実施済みとは数えない。

元の [b069ec8 の案1手順](../2026-11-17-nc-fold-stability/REPRODUCE.ja.md) と [702a19a の再実行手順](../2026-11-24-nc-pc-rerun-fix/REPRODUCE.ja.md) の既存オプションだけを使う。NC の物理・ソルバー・既定・許容値を変更しない。終了コード0はプロセスの正常終了であり、fold 確定・数値受入ではない。各失敗は保存し、後段に欠落値を代入しない。

## 1. コード・入力・出力先を固定する

計算コードの基準 HEAD は今回の作業起点 `e6e4e735495951661bcc7b1d72168c14866597b6`。転記PC結果の HEAD `702a19a7ce189f1173caed727ee2a0059298a607` と混同しない。履歴の前史は既存の `runs/NEW-nc-fold-PC/series/{128,256}/history` を使い、再求解しない。数値コードを変えていない今回の文書コミットを読む場合も、実際に計算したHEADを必ず保存する。

```powershell
Set-Location 'D:\Work\scsa_b069ec8'
git status --short
git diff --quiet
if ($LASTEXITCODE -ne 0) { throw '追跡ファイルの変更を先に保存してください' }
git diff --cached --quiet
if ($LASTEXITCODE -ne 0) { throw 'ステージ済み変更を先に保存してください' }
git fetch origin codex/plastic-spatial-convergence
if ($LASTEXITCODE -ne 0) { throw 'fetch失敗' }
git switch --detach e6e4e735495951661bcc7b1d72168c14866597b6
if ($LASTEXITCODE -ne 0) { throw 'HEAD固定失敗' }
if ((git rev-parse HEAD) -ne 'e6e4e735495951661bcc7b1d72168c14866597b6') { throw 'HEAD不一致' }
$Py = (Resolve-Path '.venv/Scripts/python.exe').Path
$TraceRun = 'runs/NEW-nc-fold-PC'
$TaskRun = 'runs/NEW-nc-fold-confirmation-PC' # 未存在の保存先。再実行時は新しい名前にする
if (Test-Path $TaskRun) { throw '新しい出力先を指定してください' }
New-Item -ItemType Directory -Path $TaskRun -ErrorAction Stop | Out-Null
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:RUST_TEST_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
git rev-parse HEAD | Out-File "$TaskRun/head.txt" -Encoding utf8
Get-Date -Format o | Out-File "$TaskRun/start.txt" -Encoding utf8
& $Py -c "import sys,numpy,scipy,platform; print(sys.version); print(numpy.__version__,scipy.__version__); print(platform.platform())" | Out-File "$TaskRun/runtime.txt" -Encoding utf8
if ($LASTEXITCODE -ne 0) { throw 'Python環境確認失敗' }
foreach ($denom in @(128,256)) {
  $History = "$TraceRun/series/$denom/history"
  if (-not (Test-Path "$History/mesh.json")) { throw "元meshを備えた既存履歴が必要: $History" }
  & $Py tools/nc_fixed_inputs.py --history $History --output "$TaskRun/fixed-input-$denom.json"
  if ($LASTEXITCODE -ne 0) { throw "固定入力不一致: $denom" }
}
& $Py tools/prepare_nc_fold_observer.py --output "$TaskRun/observer"
if ($LASTEXITCODE -ne 0) { throw '観測器ビルド失敗' }
$Probe = (Resolve-Path "$TaskRun/observer/target/release/examples/material_diagnostic.exe").Path
```

固定依存は既存 lockfile を使う。観測器は一度だけビルドし、以下の直接／再開計算で同じバイナリを使う。入力検査の合格条件は `fixed_inputs_match=true` と正規化後の実差分0。`raw_equal=false` だけでは不合格にしない。未知フィールド・物理条件・制御の実変更は元検査のまま拒否する。想定：既存履歴検査は数秒〜数分、観測器の初回ビルド2〜10分程度。

**128/256 は同一メッシュでの撤去初期刻みの分母**。`prepare_nc_fold_series.py` は `stages[2].initial_increment` 以外を変えない。今回の経路へ細かいメッシュを渡す既存オプションはなく、`FreshNative` は固定幾何の変更を拒否する。新しいメッシュ系列や一致閾値をこの手順で作らない。さらに細かい撤去初期刻み1/512等もFreshNativeの既存許可値にないため提案しない。局所ds細分化で初期刻み系列間の差まで解消したとは扱わない。

## 2. 同一前史で局所刻みをさらに半減する

既報の粗・細を同じ観測器で再取得し、5e-6 m を追加する。3本とも弧長合計0.00256 mで比較する。bracket 近傍だけを取り出して checkpoint を編集する操作は行わない。

```powershell
foreach ($denom in @(128,256)) {
  $History = "$TraceRun/series/$denom/history"
  foreach ($spec in @(@('coarse','0.00002',128),@('fine','0.00001',256),@('finer','0.000005',512))) {
    $PathOut = "$TaskRun/native$denom-$($spec[0])"
    & $Py tools/nc_fold_path.py --map native --fresh-native $History --binary $Probe --step-size $spec[1] --count $spec[2] --budget-seconds 1800 --output $PathOut
    if ($LASTEXITCODE -ne 0) { throw "fold_path異常終了: $PathOut （出力を保持）" }
    & $Py tools/analyze_nc_fold_path.py --path $PathOut --output "$PathOut-audit"
    if ($LASTEXITCODE -ne 0) { throw "analyze異常終了: $PathOut （出力を保持）" }
  }
}
```

想定：PC転記の128点14秒・256点28秒から、512点は各約1分、6本の経路計算は計約3〜5分＋解析・I/O。これは見積もりであり時間保証ではない。各呼出しの上限は既存1800秒。点数制限や時間制限で停止したことを必ず区別する。

見るファイルは `report.json`、`*-audit/fold-results.json`、`*-audit/low-stress.json`、`controls.json`、`original-load.json`、`material/initial-replay.json`、`seed-check.json`、`events.jsonl`、`samples.jsonl`、`gp-tangents.jsonl`、`linear-systems.jsonl`、`initial.json`、`restart.json`、`evidence-sha256.json`。`initial-replay.json` は fresh-native の初期場再生の記録。出力で実際の配置を確認し、存在しない値を補わない。

局所 bracket の既存合格条件：3点以上の受理列、全点の力残差1e-8以下、拡張系／力ブロック線形残差各1e-10以下、向きを継承した `dλ/ds` が正→負、bracket 周辺で拡張系rank55、横断性 `>64*eps`、剛性rank53以上。経路の受理では力残差 `<1e-8`、弧長拘束 `<1e-12 m`、元の線形検査を満たす。全成分の比例荷重、GP接線組立、凍結GP／接触、source/model/checkpoint hash検査を維持する。

刻み細分化の確認はpeakだけで終えず、bracket幅、隣接点の滑らかさ、左右弱モード（符号任意性を考慮）、特異値の分離、全GP・接触履歴、荷重を比較する。**peak差の一致許容値は既存基準未定義**なので、新しい数値閾値で合格にしない。差が未収束、監査失敗、枝切替、入力・前史不一致なら未確定。全場は異なるdsの点を同一stepとして完全一致比較しない。

## 3. 点数上限で止まった経路の継続と直接／分割一致

次はfineの256点を実checkpointから256点だけ延長し、同じ起点からの直接512点と比較する。新しい出力先・同じソース・バイナリ・runtime・入力・dsを使う。転記対象の旧702a19a checkpointを別HEAD/バイナリで再開しない。

```powershell
foreach ($denom in @(128,256)) {
  $History = "$TraceRun/series/$denom/history"
  $First = "$TaskRun/native$denom-fine"
  $Second = "$TaskRun/native$denom-fine-next256"
  $Direct = "$TaskRun/native$denom-fine-direct512"
  if (-not (Test-Path "$First/restart.json")) { throw '実行時に発行されたcheckpointが必要' }
  & $Py tools/nc_fold_path.py --map native --fresh-native $History --binary $Probe --step-size 0.00001 --count 256 --budget-seconds 1800 --resume "$First/restart.json" --output $Second
  if ($LASTEXITCODE -ne 0) { throw "診断再開失敗: $denom" }
  & $Py tools/analyze_nc_fold_path.py --path $Second --output "$Second-audit"
  if ($LASTEXITCODE -ne 0) { throw "再開監査失敗: $denom" }
  & $Py tools/nc_fold_path.py --map native --fresh-native $History --binary $Probe --step-size 0.00001 --count 512 --budget-seconds 1800 --output $Direct
  if ($LASTEXITCODE -ne 0) { throw "直接計算失敗: $denom" }
  & $Py tools/analyze_nc_fold_path.py --path $Direct --output "$Direct-audit"
  if ($LASTEXITCODE -ne 0) { throw "直接監査失敗: $denom" }
  & $Py -c "import json,pathlib; root=pathlib.Path('$TaskRun'); read=lambda p:[json.loads(s) for s in p.read_text(encoding='utf-8').splitlines()]; names=('samples.jsonl','gp-tangents.jsonl','linear-systems.jsonl'); pairs=[(n,read(root/'native$denom-fine-direct512'/n),read(root/'native$denom-fine'/n)+read(root/'native$denom-fine-next256'/n)) for n in names]; assert len(pairs[0][1])==len(pairs[0][2])==512; assert all(a==b for n,a,b in pairs),{n:a==b for n,a,b in pairs}; assert json.loads((root/'native$denom-fine-direct512/restart.json').read_text())==json.loads((root/'native$denom-fine-next256/restart.json').read_text()); print('all fields, GP tangents, actual systems and checkpoint exact:',{n:len(a) for n,a,b in pairs})" | Out-File "$TaskRun/native$denom-restart-comparison.txt" -Encoding utf8
  if ($LASTEXITCODE -ne 0) { throw "直接/分割の全保存場不一致: $denom" }
}
```

想定：上記追加分は経路計算計約3分＋解析・I/O。合格は直接512点と256＋256点の全保存行完全一致。checkpointの本体・sidecarを同じフォルダで保持する。各実線形系、GP接線、最終checkpointのpayloadも比較し、パス・elapsed等の来歴は数値場と分けて記録する。上記コマンドは全保存行・GP接線・全実線形系と終端checkpointを完全一致で照合する。継続側だけに正→負がない場合の `undetermined` は、前半のbracketの消失を意味しない。直接512点の監査も読む。

追加継続が必要なら、同じdsで直前の実 `restart.json` を使い、新出力先に同じ256点ずつ追加する。既存上限は1呼出し `count<=4096`、`budget_seconds<=1800`、`0<ds<=1e-4`。上限内でも長時間の下降枝を自動追加しない。得た範囲・停止理由・最終λを記録し、未確定なら停止してよい。局所継続は本番の全撤去／restart-to-completionの代用ではない。

## 4. 入力由来の独立実装で前史400/800と材料内4/8を照合する

nativeの変位・応力・pc/e・保持荷重・受理係数は独立側へ渡さない。元の入力から独立に前史を作る。既存RK4は独立な静的材料写像の積分であり、案2の計算ではない。

```powershell
foreach ($steps in @(400,800)) {
  $Ref = "$TaskRun/reference$steps"
  $Arc = "$TaskRun/rate-arc$steps"
  & $Py tools/low_preload_nc_reference.py --output $Ref --steps $steps --method rk4 --budget-seconds 7200
  if ($LASTEXITCODE -ne 0) { throw "独立前史失敗: $steps" }
  & $Py interim/2026-11-09-nc-pseudo-arclength/independent_arc.py --reference $Ref --output $Arc --step-size 0.00002 --count 128
  if ($LASTEXITCODE -ne 0) { throw "独立arc異常終了: $steps" }
  foreach ($spec in @(@('coarse','0.00002',64),@('fine','0.00001',128))) {
    $PathOut = "$TaskRun/independent$steps-$($spec[0])"
    & $Py tools/nc_fold_path.py --map independent --reference $Ref --arc $Arc --step-size $spec[1] --count $spec[2] --material-steps 4 --budget-seconds 1800 --output $PathOut
    if ($LASTEXITCODE -ne 0) { throw "独立fold異常終了: $PathOut" }
    & $Py tools/analyze_nc_fold_path.py --path $PathOut --output "$PathOut-audit"
    if ($LASTEXITCODE -ne 0) { throw "独立監査失敗: $PathOut" }
  }
  $PathOut = "$TaskRun/independent$steps-material8"
  & $Py tools/nc_fold_path.py --map independent --reference $Ref --arc $Arc --step-size 0.00002 --count 64 --material-steps 8 --budget-seconds 1800 --output $PathOut
  if ($LASTEXITCODE -ne 0) { throw "材料内8の異常終了: $steps" }
  & $Py tools/analyze_nc_fold_path.py --path $PathOut --output "$PathOut-audit"
  if ($LASTEXITCODE -ne 0) { throw "材料内8の監査失敗: $steps" }
}
```

想定：前史は各5〜120分、局所は各1〜10分。既存 `independent_arc.py` の180秒制限や降伏ドリフト・線形残差による停止は成功にしない。rc=0でも、`reference.json`／`report.json`／`samples.json` の実際の受理起点・停止理由を確認する。空の列・未受理起点・閉鎖固着枝の範囲外なら後段を受入にせず、失敗出力を保持する。

合格条件は節2の残差・横断性・rank・比例荷重・凍結監査と、元の独立側検査を満たすこと。400→800、材料内4→8、局所ds半減、本体との各差を別表で示す。peakだけでなく全場・モード・荷重を比較する。native対独立の旧差約1.78546e-4は削除しない。独立一致のpeak許容値も**基準未定義**であり、数値だけから自動合格にしない。

## 4.1. 独立側でも新しい実checkpointを照合する

節4で同じ400/800前史から成功した材料内4の直接64点と、32＋32点を比較する。節3のnative再開一致を独立側の代用にしない。

```powershell
foreach ($steps in @(400,800)) {
  $Ref = "$TaskRun/reference$steps"
  $Arc = "$TaskRun/rate-arc$steps"
  $First = "$TaskRun/independent$steps-restart-first32"
  $Second = "$TaskRun/independent$steps-restart-second32"
  & $Py tools/nc_fold_path.py --map independent --reference $Ref --arc $Arc --step-size 0.00002 --count 32 --material-steps 4 --budget-seconds 1800 --output $First
  if ($LASTEXITCODE -ne 0) { throw "独立前半異常終了: $steps" }
  & $Py tools/nc_fold_path.py --map independent --reference $Ref --arc $Arc --step-size 0.00002 --count 32 --material-steps 4 --budget-seconds 1800 --resume "$First/restart.json" --output $Second
  if ($LASTEXITCODE -ne 0) { throw "独立再開異常終了: $steps" }
  foreach ($PathOut in @($First,$Second)) {
    & $Py tools/analyze_nc_fold_path.py --path $PathOut --output "$PathOut-audit"
    if ($LASTEXITCODE -ne 0) { throw "独立再開監査失敗: $PathOut" }
  }
  & $Py -c "import json,pathlib; root=pathlib.Path('$TaskRun'); read=lambda p:[json.loads(s) for s in p.read_text(encoding='utf-8').splitlines()]; names=('samples.jsonl','gp-tangents.jsonl','linear-systems.jsonl'); pairs=[(n,read(root/'independent$steps-coarse'/n),read(root/'independent$steps-restart-first32'/n)+read(root/'independent$steps-restart-second32'/n)) for n in names]; assert len(pairs[0][1])==len(pairs[0][2])==64; assert all(a==b for n,a,b in pairs),{n:a==b for n,a,b in pairs}; assert json.loads((root/'independent$steps-coarse/restart.json').read_text())==json.loads((root/'independent$steps-restart-second32/restart.json').read_text()); print('all independent fields, GP tangents, actual systems and checkpoint exact:',{n:len(a) for n,a,b in pairs})" | Out-File "$TaskRun/independent$steps-restart-comparison.txt" -Encoding utf8
  if ($LASTEXITCODE -ne 0) { throw "独立の直接/分割不一致: $steps" }
}
```

想定：各前史につき直接64点1本程度に相当する追加計算時間＋解析（概ね各1〜10分、実機差あり）。合格は64点の全保存行の完全一致と、元checkpoint／sidecar・全実線形系／GP接線の照合。材料内8の再開や本番の完了再開まで確認したとは数えない。原入力・runtime・ds・材料内分割・バイナリ／ソースを途中で変更しない。後半だけの`undetermined`と、直接経路全体のbracketを区別する。

## 5. 比較・PCからの持ち帰り

```powershell
$Reports = @(
  "$TaskRun/native128-coarse-audit/fold-results.json",
  "$TaskRun/native128-fine-audit/fold-results.json",
  "$TaskRun/native128-finer-audit/fold-results.json",
  "$TaskRun/native256-coarse-audit/fold-results.json",
  "$TaskRun/native256-fine-audit/fold-results.json",
  "$TaskRun/native256-finer-audit/fold-results.json",
  "$TaskRun/independent400-coarse-audit/fold-results.json",
  "$TaskRun/independent400-fine-audit/fold-results.json",
  "$TaskRun/independent400-material8-audit/fold-results.json",
  "$TaskRun/independent800-coarse-audit/fold-results.json",
  "$TaskRun/independent800-fine-audit/fold-results.json",
  "$TaskRun/independent800-material8-audit/fold-results.json"
)
foreach ($Report in $Reports) { if (-not (Test-Path $Report)) { throw "未実施/欠落を保持: $Report" } }
& $Py tools/compare_nc_fold_reports.py --reports $Reports --output "$TaskRun/comparison.json"
if ($LASTEXITCODE -ne 0) { throw '比較器異常終了' }
Get-Date -Format o | Out-File "$TaskRun/end.txt" -Encoding utf8
```

比較器の結論は固定で未確定、`passed=false`／`nc01_complete=false`。`static_limit_confirmed`、`history_converged`、`full_release_freeze_verified`、`restart_to_completion_verified` も診断器では固定falseであり、後から編集しない。`full_removal_reached` だけは終端λ=1に応じるが、それだけで本番受入にしない。局所foldの確定レビューとNC-01の完了は別で、NC-01は未達を保持する。

PCからはpushしない。フォルダ全体をローカルに保持し、原 `evidence-sha256.json` と全sidecar、observer情報、runtime、HEAD、実コマンド、日時、終了コード、各報告・比較・全場を一緒に持ち帰る。サイズの大きいmaterial呼出しログを含むため、必要なら既存のアーカイブ手段で分割する。クラウドへ渡した原ファイルがない場合は今回同様「PC出力の転記」として扱い、SHA・全場一致を検証済みにしない。

## 6. Linux全件検証の再現

```bash
source /workspace/soil-environment/activate.sh
export RUST_TEST_THREADS=1
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo test --locked --release --all-targets
cargo test --locked --release --doc
cargo build --locked --release
for group in python-core-0 python-core-1 python-mean-evidence python-spatial-evidence; do
  python tools/run_ci_python_tests.py --group "$group" --output "work/NEW-nc-fold-confirmation-validation/$group"
done
```

各群の終了コードを保存し、失敗しても残りを実行する。IDの和集合・重複・欠落とfailure/error/skip/ignored/filteredを実数で記録する。旧NC mean/spatialの計8件は成功・失敗・setup errorによる本体未実行を区別し、過去の未実行記録も保持する。Linuxでは `--allow-platform-float-differences` を使わない。CIの時間切れ・取消は成功へ読み替えない。

## 7. 今回のLinuxの短い経路計算を再現する

PCのZIPを復元する手順ではない。既存の封印アーカイブの1/128前史を、全メンバーのSHA照合後に読み出し、既存 `--evidence` adapter で新規局所経路だけを計算する。前史の再求解はしない。bracket近傍を含む粗32点／細64点、弧長合計0.00064 mに限定する。PCの128点／256点・0.00256 mの全区間再現とは区別する。PCは `--fresh-native` のWindows履歴を使っているため、peak差を純粋なOS差や同一前史の全場一致と解釈しない。128/256の両PC履歴をLinuxで再計算したとは数えない。

```bash
source /workspace/soil-environment/activate.sh
python interim/2026-11-17-nc-fold-stability/extract.py --output work/NEW-linux-fold-original
python tools/prepare_nc_fold_observer.py --output work/NEW-linux-fold-observer
python tools/nc_fold_path.py --map native --evidence work/NEW-linux-fold-original/evidence/prior --binary work/NEW-linux-fold-observer/target/release/examples/material_diagnostic --step-size 0.00002 --count 32 --budget-seconds 60 --output work/NEW-linux-fold-native128-coarse
python tools/analyze_nc_fold_path.py --path work/NEW-linux-fold-native128-coarse --output work/NEW-linux-fold-native128-coarse-audit
python tools/nc_fold_path.py --map native --evidence work/NEW-linux-fold-original/evidence/prior --binary work/NEW-linux-fold-observer/target/release/examples/material_diagnostic --step-size 0.00001 --count 64 --budget-seconds 60 --output work/NEW-linux-fold-native128-fine
python tools/analyze_nc_fold_path.py --path work/NEW-linux-fold-native128-fine --output work/NEW-linux-fold-native128-fine-audit
```

各コマンドの終了コードを個別に保存し、異常終了した経路の後段へ架空の結果を渡さない。アーカイブ取得・検証と観測器ビルドは経路計算時間と分ける。各経路は60秒上限で、時間切れなら未確定のまま保存する。見る値と局所bracketの合否条件は節2と同じ。実際のコマンド・runtime・開始時刻・停止理由・PC差は `linux-observations.json`、新しい全場・実行列・接線・checkpointの原出力は `linux-fold-path-evidence.zip` に保存する。PC転記JSONは変更しない。
