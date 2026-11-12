# PC-4：検証済みPC-3の400系を参照して§4の残りから進める

対象はPR #13の `codex/plastic-spatial-convergence`。**PCからpushしない。PRをマージしない。NC-01・P6は未達のまま。** フォルダ名は依頼ラベルであり実行日ではない。本書は今後のPC手順で、実施済みの記録ではない。

## 1. 新しいworktree・出力先を固定する

既存PC-3（`8342762`）のworktreeと出力を保持する。新worktreeを作り、既存PC-3を更新・移動・削除しない。通常のGit/Rust/Pythonが使えるPowerShellから実行する。

```powershell
Set-Location 'D:\Work\scsa_nc_arc_PC3'
git fetch origin codex/plastic-spatial-convergence
if ($LASTEXITCODE -ne 0) { throw 'fetch失敗' }
$RunHead = (git rev-parse FETCH_HEAD).Trim()
$RunHead
# チャット最終報告のpush済みSHAと一致することを確認する。
git worktree add --detach 'D:\Work\scsa_nc_arc_PC4' $RunHead
if ($LASTEXITCODE -ne 0) { throw '新worktree作成失敗。既存先を上書きしない' }
Set-Location 'D:\Work\scsa_nc_arc_PC4'
if ((git rev-parse HEAD).Trim() -ne $RunHead) { throw 'HEAD不一致' }
git diff --exit-code
if ($LASTEXITCODE -ne 0) { throw '追跡ファイル変更あり' }
git diff --cached --exit-code
if ($LASTEXITCODE -ne 0) { throw 'ステージ済み変更あり' }

$Pc3 = 'D:\Work\scsa_nc_arc_PC3\runs\NEW-nc-fold-confirmation-PC-3'
$Pc2 = 'D:\Work\scsa_b069ec8\runs\NEW-nc-fold-confirmation-PC-2'
$TaskRun = 'runs\NEW-nc-fold-confirmation-PC-4'
if (Test-Path $TaskRun) { throw '既存PC-4も失敗記録として保持。別の未存在名を使う' }
# PC-3で実際に使った既存Pythonを読み取る。新規環境に置き換えない。
$Runtime = Get-Content -LiteralPath "$Pc3\runtime-command.json" -Raw -Encoding utf8 | ConvertFrom-Json
$Py = (Resolve-Path -LiteralPath $Runtime.executable).Path
```

`$Pc2` は実際のPC-2出力場所に合わせる。新規 `.venv` が必要な場合のみ、Python 3.13と `requirements-reference.lock.txt` の固定依存で作成し、`$Py` をその実行ファイルへ変更する。新runtimeはPC-4に記録され、PC-3と同じ環境だったという主張にはしない。

疎配置の場合、新規2テスト、`tools/nc_pc_gate.py`、新スクリプト、旧arcディレクトリのPython、tools Python、examples JSON、`reference400-fixture.zip`とそのSHA JSONが必要。検査が要求する元ファイルを固定HEADから復元し、省略・別版・ダミー値で通さない。バイナリZIPはPowerShellのテキストpipeで保存しない。浅いcloneで `8342762bd714e33523ce5fef772572699fc21881` が見つからない場合は、その完全SHAを `git fetch origin <完全SHA>` で取得する。検証不能なら停止する。

## 2. シェルの旧件数と新ゲートを確認する（小さい検査）

Windows PowerShell 5.1とPowerShell 7が両方ある場合、以下をそれぞれ実行する。いずれも合成128行の検査で数値経路を計算しない。出力先は未使用名を指定する。

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File `
  './interim/2026-12-08-pc3-rate-arc-gate-fix/test-rate-arc-gate.ps1' `
  -Python $Py -Output 'runs\PC-4-gate-check-ps51'
if ($LASTEXITCODE -ne 0) { throw '5.1ゲート検査失敗' }
pwsh -NoProfile -ExecutionPolicy Bypass -File `
  './interim/2026-12-08-pc3-rate-arc-gate-fix/test-rate-arc-gate.ps1' `
  -Python $Py -Output 'runs\PC-4-gate-check-ps7'
if ($LASTEXITCODE -ne 0) { throw '7ゲート検査失敗' }
```

`shell-result.json` の `old_wrapped_count`、`direct_count`、10ケースのrc・判定・失敗条件を保存する。旧式の件数は5.1で1、7で128、直接代入した配列の件数は両方128を検査する。新ゲートの10ケースは両方同じ結果になる。7未導入なら5.1だけでもPC-4の共通Pythonゲートは使えるが、PC実機の両版検証を実施済みと記録しない。

PC-3生ファイルについて旧式の観測だけを保存する場合は、新しい `runs` ファイルへ出力する。元PC-3には書き込まない。

```powershell
function Read-Json([string]$Path) {
  Get-Content -LiteralPath $Path -Raw -Encoding utf8 | ConvertFrom-Json
}
$Wrapped = @(Read-Json "$Pc3\rate-arc400\samples.json")
$Direct = Read-Json "$Pc3\rate-arc400\samples.json"
$Rate = Read-Json "$Pc3\rate-arc400\report.json"
[ordered]@{
  shell=$PSVersionTable.PSVersion.ToString(); wrapped_count=$Wrapped.Count
  direct_count=$Direct.Count; equilibrium_samples=$Rate.equilibrium_samples
  stop_reason=$Rate.stop_reason; passed=$Rate.passed
  maximum_bordered_true_linear_residual=$Rate.maximum_bordered_true_linear_residual
} | ConvertTo-Json -Depth 8
# 上の表示を持ち帰る。report/samples/失敗JSONの原ファイルも保持する。
```

## 3. PC-4の実行コマンド

```powershell
& './interim/2026-12-05-nc-independent-arc-residual/rerun-pc.ps1' `
  -Python $Py -Pc2 $Pc2 -ResumeFrom $Pc3 `
  -Output 'runs\NEW-nc-fold-confirmation-PC-4'
```

新出力先の存在時は開始を拒否する。PC-3の下を出力先にすることも拒否する。例外が出たらその場で停止し、出力全体を残す。既存PC-4への途中書込み再開・自動 `SkipCompleted` は行わない。本オプションは検証したPC-3の2出力だけを参照する。

実際の順序は以下のとおり。

1. 新しいHEAD・PowerShell版・日時、runtime、build、既存arc回帰6件、新ゲート15件をPC-4へ記録する。PC-2 native監査6本の元パスとSHAも保存する。
2. `pc3-reuse-validation.json` を作成する。元HEADは8342762の完全SHA、reference400/rate-arc400の元コマンドは同HEAD・rc=0・所定の400/RK4、128点/ds=2e-5を要求。400前史にexcavate_upper起点・元入力・unloadがあり、arcの128点・明示null・有限の元残差1e-10以下・boolean falseを検査する。全サンプルの形状と有限値・制御、元HEADと現ソース/入力の同一性も確認する。参照ファイルとPC-3失敗のSHAを保存する。
3. 検証を通った `reference400` と `rate-arc400` は**PC-3の元パスのまま読み取り専用参照**する。コピーや再求解でPC-4生成物に見せない。invalid/欠損/件数不足/型違い/ソース差なら `pc4-failure.json` に具体的な条件を保存して停止する。
4. §4の残り：PC-4でindependent400-coarse（ds=2e-5・64点・材料内4）、fine（ds=1e-5・128点・材料内4）、material8（ds=2e-5・64点・材料内8）と各監査を新規実行。
5. §4.1の400系：PC-4の新しい32点checkpointから後半32点を実行。直接64点と `samples.jsonl`・`gp-tangents.jsonl`・`linear-systems.jsonl` の全行、終端restartをJSON値で完全一致比較する。旧HEADの独立checkpointは使わない。
6. 800系：reference800、rate-arc800の128点と同じ厳格ゲート、independent800の3本と各監査、同じ§4.1をPC-4に新規実行する。
7. §5：PC-2 native監査6本とPC-4 independent監査6本を `compare_nc_fold_reports.py` で比較する。`comparison.json` の `passed=false`、`nc01_complete=false` を維持。PC-2元監査とPC-3参照元全ファイルのSHAが変わっていないことを最後に確認する。

検証できない場合は停止が既定。400系から再計算する場合は、PC-3/失敗PC-4とは**別の未存在出力先**を指定し、`-ResumeFrom` を付けず同スクリプトを実行する。元PC-3を修正して再利用検査に合わせない。

## 4. 持ち帰るもの・判定範囲

PC-4出力全体を保持する。HEAD/runtime/PowerShell版、コマンドJSONとrc、全ログ、ゲートJSON、元参照のSHA manifest、前史、全場、GP接線、実線形系、restart/sidecar/hash、監査、直接/分割比較、comparison、失敗JSONも含む。`pc4-failure.json` の `failure_stage` と `failed_conditions`、`gate_diagnostic` を確認する。PC-3の `pc3-failure.json` は過去の停止としてそのまま持ち帰る。

400系はPC-3由来、800系と独立補正・新restartはPC-4由来、nativeはPC-2由来として区別する。前史400/800は積分分割数、native128/256は撤去初期刻みの分母であり、メッシュ細分ではない。後半32点だけのundeterminedと直接経路の局所bracketも区別する。

既存の力1e-8、真の線形1e-10、弧長1e-12 m等の基準・元荷重成分制約・材料写像・保存hash検査は変更しない。peak一致の許容値は未定義のまま。局所bracketや全手順の終了を、全撤去・物理限界・NC-01/P6完了に読み替えない。全解放・全解放後凍結・本番完了までの再開は別の残件である。

想定時間はビルド数分、800前史は既存7200秒上限、rate arcは既存180秒上限、各独立補正は既存1800秒上限。既存上限を増やさず、時間切れは成功にしない。クラウドでこの一連の大規模計算は行わない。

## 5. クラウドでの小さい回帰検査

```bash
source /workspace/soil-environment/activate.sh
python -m unittest discover -s tests -p test_nc_pc_gate.py -v
python -m unittest discover -s tests -p test_nc_pc_gate_schema.py -v
python -m unittest discover -s tests -p test_nc_independent_arc_residual.py -v
python -m unittest discover -s tests -p test_nc_fold_stability.py -v
git diff --check
```

PowerShellがない場合も最初のモジュールのPython14件は実行でき、シェル連携1件だけskipとなる。実生成・再利用モジュールには固定依存と元fixtureが必要で、PowerShellなしでは実スクリプト連携1件だけskipとなる。疎配置ではfold既存試験が使う `interim/2026-11-08-low-preload-nc-followup/native-evidence.zip` も固定HEADから復元する。新しい参照値や代替アーカイブを作らない。
