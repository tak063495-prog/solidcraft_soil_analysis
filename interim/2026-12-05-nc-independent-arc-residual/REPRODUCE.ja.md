# 独立 arc 残差修正後の PC-3 再実行

対象：PR #13、`codex/plastic-spatial-convergence`。起点は `e278dbf17c771e87098b6180e9c65d9e4a6ad2e7`。本書は**今後の Windows 実行手順**。PC-2（2026-11-12 13:03 JST 終了）の転記や今回の Linux 測定を PC-3 の実施結果にしない。**PCからpushしない。PRをマージしない。NC-01は未達のまま。**

## 1. HEAD・環境・保存先

新規 worktree を推奨する。通常の Git/Rust/Python が利用できる PowerShell で実行する。既存 PC-2 出力と前史は移動・削除・上書きしない。

```powershell
Set-Location 'D:\Work\scsa_b069ec8'
git fetch origin codex/plastic-spatial-convergence
if ($LASTEXITCODE -ne 0) { throw 'fetch失敗' }
$RunHead = (git rev-parse FETCH_HEAD).Trim()
# 最終報告のpush先端と同じSHAであることをここで照合する。
$RunHead
# 実装固定コミットを含むPR先端だけを使用する。
git merge-base --is-ancestor 10b81658985ba6bbd8ddecb647e991859ed99d52 $RunHead
if ($LASTEXITCODE -ne 0) { throw '今回の実装修正を含まないHEAD' }
git worktree add --detach 'D:\Work\scsa_nc_arc_PC3' $RunHead
if ($LASTEXITCODE -ne 0) { throw '新worktree作成失敗。既存先は上書きしない' }
Set-Location 'D:\Work\scsa_nc_arc_PC3'
if ((git rev-parse HEAD).Trim() -ne $RunHead) { throw 'HEAD不一致' }
git status --short
```

`$RunHead` は fetch 直後に得た不変のコミット値。実装コミットは `10b81658985ba6bbd8ddecb647e991859ed99d52`。本書を追加した文書・測定コミットを含む最終 push 先端の SHA は最終報告と照合する。計算前にローカルの追跡ファイル変更がないことを確認する。

既存 `D:\Work\scsa_b069ec8` を更新する場合は、新 worktree 作成の代わりに下記を実行する。未追跡の `runs/` や `.venv/` を消さない。

```powershell
Set-Location 'D:\Work\scsa_b069ec8'
git diff --quiet
if ($LASTEXITCODE -ne 0) { throw '追跡ファイル変更を先に保存' }
git diff --cached --quiet
if ($LASTEXITCODE -ne 0) { throw 'ステージ済み変更を先に保存' }
git switch --detach $RunHead
if ($LASTEXITCODE -ne 0) { throw 'HEAD固定失敗' }
if ((git rev-parse HEAD).Trim() -ne $RunHead) { throw 'HEAD不一致' }
```

疎配置のworktreeで新規テスト用fixtureだけが省略されている場合は、固定HEADからそのファイルを復元する。バイナリZIPをPowerShellのテキストpipeで保存しない。

```powershell
$Fixture = 'interim/2026-12-05-nc-independent-arc-residual/reference400-fixture.zip'
if (-not (Test-Path $Fixture)) {
    git restore --source=$RunHead --ignore-skip-worktree-bits -- $Fixture
    if ($LASTEXITCODE -ne 0) { throw '固定fixtureの取得失敗' }
}
```

既存の Python 環境を使える。新 worktree で仮想環境を作る場合だけ以下を行う（既存 `.venv` がある場合は作り直さない）。依存は `requirements-reference.lock.txt` の固定版を使う。

```powershell
# 新worktreeのみ。既存Pythonを直接使う場合、この2行は不要。
py -3.13 -m venv .venv
& '.venv/Scripts/python.exe' -m pip install -r requirements-reference.lock.txt
# 失敗時は先へ進まない。
if ($LASTEXITCODE -ne 0) { throw 'Python環境準備失敗' }
$Py = (Resolve-Path '.venv/Scripts/python.exe').Path
# 既存環境を直接使う場合の例：
# $Py = 'D:\Work\scsa_b069ec8\.venv\Scripts\python.exe'
$Pc2 = 'D:\Work\scsa_b069ec8\runs\NEW-nc-fold-confirmation-PC-2'
$TaskRun = 'runs/NEW-nc-fold-confirmation-PC-3'
if (Test-Path $TaskRun) { throw '別の未存在名を使う。PC-3失敗も消さない' }
```

## 2. §4・材料内4/8・§4.1・§5を実行する

このフォルダの [`rerun-pc.ps1`](rerun-pc.ps1) は下記の全工程を順に実行する。出力ディレクトリを新規作成し、HEAD、runtime、各コマンドの引数・開始終了時刻・rc、ログ、失敗報告を保存する。Rust は `cargo build --locked --release`、回帰は追加した Python 6試験を実行する。通常の JSON 出力の `passed=false` は保持する。

```powershell
& './interim/2026-12-05-nc-independent-arc-residual/rerun-pc.ps1' `
    -Python $Py -Output $TaskRun -Pc2 $Pc2
# 例外・rc非0・pc3-failure.jsonがあれば途中失敗。全出力をそのまま保持する。
```

実行する数値コマンドは旧手順のまま。以下はスクリプト内部の §4 の具体例で、スクリプト実行後に同じ出力へもう一度実行しない。

```powershell
# $Steps は 400 と 800 の両方
$Steps = 400
$Ref = "$TaskRun/reference$Steps"
$Arc = "$TaskRun/rate-arc$Steps"
& $Py tools/low_preload_nc_reference.py --output $Ref --steps $Steps --method rk4 --budget-seconds 7200
& $Py interim/2026-11-09-nc-pseudo-arclength/independent_arc.py --reference $Ref --output $Arc --step-size 0.00002 --count 128
# coarse: ds=2e-5,64点,材料内4。fine: ds=1e-5,128点,材料内4。
# material8: ds=2e-5,64点,材料内8。400/800それぞれ3本。
& $Py tools/nc_fold_path.py --map independent --reference $Ref --arc $Arc --step-size 0.00002 --count 64 --material-steps 4 --budget-seconds 1800 --output "$TaskRun/independent$Steps-coarse"
& $Py tools/analyze_nc_fold_path.py --path "$TaskRun/independent$Steps-coarse" --output "$TaskRun/independent$Steps-coarse-audit"
```

各数値呼出しの直後に rc を確認する。スクリプトは非0なら停止する。前史の rc=0 は受入を意味せず、旧剛性系残差で停止した `reference.json/failure_reason` と `last_accepted.name=excavate_upper` を確認する。`rate-arc*/report.json` は128点、`samples.json` も128点、`stop_reason=null`、元の拡張系残差 `<=1e-10` を必要とする。0点・途中停止・180秒上限・残差超過は後段へ進めない。

§4.1 は各400/800について、同じ新規入力前史・同じ arc・材料内4・ds=2e-5で32点＋32点を新しい実 `restart.json` から実行する。直接64点の `independent*-coarse` に対し、`samples.jsonl`、`gp-tangents.jsonl`、`linear-systems.jsonl` の全行と終端 `restart.json` を Python の JSON 値で**完全一致**比較し、`independent*-restart-comparison.log` に保存する。checkpoint／sidecar のhash検査も元のまま実行する。旧HEADのcheckpointを新ソースで再開しない。後半のみの `undetermined` と、直接経路の局所bracketを区別する。

§5 は PC-2 の native 6本を読み取り専用で参照し、PC-3 の independent 6本と一緒に `compare_nc_fold_reports.py --reports ... --output .../comparison.json` に渡す。PC-2 の native を PC-3 で再求解したとは数えない。`pc2-native-source.json` に元ファイルのパス・SHAを保存し、前後で一致を検査する。native のソースbindingと今回の独立側bindingは別のまま比較表へ保存する。PC-2 との同一HEAD再開や全場完全一致を主張しない。

想定時間：初回ビルド2〜10分程度、前史400/800は各5〜120分（既存上限7200秒）、rate arc各数十秒〜180秒、補正経路は各1〜10分（各上限1800秒）。全体で数十分〜数時間を見込む。実機・BLAS・I/O依存の見積もりであり、上限を増やしたり、時間切れを成功扱いにしたりしない。今回 Linux の前史400は約18.3秒・arc8点は約5秒だったが、PCの全128点／全工程の時間保証には使わない。

## 3. 見る値と判定

- `rate-arc*/linear-systems.json`：元の A,b、正規化前の実解、`true_relative_residual`、`solver.unscaled_direct_relative_residual`、平衡化後と4補正すべての残差、尺度・pivot・補正精度。**同一保存行列での前後**と、経路全体の前後最大値を区別する。元の `norm(A@x-b)/max(norm(b),1e-30)` で再計測する。Windowsの `longdouble` が倍精度と同じ場合は Python 3.12以降では `correction_precision=sumprod_extended`（拡張精度内積）、未対応Pythonでは `decimal50`。これは補正右辺だけの高精度化で、受入ノルムは変えない。
- `report.json`、`failure.json`、`failure.log`、各 `*-command.json` と `*.log`：0点、欠損、予算超過、ドリフト、rc、停止理由。今回からrate arcの数値停止はrc=1。空 `samples.json=[]` も必ず失敗。後段の欠損・空列は `failure_stage=independent_arc_seed` として保存し、FileNotFoundErrorで終わらせない。
- `*-audit/fold-results.json` と `low-stress.json`：受理点3点以上、力残差1e-8以下、拡張系／力ブロック線形残差各1e-10以下、向きを継承する `dλ/ds` の正→負、拡張系rank55、横断性 `>64*eps`、剛性rank53以上。実経路受理は力 `<1e-8`、弧長 `<1e-12 m`。元の比例荷重・GP接線・凍結・hash検査を維持する。
- 前史400→800、材料内4→8、ds半減はそれぞれ別の感度。peak λ、bracket幅、全節点変位・全GP材料場、左右弱モード（符号任意性を考慮）、特異値の分離、全成分荷重、接触・除去材料の凍結をレビューする。`comparison.json` は表とpeak差を作るだけで、全場収束を自動判定しない。独立対nativeの旧peak差約1.78546e-4も保持する。

**peak位置の一致基準は既存未定義**。新しい許容差で合格にしない。`comparison.json` は `conclusion_ja=未確定`、`passed=false`、`nc01_complete=false`。全撤去、全解放後凍結、本番の完了までの再開は今回の局所診断と別の未達項目。8点のLinux計算を局所foldの確認として使わない。

PC-2転記JSONは編集しない。PC-3のフォルダ全体・実際のHEAD/runtime/コマンド/rc/日時、前史、全場、全線形系、GP接線、全sidecar・hash、PC-2参照元を保持して持ち帰る。失敗出力も削除しない。原ファイルが渡せない場合は「PC出力の転記」と明記する。

## 4. Linux小規模再計算と全件検証

`128` はこの節では**独立前史の積分分割数**。メッシュは固定54節点・12 Tet10である。PCのnative系列にある128/256は撤去初期刻みの分母で、どちらもメッシュ細分ではない。

```bash
source /workspace/soil-environment/activate.sh
export RUST_TEST_THREADS=1
mkdir -p work/NEW-nc-independent-arc-residual
python tools/low_preload_nc_reference.py --output work/NEW-nc-independent-arc-residual/reference400 --steps 400 --method rk4 --budget-seconds 60
python interim/2026-11-09-nc-pseudo-arclength/independent_arc.py --reference work/NEW-nc-independent-arc-residual/reference400 --output work/NEW-nc-independent-arc-residual/after400 --step-size 0.00002 --count 8
python -m unittest discover -s tests -p test_nc_independent_arc_residual.py -v
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo test --locked --release --all-targets
cargo test --locked --release --doc
cargo build --locked --release
for group in python-core-0 python-core-1 python-mean-evidence python-spatial-evidence; do
  python tools/run_ci_python_tests.py --group "$group" --output "work/NEW-nc-independent-arc-residual/validation/$group"
done
```

疎配置でアーカイブを省略している場合は、全件試験の前に次の元ファイルも取得する。SHA検査は元の試験が行う。通常の完全checkoutでは不要。今回の初回Python coreの6エラー（実試験5件＋1クラス準備失敗）は旧wet-JOINTアーカイブと新規reference400-fixtureの配置不足で、元blob復元後に影響した4モジュール全21件を再試験した。初回のログ・rc・未実行6件と、復元後の実行を別記し、失敗を消さない。

```bash
# core.sparseCheckout=true の場合だけ実行
git sparse-checkout add \
  '/interim/2026-12-05-nc-independent-arc-residual/*' \
  '/interim/2026-11-08-low-preload-nc-followup/native-evidence.zip' \
  '/interim/2026-11-10-nc-mean-stress-release/*.part*' \
  '/interim/2026-11-10-nc-spatial-release/*.part*' \
  '/interim/2026-11-08-wet-joint-cut-volume/evidence.zip' \
  '/interim/2026-11-08-wet-joint-separated-convergence/evidence.zip' \
  '/interim/2026-11-10-wet-joint-small-refinement/evidence.zip' \
  '/interim/2026-11-15-wet-joint-convergence/evidence.zip.part*' \
  '/interim/2026-11-18-wet-joint-step-convergence/evidence.zip.part*'
```

各コマンドのrcを別々に保存し、失敗しても無関係な群は全件実行する。PythonのID和集合・重複・欠落、failure/error/skip/ignored/filtered、setup errorで本体未実行の件数を区別する。旧mean/spatialは必要な封印アーカイブをすべて取得してから実行し、Linuxで `--allow-platform-float-differences` を使わない。既知の2群8件が厳密検査で未実行なら実数をそのまま保存する。旧記録・snapshot/hash/参照値・試験を変更しない。

前後比較を含む小規模実験を新規出力へまとめて実行する場合は、次を使う（上の単独計算と重複して数えない）。旧版コードは指定起点のGit blobを読み取り、NumPyの求解結果を観測するだけで旧ファイルを変更しない。各前史60秒上限、arcは128前史で16点、400前史で修正前後各8点に限定する。

```bash
python interim/2026-12-05-nc-independent-arc-residual/run-small.py --output work/NEW-nc-independent-arc-before-after
python interim/2026-12-05-nc-independent-arc-residual/measure.py --work work/NEW-nc-independent-arc-before-after --output work/NEW-nc-independent-arc-before-after/measurements.json
```

保存した `bordered-regression.json` は修正前の短い400分割で実際に解いた最大残差の系。`solve_bordered` で元の系を再評価できる。`reference400-fixture.zip` は新規Linux前史の5入力／状態ファイルで、旧封印証拠の置換ではない。`cloud-evidence.zip` の原出力を新しい `work/` へ展開すれば、`measure.py --work <展開先> --output <新JSON>` で測定を再生成できる（保存版は上書きしない）。この計測器自体は物理経路を再求解しない。
