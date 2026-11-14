# A-02原因診断の再現とPC追加確認

PR #13はマージしない。2026-12-21は指定ラベル。PC実行はしていない。
元の受入基準・定義・許容値・拒否範囲を維持する。新しい空の出力先だけを使う。
GPD WIN Max 2、利用可能約48GBを前提とし、**768/6144要素の新規手順、40GB超、
1経路1日超の見込みの計算は計画しない**。その規模を必要とする確認は
**規模上限のため確認不可・未達**。許容値や合格条件は変えない。

## 準備と固定先端

リポジトリルートでPowerShell 7を開く。元の作業ツリーに変更があれば保存したまま別worktreeを使う。

```powershell
git status --short
git fetch --filter=blob:none origin codex/plastic-spatial-convergence
git -c core.autocrlf=false worktree add --detach ../NEW-a02cause-repro FETCH_HEAD
if ($LASTEXITCODE -ne 0) { throw '新規worktree作成失敗。既存の保存先を上書きしない' }
Set-Location ../NEW-a02cause-repro
git -c core.autocrlf=false sparse-checkout disable
if ($LASTEXITCODE -ne 0) { throw '元アーカイブの配置失敗' }
git rev-parse HEAD
git merge-base --is-ancestor 7b254b1cdbe9c288989631944cd698aa861f83ac HEAD
if ($LASTEXITCODE -ne 0) { throw '指定起点の子孫ではありません' }
py -3.13 -m venv .venv
if ($LASTEXITCODE -ne 0) { throw 'Python 3.13の仮想環境作成失敗' }
$Py = '.venv/Scripts/python.exe'
$Study = 'interim/2026-12-21-a02-joint-birth-cause-diagnosis'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
& $Py -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw '固定依存の配置失敗' }
cargo build --locked --release
if ($LASTEXITCODE -ne 0) { throw 'CLI build失敗' }
cargo build --locked --release --examples
if ($LASTEXITCODE -ne 0) { throw '観測example build失敗' }
```

新規worktreeをautocrlf=falseで作り、旧JSON/固定SHAのバイトを維持する。
元の作業ツリーと設定は保持する。全保存証拠の取得には別途数十GBのディスクと回線時間を見込む。
旧証拠を読む全件試験には、sparse除外した元アーカイブも必要。
必要な元blobを取得し、コミットのGit blob SHAとバイト照合する。旧snapshotや参照を再生成しない。
ビルドは通常数分〜30分、メモリ2〜8GBの計画枠（PC未測定）。診断は実行前にバイナリを固定コピーする。

## 今回と同じ小型行列

```powershell
& $Py "$Study/a02cause_driver.py" --output runs/NEW-a02cause-n2-4-8 --steps 2 4 8
if ($LASTEXITCODE -ne 0) { throw '診断未完了。failure.jsonと途中証拠を保持' }
```

既存81節点・24 Tet10・648材料GP・32接触GP、3つの同じ乾燥経路だけ。
本体は各Nで観測あり/なし、独立側は入力だけからの各Nの全体平衡、さらにN=4で
接触上限1/4とGL7の独立全体対照を実行する。N=4の実checkpoint12地点、改変18件、
CLI拒否9件を検査する。局所再生は求解数へ加算しない。
PC計画枠は10〜45分、メモリ2〜8GB、追加ディスク1〜4GB（PC未測定、24 Tet10のみ）。
プロセスの実測時間は各 `*-process.json`、未完了は `failure.json` に残る。
既定の各子プロセス時間予算は900秒（旧ヘルパー由来のnativeは300秒）。
予算超過をscope拒否・物理破壊・合格へ変換しない。

見る値は `measurements.json` の同格子全場差、共通粗座標の本体/独立それぞれの
刻み差、接触上限/体積積分対照、局所再生と元の保存仕事との差。
`projection-observations/diagnosis.json` で行列の階数/零行、節点再現誤差、厳密cone条件、
順序/丸め感度、生のtractionを確認する。cone境界対照の分類は観測であり、丸め許容値を追加しない。
`geometry-audit.json` の全GP ID/所有/形状/重み/法線、`equilibrium-audit.json` の
全自由節点残差・局所仕事恒等式を確認する。全Runの観測あり/なしの差は時間だけであることを確認する。

## 同じN4での接触分割の追加対照

```powershell
& $Py "$Study/a02cause_active_control.py" --run runs/NEW-a02cause-n2-4-8 --output runs/NEW-a02cause-active
if ($LASTEXITCODE -ne 0) { throw '接触上限対照未完了' }
& $Py "$Study/a02cause_partition_control.py" --run runs/NEW-a02cause-n2-4-8 --output runs/NEW-a02cause-min2 --minimum-parts 2
if ($LASTEXITCODE -ne 0) { throw '接触最小分割の独立平衡未完了' }
& $Py "$Study/a02cause_report.py" --run runs/NEW-a02cause-n2-4-8 --active runs/NEW-a02cause-active --partitions runs/NEW-a02cause-min2 --output runs/NEW-a02cause-report
if ($LASTEXITCODE -ne 0) { throw '集計失敗' }
```

上限0.0005 mと最小2分割を別々に診断する。各3本の入力のみの独立全体平衡で、
合計数分〜15分、メモリ2〜4GB、ディスク0.1〜0.5GBの計画枠（PC未測定）。
見る値は `active-contact-controls.json` / `partition-global-controls.json` の実分割数、
自由力残差、全場差、特に元定義の仕事と数値散逸。
上限変更だけが開口中の分割を増やした場合、閉鎖中の仕事への有効な摂動にはならない。
最小分割の対照は旧参照の1箇所だけを隔離して変え、既定1は旧関数をそのまま返す。
旧参照・生産側の物理・定義・設定は変更しない。
集計は新規求解ではない。配布アーカイブを `--extract runs/NEW-a02cause-extracted` で検証展開した場合も、
同じreportへ `--run runs/NEW-a02cause-extracted/formal`、
`--active runs/NEW-a02cause-extracted/formal/supplementary/active`、
`--partitions runs/NEW-a02cause-extracted/formal/supplementary/minimum2` を渡せる。

## 追加で有用な長系列（PCのみ、今回未実行）

同じ24 Tet10のN=16/32だけを追加候補とする。メッシュや支持を増やさず、
元の経路と初期力の定義を保つ。まず16を完了し、実時間とメモリが枠内の時だけ32へ進む。

```powershell
$env:A02_PROCESS_BUDGET_SECONDS='3600'  # 各子プロセス上限1時間。未完了は失敗記録
try {
  & $Py "$Study/a02cause_driver.py" --output runs/NEW-a02cause-n16 --steps 16
  if ($LASTEXITCODE -ne 0) { throw 'N16未完了。N32は実行せず証拠を保持' }
  & $Py "$Study/a02cause_driver.py" --output runs/NEW-a02cause-n32 --steps 32
  if ($LASTEXITCODE -ne 0) { throw 'N32未完了。証拠を保持' }
} finally { Remove-Item Env:A02_PROCESS_BUDGET_SECONDS -ErrorAction SilentlyContinue }
```

計画枠：N16は20〜90分、N32は1〜4時間、いずれもメモリ4〜12GB、ディスク2〜8GB。
PC未測定の保守的な予算であり保証値ではない。1物理経路の観測あり/なし/独立は
各1時間上限、失敗時は進まない。1日超・40GB超への自動延長はしない。
出力は上記の新フォルダー。各Nで `completion.json` の診断完了と `passed=false` を区別する。
N16/32単独指定ではN4専用の拒否/復旧/投影/GL7対照を重複実行しない。

共通 `i/2` 座標だけの比較コマンド（補間しない）は次のとおり。

```powershell
& $Py -c "import sys;from pathlib import Path;sys.path.insert(0,'interim/2026-12-21-a02-joint-birth-cause-diagnosis');from a02cause_binding import fields;from a02cause_metrics import compare,common;from a02cause_driver import write;roots=[Path('runs/NEW-a02cause-n16/n16'),Path('runs/NEW-a02cause-n32/n32')];result={name:{side:compare(common(list(fields.load_fields(roots[0]/name/(side+'-fields.jsonl'))),2),common(list(fields.load_fields(roots[1]/name/(side+'-fields.jsonl'))),2),fields.CATEGORIES) for side in ('native','independent')} for name in ('transfer_force','stress_free','stress_free-cycle')};write(Path('runs/NEW-a02cause-16-to-32.json'),result)"
if ($LASTEXITCODE -ne 0) { throw '座標/全場比較失敗' }
```

変位/応力/traction/塑性すべり/累積長/仕事/散逸/誕生参照/凍結履歴を全点で見る。
差が残る場合は刻み収束未達。新しい許容値、Richardson補正、外挿だけで合格を作らない。
N32を連続極限の参照としない。別のh系列の独立参照は未実装・未確認であり、
24 Tet10の時間/荷重格子対照を空間収束に読み替えない。
より大きな空間確認が上記上限を必要とする場合は「規模上限のため確認不可・未達」。
小型のままで可能な局所分割方針・枝切替・投影境界感度の診断はクラウド残件である。

## 生証拠の検証と全件試験

```powershell
& $Py "$Study/a02cause_evidence.py" --verify $Study
if ($LASTEXITCODE -ne 0) { throw '新しい封印証拠のSHA不一致' }
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release
cargo test --locked --release
cargo test --locked --release --examples
foreach ($group in @('python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence','python-c01-switch','python-c01-path','python-c01-cause')) {
  & $Py tools/run_ci_python_evidence_shards.py --group $group --output "runs/NEW-a02cause-validation/$group"
  Write-Output "$group exit_code=$LASTEXITCODE"
}
& $Py interim/2026-12-06-p6-01-status-reconciliation/audit_status.py --help
```

全件の所要は数十分〜数時間、メモリ2〜12GBと元アーカイブの展開用ディスクを見込む。
全群を実行し、全発見IDの和集合・重複/欠落/skip0を確認する。失敗しても他群を実行し、
成功/失敗/error/未実行/skip/ignore/filter除外を区別して保存する。
新しい文書末尾検査はuniversal newlineで読み、実CRLFを作る試験を含む。
既存の分割方式を使い、重い新計算をCIへ追加しない。
元のNC厳密差が出た場合も基準・コード・参照・封印を変更しない。
今回の結果はREADME/validation.jsonを参照。**push後Windows CIは未確定**。
