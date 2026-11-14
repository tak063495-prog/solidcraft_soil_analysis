# Windows再現手順（PCからpushしない）

今回の仮説の導出は手順2で不支持となった。手順3の材料点／円錐LP・SOCP／上界と、手順4の36条件は停止したため、ここで追加計算を求めるものではない。以下は保存した小部分の再測定と、未取得の1/256 fold場を将来レビューする場合の手順。

## 固定worktreeと実行環境

PowerShell 7またはWindows PowerShell。全パスは例なので未存在の作業先にする。`$Target`には**今回の最終報告に記載したpushコミットの40桁SHA**を入れる。ブランチ名だけで再現点を決めない。

```powershell
$Repo = 'D:/Work/solidcraft_soil_analysis'
$Target = '<今回の最終報告の40桁コミットSHA>'
$Work = 'D:/Work/scsa-nc-physical-limit'
git -C $Repo fetch origin codex/plastic-spatial-convergence
git -C $Repo worktree add --detach $Work $Target
Set-Location $Work
git rev-parse HEAD
git merge-base --is-ancestor e541fe1d29d05f7301968555845accd628e5e76f HEAD
if ($LASTEXITCODE -ne 0) { throw 'wrong base' }
py -3.13 -m venv .venv
$Py = (Resolve-Path '.venv/Scripts/python.exe').Path
& $Py -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw 'dependency install' }
$env:OPENBLAS_NUM_THREADS='1'; $env:MKL_NUM_THREADS='1'; $env:OMP_NUM_THREADS='1'
$Out = 'runs/NEW-nc-physical-review'  # 必ず未存在
New-Item -ItemType Directory $Out | Out-Null
# 次のPythonで分割をSHA照合して結合する。改行処理やGet-Contentでバイナリを扱わない。
$env:NC_REVIEW_OUT = $Out
& $Py -c 'import hashlib,json,os,pathlib; r=pathlib.Path("interim/2026-11-13-nc-lambda-physical-limit"); m=json.loads((r/"evidence-manifest.json").read_text()); chunks=[(r/p["name"]).read_bytes() for p in m["parts"]]; assert all(len(b)==p["bytes"] and hashlib.sha256(b).hexdigest()==p["sha256"] for b,p in zip(chunks,m["parts"])); b=b"".join(chunks); assert hashlib.sha256(b).hexdigest()==m["archive_sha256"]; (pathlib.Path(os.environ["NC_REVIEW_OUT"])/m["archive_name"]).write_bytes(b)'
if ($LASTEXITCODE -ne 0) { throw 'archive parts differ' }
$Archive = "$Out/diagnostic-evidence.zip"
Expand-Archive $Archive "$Out/evidence"
```

`evidence-manifest.json`のarchive SHAと全entry SHAを照合する。既存のmanifestを再生成して通さない。

```powershell
$env:NC_REVIEW_OUT = $Out
& $Py -c 'import hashlib,json,os,pathlib,zipfile; r=pathlib.Path("interim/2026-11-13-nc-lambda-physical-limit"); m=json.loads((r/"evidence-manifest.json").read_text()); p=pathlib.Path(os.environ["NC_REVIEW_OUT"])/m["archive_name"]; assert hashlib.sha256(p.read_bytes()).hexdigest()==m["archive_sha256"]; z=zipfile.ZipFile(p); assert set(z.namelist())==set(m["entries"]); assert all(hashlib.sha256(z.read(n)).hexdigest()==s for n,s in m["entries"].items()); print(len(m["entries"]))'
if ($LASTEXITCODE -ne 0) { throw 'archive differs' }
```

## 保存済み1/128の全GP・荷重と1/256の実受理2行

求解しない。想定は数秒〜数十秒。Linux保存checkpointをWindows本番ソルバーで再開する操作ではない。

```powershell
$E = "$Out/evidence"
& $Py tools/analyze_nc_fold_path.py --path "$E/native128" --output "$Out/reanalyzed128"
if ($LASTEXITCODE -ne 0) { throw 'fold audit' }
& $Py tools/nc_physical_limit_audit.py --path "$E/native128" --origin-checkpoint "$E/origin/prior/fresh-head-midpoint/job/recovery-checkpoint-0.json" --origin-mesh "$E/origin/prior/fresh-head-midpoint/input/mesh.json" --output "$Out/physical128"
if ($LASTEXITCODE -ne 0) { throw 'physical GP audit' }
& $Py tools/nc_physical_limit_audit.py --partial-history "$E/series/256" --accepted-rows "$E/series/256/history/accepted-tail.jsonl" --output "$Out/physical256-partial"
if ($LASTEXITCODE -ne 0) { throw 'partial GP audit' }
```

出力先の`load-audit.json`、`gp-audit.json`、`weak-mode.json`を読む。全162活動点、荷重の成分、せん断込みq、dry/wet、塑性体積ひずみの符号、楕円の降伏値を比較する。新規参照値／一致許容値を作らない。Linux/Windowsの丸め差は両実値を保存し、元の力1e-8・真の線形1e-10を緩めない。

1/256の2行は元の501行の499・500番目（0始まり）をそのまま切り出したもの。`trace-slice-receipt.json`に元全体SHAと各行SHAを保存した。残りの前史や新規fold場がZIPにあるとは扱わない。最大λ=0.951171875は時間切れ前の受理点であり、fold荷重・崩壊・完了を示さない。生成されたJSONのソースパスや部分配列の選択情報は実行場所によって変わるので、数値場と来歴を分けて比較する。

## 将来1/256のfold場を追加監査する場合だけ

PC-2の実際の`history/run.json`・`accepted.jsonl`・同じ入力がPCにある場合は、それを編集・再封印せず既存手順を使う。今回のクラウドでは未取得。元ファイルがなければWindowsで最初から作る。前史はPCで数分〜数十分以上かかり得る。クラウドの45秒を完了時間の予測にしない。長いM・側圧・自重系列へ自動的に進めない。

```powershell
cargo build --locked --release --example nc_history_trace
if ($LASTEXITCODE -ne 0) { throw 'trace build' }
& $Py tools/prepare_nc_fold_observer.py --output "$Out/observer"
if ($LASTEXITCODE -ne 0) { throw 'observer build' }
$Probe = (Resolve-Path "$Out/observer/target/release/examples/material_diagnostic.exe").Path
& $Py tools/prepare_nc_fold_series.py --output "$Out/series"
# 必要な分母だけ。物理・力/線形許容値・入力既定は変えない。
$Input256 = "$Out/series/256"
$History256 = "$Input256/history"
& target/release/examples/nc_history_trace.exe trace "$Input256/mesh.json" "$Input256/ground.json" "$Input256/stages.json" $History256
if ($LASTEXITCODE -ne 0) { throw 'trace process failed; retain logs and partial rows' }
Copy-Item "$Input256/mesh.json" "$History256/mesh.json"
& $Py tools/nc_fixed_inputs.py --history $History256 --output "$Out/fixed256.json"
if ($LASTEXITCODE -ne 0) { throw 'fixed inputs' }
& $Py tools/nc_fold_path.py --map native --fresh-native $History256 --binary $Probe --step-size 2e-5 --count 32 --budget-seconds 60 --output "$Out/native256"
if ($LASTEXITCODE -ne 0) { throw 'local path failed; not physical collapse' }
& $Py tools/analyze_nc_fold_path.py --path "$Out/native256" --output "$Out/native256-audit"
if ($LASTEXITCODE -ne 0) { throw 'local audit failed' }
& $Py tools/nc_physical_limit_audit.py --path "$Out/native256" --origin-checkpoint "$History256/checkpoint-stage-3.json" --origin-mesh "$Input256/mesh.json" --output "$Out/physical256-fold"
```

局所32点はPCで概ね数秒〜数十秒の目安だが、`local_fold_bracket`が出る保証はない。点数不足・時間切れなら未確定を保持する。`static_limit_confirmed`・`passed`はfalseのまま。基準外残差、入力不一致、欠損、棄却や時間切れを成功としない。新規上下界・極限定理の適用可否を未検証のままNC-01の受入へ進めない。PCからpushしない。
