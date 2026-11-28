# Windows再現・長い系列（PCからpushしない）

未取得のfold・材料点限界・連続体上下界を成功として補完しない。PCでもNC-01入力、src、参照値、許容値、旧証拠は編集しない。自重オフは現行入力で拒否されることが今回の結果であり、0を微小値へ置換しない。

## 固定checkoutと準備

PowerShell。`$Target`は最終報告の40桁pushコミットSHA。空のworktreeと未存在の出力先を用いる。

```powershell
$Repo = 'D:/Work/solidcraft_soil_analysis'
$Target = '<最終報告の40桁SHA>'
$Work = 'D:/Work/scsa-nc-collapse-review'
git -C $Repo fetch origin codex/plastic-spatial-convergence
git -C $Repo worktree add --detach $Work $Target
Set-Location $Work
git merge-base --is-ancestor 2a826c48a73049f53d39d520ca8fc64bb93f440f HEAD
if ($LASTEXITCODE -ne 0) { throw 'wrong base' }
py -3.13 -m venv .venv
$Py = (Resolve-Path '.venv/Scripts/python.exe').Path
& $Py -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw 'existing locked dependencies' }
$env:OPENBLAS_NUM_THREADS='1'; $env:MKL_NUM_THREADS='1'; $env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'; $env:PYTHONIOENCODING='utf-8'
$Out = 'runs/NEW-nc-collapse-pc'   # 必ず未存在
New-Item -ItemType Directory $Out | Out-Null
cargo build --locked --release --example nc_history_trace
if ($LASTEXITCODE -ne 0) { throw 'trace build' }
& $Py tools/prepare_nc_fold_observer.py --output "$Out/observer"
if ($LASTEXITCODE -ne 0) { throw 'isolated observer build' }
$Trace = (Resolve-Path 'target/release/examples/nc_history_trace.exe').Path
$Observer = (Resolve-Path "$Out/observer/target/release/examples/material_diagnostic.exe").Path
```

新しいSOCP依存は不要。LPは既存SciPy/HiGHS。連続体FELA用の別応力空間・別最適化依存を導入するなら、元問題・境界traction・全評価点保証・固定pcと履歴MCCの関係・依存固定と再現手順を別途レビューする。この作業には導入していない。

## 保存証拠のSHA照合と展開

```powershell
$env:NC_COLLAPSE_OUT=$Out
& $Py -c 'import hashlib,json,os,pathlib,zipfile,io; r=pathlib.Path("interim/2026-11-13-nc-lambda-collapse-classification"); m=json.loads((r/"evidence-manifest.json").read_text()); parts=[(r/p["name"]).read_bytes() for p in m["parts"]]; assert all(len(b)==p["bytes"] and hashlib.sha256(b).hexdigest()==p["sha256"] for b,p in zip(parts,m["parts"])); b=b"".join(parts); assert hashlib.sha256(b).hexdigest()==m["archive_sha256"]; z=zipfile.ZipFile(io.BytesIO(b)); assert set(z.namelist())==set(m["entries"]); assert all(hashlib.sha256(z.read(n)).hexdigest()==s for n,s in m["entries"].items()); z.extractall(pathlib.Path(os.environ["NC_COLLAPSE_OUT"])/"evidence"); print(len(m["entries"]))'
if ($LASTEXITCODE -ne 0) { throw 'archive hash differs' }
$Saved = "$Out/evidence/previous/native128"
& $Py tools/nc_collapse_classify.py --path $Saved --output "$Out/recomputed-classification"
if ($LASTEXITCODE -ne 0) { throw 'GP/bounds/stability diagnostic failed' }
& $Py tools/nc_collapse_material.py --path $Saved --binary $Observer --output "$Out/recomputed-material"
if ($LASTEXITCODE -ne 0) { throw 'material diagnostic failed' }
```

GP、bounds、stabilityの数値場を保存JSONと比較する。Linux/Windowsの浮動小数差と来歴パス・実行時間を区別して双方を保存し、元の求解基準を緩めない。`bounds.json`の上下側値はGP弱釣り合い・固定pcの離散問題。連続体FELA・硬化MCCの物理的上下界へ昇格させない。材料点の局所逆解法停止も構成則の限界へ読み替えない。

## 長い36条件前史・局所経路

以下はクラウドの20秒では未取得だった前史等をPCで観測する手順。1系列15分＋局所5分が上限で、完了時間の保証ではない。時間切れを正常な解・CI成功・崩壊にしない。自重オフ18条件は同じ入力拒否として保存される。

```powershell
& $Py tools/nc_collapse_sweep.py --output "$Out/series" --trace-binary $Trace --observer $Observer --budget-seconds 900 --fold-budget-seconds 300
if ($LASTEXITCODE -ne 0) { throw 'sweep driver error; retain all partial files and logs' }
```

`result.json`の`timed_out`・実終了コード・`full_removal_reached`・`lambda_fold`を読む。ドライバ終了0は全系列の完了を意味しない。全撤去到達もNC-01の独立全場・刻み収束・完了再開の受入ではない。

局所128点の点数上限が残った場合だけ、同じ保存前史から既存弧長法の長い局所観測を別フォルダへ作る。元checkpointの編集・再封印・既存foldの上書きはしない。

```powershell
$env:NC_COLLAPSE_SERIES="$Out/series"
$env:NC_COLLAPSE_OBSERVER=$Observer
& $Py -c 'import sys,os,pathlib; sys.path.insert(0,"tools"); from nc_collapse_sweep import follow; from nc_collapse_classify import read,write; base=pathlib.Path(os.environ["NC_COLLAPSE_SERIES"]); binary=pathlib.Path(os.environ["NC_COLLAPSE_OBSERVER"]); folders=[p for p in base.iterdir() if p.is_dir() and (p/"result.json").exists()]; [(write(p/"pc-long-result.json",follow(p.resolve(),binary.resolve(),1800.,4096,"fold-PC-long"))) for p in folders if read(p/"result.json").get("last_stage")==3 and read(p/"result.json").get("last_progress",0)>.85 and not read(p/"result.json")["full_removal_reached"] and read(p/"result.json")["lambda_fold"] is None]'
if ($LASTEXITCODE -ne 0) { throw 'local path failed; inconclusive, retain evidence' }
```

最大点数4096、最大30分はPC専用。λ増分／経路接線の正→負、元の釣り合い・真の線形残差、border正則性を実際の保存行列で検査する。単に最大サンプルをλ_foldにしない。新系列の正確な点数・ds・前史入力SHAを比較し、基点1/128の0.9514000964へ勝手に統一しない。

## 96 Tet10の前史（PCで別入力）

```powershell
$env:NC_COLLAPSE_FINE="$Out/refined96"
& $Py -c 'import os,sys,pathlib; sys.path.insert(0,"tools"); from nc_collapse_sweep import profile; from refine_tet10 import refine; from nc_collapse_classify import write; p=pathlib.Path(os.environ["NC_COLLAPSE_FINE"]); p.mkdir(); mesh,g,c=profile(1.4,100,True); fine,_=refine(mesh); c["settings"]["maximum_free_dofs"]=3*len(fine["nodes"]); [write(p/(n+".json"),v) for n,v in zip(("mesh","ground","stages"),(fine,g,c))]'
if ($LASTEXITCODE -ne 0) { throw 'refined input generation' }
& $Trace trace "$Out/refined96/mesh.json" "$Out/refined96/ground.json" "$Out/refined96/stages.json" "$Out/refined96/history" *> "$Out/refined96/trace.log"
if ($LASTEXITCODE -ne 0) { throw 'refined trace process error; not collapse' }
```

現在の掃引の上下界集計は元Tet10との同一形状を厳密に要求する。96 Tet10へ元54自由度の弱モードやRを流用して通さない。細メッシュの独自R/Bと実場を用いたfold監査は追加作業で、クラウドで比較完了したとはしない。長時間計算を止める場合は実プロセスの終了と受理ファイルを確認し、部分結果を保持する。

全てPCのローカル記録として保存し、PCからpushしない。新しい合格基準・参照値・完了分類は作らない。
