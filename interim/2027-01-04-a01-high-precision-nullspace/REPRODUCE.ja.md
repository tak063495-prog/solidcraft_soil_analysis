# 48DOF高精度診断の任意再現

PR #13、codex/plastic-spatial-convergence、起点5984f0955c37542f808335cb17cae5f766356945以降。マージしない。
保存点・旧世代・変位を使用する観測専用手順。新しい全体求解は0件。
GPD WIN Max 2、使用可能メモリ約48GBの範囲で、18 Tet4／48縮約DOFだけを扱う。
Windows PCでは未実行。保守的予算は高精度5対照合わせて30分以内・メモリ2GB以内。
48DOFの行列・入力規模に基づく予算でありWindows実測保証ではない。実測秒数は封印内各raw JSONに記録する。
768/6144要素、40GB超、1経路1日超の新しい手順は作らない。
上限を超える確認は「規模上限のため確認不可・未達」とし、元の許容値を維持する。

## 必要な新フォルダだけをfiltered取得

未存在の取得先を使う。初回全履歴転送を避け、過去の作業や証拠を上書きしない。

```powershell
git clone --filter=blob:none --single-branch --branch codex/plastic-spatial-convergence --depth=1 --no-checkout https://github.com/tak063495-prog/solidcraft_soil_analysis.git NEW-a01hp
Set-Location NEW-a01hp
git config --local core.autocrlf false
git sparse-checkout set --no-cone '/interim/2027-01-04-a01-high-precision-nullspace/' '/interim/2026-12-28-a01-joint-combo-cloud-followup/a01follow_contact.py'
git switch codex/plastic-spatial-convergence
$D='interim/2027-01-04-a01-high-precision-nullspace'
py -3.13 -m venv .venv
$Py='.venv/Scripts/python.exe'
New-Item -ItemType Directory -Force work | Out-Null
& $Py "$D/a01hp_evidence.py" --verify $D --extract work/NEW-a01hp-reopened
if ($LASTEXITCODE -ne 0) { throw '封印SHA不一致。先へ進まない' }
```

新アーカイブの全part/member SHAと、封印された入力の旧part/member出典はinput-binding.jsonで照合する。
本診断に必要な停止点入力は新封印のsource/以下にも同一バイトで保存してあり、旧100MB封印を変更していない。
任意精度再計算には次の固定追加依存だけを使う。CIでこれを必須にしない。

```powershell
& $Py -m pip install mpmath==1.3.0
if ($LASTEXITCODE -ne 0) { throw '任意精度依存の準備未完了' }
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
$S='work/NEW-a01hp-reopened/source'
New-Item -ItemType Directory -Force work/NEW-a01hp-results | Out-Null
& $Py "$D/a01hp_run.py" --source $S --output work/NEW-a01hp-results/rebuilt50.json --dps 50
if ($LASTEXITCODE -ne 0) { throw '50桁未完了。ログと部分出力を保持' }
& $Py "$D/a01hp_run.py" --source $S --output work/NEW-a01hp-results/rebuilt80.json --dps 80
if ($LASTEXITCODE -ne 0) { throw '80桁未完了。ログと部分出力を保持' }
& $Py "$D/a01hp_run.py" --source $S --output work/NEW-a01hp-results/binary-rebuilt50.json --dps 50 --binary64-inputs
if ($LASTEXITCODE -ne 0) { throw 'binary64入力50桁対照未完了' }
& $Py "$D/a01hp_run.py" --source $S --output work/NEW-a01hp-results/binary-rebuilt80.json --dps 80 --binary64-inputs
if ($LASTEXITCODE -ne 0) { throw 'binary64入力対照未完了' }
& $Py "$D/a01hp_run.py" --source $S --output work/NEW-a01hp-results/frozen80.json --dps 80 --frozen-local-inputs
if ($LASTEXITCODE -ne 0) { throw '局所固定入力対照未完了' }
```

各実行が解析とFD h/h2/h4、保存予測子前と解法行列の直接SVDを記録する。
全スペクトル、左右2ベクトル、残差・誕生保持力・外力の生の成分はraw JSON。
基底の符号固定は右ベクトル最大成分を正とする。縮退2方向は個別の向きではなく射影行列で比較する。
約1e-47から約1e-76へ下がる値は「0の数値残り」。物理的な微小支持や受入許容値にしない。
生成時の実測秒数・数学的に任意な縮退基底・最下位桁までの別環境バイト一致を要求していない。
固定順演算と保存結果のSHAは区別し、集計順差を許容値で吸収しない。

## 封印からの軽い観測再現

高精度計算なしで同じ観測JSONを再生成できる。これは新しい原因同定や新しい求解ではなく、証拠の再読込検査である。

```powershell
& $Py "$D/a01hp_report.py" --raw work/NEW-a01hp-reopened/raw --output work/NEW-a01hp-measurements.json
if ($LASTEXITCODE -ne 0) { throw '観測再現未完了' }
(Get-FileHash work/NEW-a01hp-measurements.json -Algorithm SHA256).Hash
(Get-FileHash "$D/measurements.json" -Algorithm SHA256).Hash
```

上の2 SHAは一致する。原資料のHTTP取得記録は封印内preparation/sources/。
Wriggers等の書籍本文は未取得・書誌だけ。本文が取得できた資料と混同しない。
倍精度再組立を任意再実行する場合は、リポジトリの固定requirements-reference.lock.txtに従ってnumpyも準備し、
a01hp_double.pyを同じsourceへ向ける。数値SVDは計算方法・ライブラリにより丸め域の最小値が変わり得る。
元行列の高精度SVDと高精度再組立を取り違えない。

## 全件検証・台帳用の取得

全件試験は本診断だけのsparse取得とは別で、既存アーカイブ全体を必要とする。
A-02〜C-04のgit-recoveryと今回のvalidation/git-recoveryを参照する。
固定P6コミットとP6-03 BASE/tree/29固有blobはdepth=1/blob:noneで追加取得する。
アーカイブは固定OID単位で取得し、実ファイルのGit blob/SHA-256一致を確認してから、
再取得できる単一blob promisorキャッシュだけを削減する。HOMEや既存科学ファイルは変更しない。
最終sparse除外は既存CIのpostfold74パスだけ。必要な新封印を追加除外しない。
このため一時ディスクには旧アーカイブ実ファイル約30GBとビルド・展開の余裕が必要。
一括取得途中の重複キャッシュを残してディスクを枯渇させない。

```powershell
$env:CARGO_BUILD_JOBS='2'
$env:RUST_TEST_THREADS='2'
& $Py "$D/a01hp_validation.py" rust --output work/NEW-a01hp-validation/rust
if ($LASTEXITCODE -ne 0) { throw 'Rust検証器未完了。各process JSONも確認' }
& $Py "$D/a01hp_validation.py" python --output work/NEW-a01hp-validation/python
if ($LASTEXITCODE -ne 0) { throw 'Python検証器未完了。各process JSONも確認' }
& $Py "$D/a01hp_validation.py" collect --output work/NEW-a01hp-validation
```

各processのexit_codeとsummaryを確認する。collectのpassed=falseを成功に読み替えない。
実施済み全件の確定件数はvalidation.json。単独スモーク・封印再読込・登録後検査は重複加算しない。
PowerShellを使う既存試験にはPowerShellが必要。読取専用HOMEではwork内にXDG_CACHE_HOME/CONFIG_HOME/DATA_HOMEだけを設定する。
mpmathを入れないCIでも軽い新規試験は動作し、高精度再計算をrust-integration-1へ追加していない。
文書比較はuniversal newline、封印バイトは改行変換しない。
push後のWindows CIは未確定。

停止した枝の荷重不整合とFD微小剛性・丸めの分離だけを限定確認した。
(d)物理的不安定は判定不能。PCへの必須移管はなく、現行方針で追加実行すべき既知のクラウド診断は残していない。
湿潤先端・複数先端/部材複合の仕様と数値方針はユーザー判断待ちである。

A-01は未達、passed=false。
