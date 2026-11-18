# 不等厚 × 明示異方的初期応力の任意再現

PR #13、codex/plastic-spatial-convergence、起点484f007490d4a042c7ad330de75881025c76a8df。マージしない。
2027-01-06は指定保存先ラベル、実施日2026-11-18 JST。Windows PCは未実行、push後のWindows CIは未確定。

GPD WIN Max 2（使用可能メモリ約48GB）用に、同じ54節点・12 Tet10だけを扱う。
保守的な目安はビルド8GB以内、求解/観測12GB以内、ビルド10〜40分、全行列3時間程度。
Windowsの実測保証ではない。各子プロセス900秒上限で停止し、数値未完了ならprefixとログを保持する。
768/6144要素・40GB超・1経路1日超の手順は追加しない。一般空間確認は規模上限のため確認不可・未達。

## 新規保存先へfiltered取得

```powershell
git clone --filter=blob:none --single-branch --branch codex/plastic-spatial-convergence --depth=1 --no-checkout https://github.com/tak063495-prog/solidcraft_soil_analysis.git NEW-a03aniso
Set-Location NEW-a03aniso
git config --local core.autocrlf false
git sparse-checkout set --no-cone '/Cargo.toml' '/Cargo.lock' '/rust-toolchain.toml' '/requirements*.txt' '/src/' '/examples/' '/tools/' '/interim/2026-11-26-wet-soil-birth-darcy-a03/*.py' '/interim/2026-12-30-a03-wet-soil-birth-cloud-followup/*.py' '/interim/2027-01-06-a03-unequal-thickness-anisotropic-k/'
git switch codex/plastic-spatial-convergence
py -3.13 -m venv .venv
$Py='.venv/Scripts/python.exe'
& $Py -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw '依存準備未完了' }
$env:CARGO_BUILD_JOBS='2'
$env:RUST_TEST_THREADS='2'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
$D='interim/2027-01-06-a03-unequal-thickness-anisotropic-k'
New-Item -ItemType Directory -Force work | Out-Null
& $Py "$D/a03aniso_evidence.py" --verify $D --extract work/NEW-a03aniso-reopened
if ($LASTEXITCODE -ne 0) { throw '封印SHA不一致。先へ進まない' }
```

全part、結合ZIP、全member SHAを照合する。旧全アーカイブの転送は今回の観測再現には不要。
封印のformalには新規の入力・全受理場・残差・Newton回数・本体全Run・独立結果・実行版/ソースSHAを含む。
retainedには準備失敗、旧CLIによる補助対照の失敗とその後の新規再求解を区別して残す。
旧保存値を新規求解数に数えない。

## 封印からの集計再現

```powershell
& $Py "$D/a03aniso_report.py" --run work/NEW-a03aniso-reopened/formal --output work/NEW-a03aniso-report
if ($LASTEXITCODE -ne 0) { throw '観測再現未完了' }
(Get-FileHash work/NEW-a03aniso-report/measurements.json -Algorithm SHA256).Hash
(Get-FileHash "$D/measurements.json" -Algorithm SHA256).Hash
```

同じLinux環境では2つのSHAとdifference-vectors.jsonのバイト一致を確認した。
WindowsではこのSHA比較自体も未確認であり、不一致なら差分と実行環境を保存して一致確認は未達とする。
再組立/拘束空間SVD等を含む最終ビットのOS間一致を保証せず、許容値で差を吸収しない。
これは観測の再読込であり、新しい全体求解ではない。
共通i/4の16座標、工程→物理接続順→成分、C-order、math.fsum、left−rightを固定する。
Windowsの集計順差を許容値で吸収しない。

## 同じ固定行列を新規求解

```powershell
cargo build --locked --release --bin soil-fem-3d
if ($LASTEXITCODE -ne 0) { throw '通常CLIビルド未完了。旧実行版を使用しない' }
& $Py "$D/a03aniso_build.py" --output work/NEW-a03aniso-build
if ($LASTEXITCODE -ne 0) { throw '隔離ビルド未完了' }
& $Py "$D/a03aniso_run.py" --output work/NEW-a03aniso-formal --build work/NEW-a03aniso-build
if ($LASTEXITCODE -ne 0) { throw '数値未完了あり。failure/process/prefixを保持' }
& $Py "$D/a03aniso_report.py" --run work/NEW-a03aniso-formal --output work/NEW-a03aniso-new-report
if ($LASTEXITCODE -ne 0) { throw '新規求解の観測未完了' }
```

h=0.4/0.5/0.6 × K55/K58/K85 × N4/8/16を両実装で54経路2016状態、独立N8/GL7を9経路288状態。
正式63経路2304状態。本体の観測なしwhole27経路1008状態、通常CLI9経路288状態を正式数と分ける。
生産Rustは読み取り/コピーだけ。専用観測probeを明示実行するため既定オフ。
初期応力はボディ別定数・せん断0で明示し、K0自動生成・深度連続分布・新しい支持/接続を追加しない。
独立側は入力のみから全体直接線形求解（Newton0回、直接求解1回/状態）。本体Newton回数は実attemptから保存。
全点で弾性MCと閉鎖固着を検査し、参照の範囲外や不収束を合格へ読み替えない。
12-30の2形状×2応力、E/F/Q、N×Mおよび12-22の原因診断は再実行しない。

## 全件検証とP6固定Git出典

上の最小sparseとは別に全件用checkoutが必要。旧証拠実ファイル約30GBとビルド/展開のディスク余裕を確保する。
A-02 01-05・A-01 01-04のgit-recovery、および今回validation/git-recoveryを参照する。
初回全履歴転送は避け、固定先端のアーカイブはパス/Git blob単位でfiltered取得し、元blob・SHA-256を照合する。
最終sparse除外は既存CIのNC postfold74partのみ。新封印はpython-core-0/1の実sparse規則で選択される。
P6-01の固定出典コミット、P6-03の固定BASE/treeと29固有blobも限定取得する。
実ファイルと固定Git出典は別々に検査し、欠落をskipにしない。HOMEは変更せず、必要ならXDG各キャッシュをwork内へ指定する。

```powershell
# 全件用checkoutで必要な旧アーカイブ・固定Git出典を配置/SHA照合した後に実行する。
& $Py "$D/a03aniso_validation.py" rust --output work/NEW-a03aniso-full/rust
if ($LASTEXITCODE -ne 0) { throw 'Rust検証器未完了' }
& $Py "$D/a03aniso_validation.py" python --output work/NEW-a03aniso-full/python
if ($LASTEXITCODE -ne 0) { throw 'Python検証器未完了' }
& $Py "$D/a03aniso_validation.py" collect --output work/NEW-a03aniso-full
```

PowerShellを使う既存試験にはPowerShellが必要。collectのpassed=falseを成功へ読み替えない。
成功/失敗/error/未実行/skipの固有IDを照合する。単独試験・登録後の確認は全件数へ二重加算しない。
新試験は入力・固定順・CRLF・SHAだけで、重いFEM再計算をrust-integration-1へ追加しない。
文書はuniversal newlineで比較し、封印バイトには改行変換をしない。

S1 Darcy後のσ'/p・所有/活動面・再平衡、S2自由水面/水圧割当て、S3 MCC形成史/OCR/pc/e、
S4履歴・配置・保持力、S5誕生水量・貯留/流量の契約はユーザーの仕様判断待ち。
PCへ未定義仕様を渡して解決済みにしない。現行拒否とnullを保持する。
A-01〜A-04、B-01〜B-03、C-01〜C-04、NC-01、P6は全体未達。push後のWindows CIは未確定。

A-03は未達、passed=false。
