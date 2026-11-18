# 明示freshせん断の任意再現

PR #13、codex/plastic-spatial-convergence、起点cdaf850c0347c8d841decce7ddff608d03ae7933。マージしない。
2027-01-07は指定保存先ラベル。Windows PC未実行、push後のWindows CIは未確定。

GPD WIN Max 2（使用可能メモリ約48GB）で同じ54節点・12 Tet10だけを扱う。
保守的目安はビルド8GB以内、求解/観測12GB以内、全行列3時間程度。Windows実測保証ではない。
1子プロセス900秒で停止しprefix/ログを残す。不収束は数値未完了として扱う。
768/6144要素・40GB超・1経路1日超の計算手順は追加しない。
一般空間確認のうち上限超が必要な点は、規模上限のため確認不可・未達。

## 新規保存先への取得と封印確認

```powershell
git clone --filter=blob:none --single-branch --branch codex/plastic-spatial-convergence --depth=1 --no-checkout https://github.com/tak063495-prog/solidcraft_soil_analysis.git NEW-a03shear
Set-Location NEW-a03shear
git config --local core.autocrlf false
git sparse-checkout set --no-cone '/Cargo.toml' '/Cargo.lock' '/rust-toolchain.toml' '/requirements*.txt' '/src/' '/examples/' '/tools/' '/interim/2026-11-26-wet-soil-birth-darcy-a03/*.py' '/interim/2026-12-30-a03-wet-soil-birth-cloud-followup/*.py' '/interim/2027-01-07-a03-explicit-shear-initial-stress/'
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
$D='interim/2027-01-07-a03-explicit-shear-initial-stress'
New-Item -ItemType Directory -Force work | Out-Null
& $Py "$D/a03shear_evidence.py" --verify $D --extract work/NEW-a03shear-reopened
if ($LASTEXITCODE -ne 0) { throw '封印SHA不一致' }
```

全part・結合ZIP・全memberのサイズ/SHAを照合する。旧全アーカイブの転送はこの観測再現には不要。
formalに新しい入力・全受理場・残差・Newton回数・本体全Run・独立結果・ソース/実行版SHAを含む。
retainedに一次資料取得失敗と今回の求解前監査・報告の生差ベクトルを保持する。

## 保存場からの同一集計

```powershell
& $Py "$D/a03shear_report.py" --run work/NEW-a03shear-reopened/formal --output work/NEW-a03shear-report
if ($LASTEXITCODE -ne 0) { throw '観測再現未完了' }
(Get-FileHash work/NEW-a03shear-report/measurements.json -Algorithm SHA256).Hash
(Get-FileHash "$D/measurements.json" -Algorithm SHA256).Hash
```

工程→共通i/4の16座標→物理接続順→成分、C-order、math.fsum、left−rightを固定する。
差が出ればOS/バージョン/差分を保存し、一致確認は未達とする。Windowsの集計順差を許容値で吸収しない。
再組立/SVDを含む最終ビットのOS間一致を保証しない。保存場の集計は新しい全体求解に数えない。

## 固定行列を入力から新規求解

```powershell
cargo build --locked --release --bin soil-fem-3d
if ($LASTEXITCODE -ne 0) { throw '通常CLIビルド未完了。古い実行版を使用しない' }
& $Py "$D/a03shear_build.py" --output work/NEW-a03shear-build
if ($LASTEXITCODE -ne 0) { throw '隔離ビルド未完了' }
& $Py "$D/a03shear_run.py" --output work/NEW-a03shear-formal --build work/NEW-a03shear-build
if ($LASTEXITCODE -ne 0) { throw '数値未完了又は事前拒否あり。failure/process/prefixを保持' }
& $Py "$D/a03shear_report.py" --run work/NEW-a03shear-formal --output work/NEW-a03shear-new-report
if ($LASTEXITCODE -ne 0) { throw '観測未完了' }
```

h0.5、せん断0/τxz/τyz各±1%/±5%の9入力×N4/8/16を両実装で54経路2016状態。
独立N8/GL7が9経路288状態、正式計63経路2304状態。
本体の観測なしwhole27経路1008状態と通常CLI9経路288状態は別集計。
生産Rustを読み取りコピーするだけで専用probeは明示実行、通常CLI既定オフ。
Lower・非活動Upperは明示対角定数、fresh Upperだけ明示せん断を変える。
全入力の初期降伏値・種応力の等価節点力を最初の求解前に保存する。定数応力の内部div=0と境界内力を区別する。
配置後の不釣合いと保持力は本体prepared/transition、独立preparedに保存する。
新しい再配分・補正はせず、現行の保持力解放/全体平衡を両実装が新規に解く。
独立側は直接線形求解1回/状態、Newton0回。本体は実attemptのNewton回数を保存する。
各受理点でMC弾性と閉鎖固着を検査する。途中失敗は数値未完了、許可域外は求解せず拒否する。
12-22原因診断、12-30(1)〜(3)、01-06不等厚×Kx≠Ky行列を再実行する手順ではない。

## 全件検証とP6登録

全件用には別checkoutと旧実証拠約30GB＋ビルド/展開余裕が必要。
初回全履歴転送を避けfiltered取得し、固定先端のパス/Git blob/SHAを照合する。
今回のwork側復旧スクリプトと検証ログはvalidationの封印に保持する。
アーカイブは一時のbare transportでblobだけを取得して照合し、作業用コピー後に一時transportを片付けると
Git objectと作業ファイルの二重保存を避けられる。読取専用HOMEを変更せずキャッシュはwork側へ指定する。
既存CIで除外するNC postfold74part以外の欠落をskipにしない。
新封印はpython-core-0/1の実sparse規則で選択される。パスはPOSIX形式/NUL区切りで検査する。

```powershell
# 旧証拠とP6固定Git出典を配置/SHA確認した全件用checkoutで実行。
& $Py "$D/a03shear_validation.py" rust --output work/NEW-a03shear-full/rust
if ($LASTEXITCODE -ne 0) { throw 'Rust検証器未完了' }
& $Py "$D/a03shear_validation.py" python --output work/NEW-a03shear-full/python
if ($LASTEXITCODE -ne 0) { throw 'Python検証器未完了' }
& $Py "$D/a03shear_validation.py" collect --output work/NEW-a03shear-full
```

Python全7群の発見/割当/実行IDを照合し、成功/失敗/error/未実行/skipを分ける。
PowerShellを使う既存試験にはPowerShellが必要。単独確認や登録後監査を全件固有数へ二重加算しない。
P6-01へ科学記録と確定validationを別の固定コミット・Git blob・SHA-256・型付きboolで追記する。
文書末尾比較はuniversal newlineを使い、封印バイトは改行変換しない。
重い再計算を新試験/rust-integration-1へ追加しない。ソフトウェア成功をA-03受入へ転用しない。

S1 Darcy後BIRTHのσ'/p・活動面/所有/誕生水量、S2部分水没・水圧割当て、S3 MCC形成史/OCR/pc/e、
S4旧履歴・配置・保持力、S5水量/収支nullはユーザーの仕様判断待ち。拒否とnullを維持する。
A-01〜A-04、B-01〜B-03、C-01〜C-04、NC-01、P6は全体未達。push後のWindows CIは未確定。

A-03は未達、passed=false。
