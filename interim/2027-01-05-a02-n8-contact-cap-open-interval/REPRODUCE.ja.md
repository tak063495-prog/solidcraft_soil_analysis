# N=8接触上限・open最初の区間の任意再現

PR #13、`codex/plastic-spatial-convergence`、起点 `612f0f64c704f0ef0dee950a81a10f7e312a745c` 以降。マージしない。
2027-01-05はユーザー指定の保存先ラベル。Windows PCでは未実行、push後のWindows CIは未確定。

GPD WIN Max 2・使用可能メモリ約48GBの範囲で、同じ81節点・24 Tet10・648材料GP・32接触GPだけを扱う。
768/6144要素、メモリ40GB超、1経路1日超の新しい手順は作らない。
上限を超える一般化は「規模上限のため確認不可・未達」。許容値・参照値を変えない。
保守的予算はビルド10〜40分/8GB以内、1主経路16分以内、全28主経路と14観測なし対照で3時間/8GB以内。
Windows実測保証ではなく、小型入力と既存プロセス予算に基づく目安。時間切れは数値未完了として保存する。

## 新しい封印と必要なコードだけをfiltered取得

未存在の取得先を使い、過去の出力は上書きしない。初回全履歴転送を避ける。

```powershell
git clone --filter=blob:none --single-branch --branch codex/plastic-spatial-convergence --depth=1 --no-checkout https://github.com/tak063495-prog/solidcraft_soil_analysis.git NEW-a02n8
Set-Location NEW-a02n8
git config --local core.autocrlf false
git sparse-checkout set --no-cone '/Cargo.toml' '/Cargo.lock' '/rust-toolchain.toml' '/requirements*.txt' '/src/' '/examples/' '/tools/' '/interim/2026-11-25-joint-birth-projection-a02/*.py' '/interim/2026-12-21-a02-joint-birth-cause-diagnosis/*.py' '/interim/2026-12-29-a02-joint-birth-cloud-followup/*.py' '/interim/2026-12-29-a02-joint-birth-cloud-followup/a02follow_service.rs' '/interim/2027-01-05-a02-n8-contact-cap-open-interval/'
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
$D='interim/2027-01-05-a02-n8-contact-cap-open-interval'
New-Item -ItemType Directory -Force work | Out-Null
& $Py "$D/a02n8_evidence.py" --verify $D --extract work/NEW-a02n8-reopened
if ($LASTEXITCODE -ne 0) { throw '封印SHA不一致。先へ進まない' }
```

全part・結合ZIP・全member SHAを照合する。新封印の`old-controls/`は12-29の選択memberを同一バイトで保持し、
`input-binding.json`に元part/結合ZIP/member SHAを記録している。旧封印を書き換えず、旧値を新規求解に数えない。
本作業の観測再現に旧アーカイブ全体の転送は不要。

## 封印からの軽い観測再現

```powershell
& $Py "$D/a02n8_report.py" --run work/NEW-a02n8-reopened/formal --old work/NEW-a02n8-reopened/old-controls --output work/NEW-a02n8-reported
if ($LASTEXITCODE -ne 0) { throw '観測再現未完了' }
(Get-FileHash work/NEW-a02n8-reported/measurements.json -Algorithm SHA256).Hash
(Get-FileHash "$D/measurements.json" -Algorithm SHA256).Hash
```

封印の同じ入力からは2つのSHAが一致する。新しい原因診断や全体求解ではなく、証拠の再読込検査である。
`difference-vectors.json`、`branch-switches.json`、`accepted-traces.json`も同じ封印から再現する。
集計は工程→共通`i/2`→GP/節点→成分のC順、`math.fsum`。Windowsの集計順差を許容値で吸収しない。

## 指定2件だけの新規全体平衡

```powershell
& $Py "$D/a02n8_build.py" --output work/NEW-a02n8-build
if ($LASTEXITCODE -ne 0) { throw '隔離ビルド未完了' }
& $Py "$D/a02n8_run.py" batch --output work/NEW-a02n8-formal --build work/NEW-a02n8-build
if ($LASTEXITCODE -ne 0) { throw '数値未完了あり。process/failure/受理prefixを保持' }
& $Py "$D/a02n8_report.py" --run work/NEW-a02n8-formal --old work/NEW-a02n8-reopened/old-controls --output work/NEW-a02n8-new-report
if ($LASTEXITCODE -ne 0) { throw '新規求解の観測未完了' }
```

N=8/M=1で、接触上限0.0005 mと通常0.05 mを両実装のP×R全4通りで比較する16主経路。
N=4/8・通常上限・M=4で、openの最初の区間だけを2/4/8分割する両実装12主経路。
本体14主経路は観測なしwholeも別求解し、elapsed_secondsだけを除いた全Runを厳密比較する。
観測付きと観測なしを含む42経路を、主経路28本と混同しない。
N=4のP×R、M×N全交差、ULP/cone、12-21の原因診断は新規実行しない。

`global-input-endpoints.json`は新規求解前に作る入力で、保存jumpや保存受理列から作らない。
N=4は0〜0.25、N=8は0〜0.125の区間を等分し、それ以後のopen格子と他工程の格子を保持する。
Rustは専用ソースコピーで、この入力の次のprogressを選ぶ箇所だけを明示切替する。
Python独立側も同じ入力progressで全拘束・全体残差・Newtonを新規に解く。
最小増分の入力下限を子区間幅へ下げるが、力・線形残差・接触則・仕事定義・接線差分幅は変えない。
両実装とも不収束を受理しない。未完了を範囲拒否・物理破壊・合格へ変換しない。
区間上端は観測値であり、連続の厳密切替根ではない。

元の0.157080515は12-29のN4→8/M=4・共通i/2の保存値として保持する。
新しい差ベクトルと比較点集合を別記し、旧数値の再集計だけから新しい結論を作らない。
新規求解はOS間の最終ビット一致を保証せず、生差と実行環境を保存する。

## 全件検証と固定Git出典

全件試験は上の最小sparseとは別で、旧アーカイブ実ファイル約30GBとビルド・展開のディスク余裕が必要。
12-29〜01-04のgit-recoveryと本フォルダの`validation/git-recovery/`を参照する。
まずtextだけをfiltered取得し、必要なアーカイブを固定OID単位で取得・照合する。
確認済みの単一blob promisorキャッシュだけを削減し、実ファイルや旧証拠は消さない。
最終sparse除外は既存CIのpostfold74パスだけ。新しい封印を除外しない。
P6-01の固定出典コミットと、P6-03の固定BASE/tree・29固有blobもdepth=1/blob:noneで限定取得する。
生ファイルがあることとGit出典が揃っていることは別々に検査する。初回全履歴転送で代用しない。

```powershell
# 必要アーカイブと固定Git出典をすべて配置・SHA検証した全件用checkoutで実行する。
& $Py "$D/a02n8_validation.py" rust --output work/NEW-a02n8-full/rust
if ($LASTEXITCODE -ne 0) { throw 'Rust検証器未完了。process JSONを確認' }
& $Py "$D/a02n8_validation.py" python --output work/NEW-a02n8-full/python
if ($LASTEXITCODE -ne 0) { throw 'Python検証器未完了。process JSONを確認' }
& $Py "$D/a02n8_validation.py" collect --output work/NEW-a02n8-full
```

collectのpassed=falseを成功に読み替えない。成功/失敗/error/未実行/skipを固有IDと照合する。
PowerShellを使う既存試験にはPowerShellが必要。読取専用HOMEではwork内へXDG_CACHE_HOME/CONFIG_HOME/DATA_HOMEを設定する。
HOMEは変更しない。単独スモーク・封印再読込・登録後検査を全件件数へ二重加算しない。
重い再計算は封印に保存し、rust-integration-1への重い追加試験は0件。
文書比較はuniversal newlineで、実CRLF版も試験する。封印バイトには改行変換をしない。

S1/S2の初期応力・参照配置/所有/世代/独立平衡、S3/S4の不適格初期力の物理過程、
S5の先端力の受持ち・一意な投影は仕様・実装の判断待ち。現行拒否を保持し、初期力クリップを追加しない。
A-01〜A-04、B-01〜B-03、C-01〜C-04、NC-01、P6は全体未達のまま。
push後のWindows CIは未確定。

A-02は未達、passed=false。
