# A-04：再現手順（Windows PC／Linux）

A-04全体の`passed=false`、双方向/湿潤Richards JOINTの拒否を保持する。
この手順は小型の既存受理経路を新規求解する。保存ZIPの監査と求解を区別する。
既存の失敗・参照値・許容値・NC・旧snapshot/hash・試験のskip設定を変更しない。

## Windowsの準備と全件検証

Rustはリポジトリの`rust-toolchain.toml`（1.99.1）、PythonはCIと同じ3.13を使う。
短いASCIIパスにcheckoutし、旧封印証拠も元Git blobのまま配置する。
全件試験には旧NC等のアーカイブ約8.26 GBと展開領域が必要で、空き50 GB以上を目安にする。
今回の疎checkoutでは必要466 blobを元SHA/SHA256と照合して配置した。
省略したアーカイブを空ファイル・再生成参照値で置き換えない。

```powershell
git fetch origin codex/plastic-spatial-convergence
git switch codex/plastic-spatial-convergence
git status --short
# READMEに示す検証コミットを含む先端であることを確認。force/reset/mergeは不要。
py -3.13 -m venv .venv
$Py = (Resolve-Path .venv/Scripts/python.exe).Path
& $Py -m pip install -r requirements-reference.lock.txt
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:RUST_TEST_THREADS='2'
$env:CARGO_BUILD_JOBS='2'
$env:PYTHONUTF8='1'
$Study = 'interim/2026-11-29-richards-wet-joint-a04'
cargo build --locked --release
if ($LASTEXITCODE -ne 0) { throw '本体ビルド失敗' }
cargo build --locked --release --examples
if ($LASTEXITCODE -ne 0) { throw '観測器/再hash補助のビルド失敗' }
& $Py "$Study/validate.py" --output runs/NEW-a04-PC-tests
if ($LASTEXITCODE -ne 0) { throw '全件成功ではない。全ログ/実件数を保持して原因を確認' }
```

`validate.py`はfmt、`clippy --locked --all-targets -- -D warnings`、
`cargo test --locked --release`、`cargo test --locked --release --examples`、Pythonの4群を実行する。
終了コードだけでなく、各`*-process.json`とPython inventoryの発見/割当/実行IDを照合する。
失敗した群のログを残し、必要な元blobの配置後は同じ群全体を新規フォルダーで再実行する。
同一IDの重複成功を加算しない。初回ビルド/全Rustは10〜60分、Python通常群は10〜40分、
旧NC2群はさらに20〜90分をPCの計画枠とする（Windows未測定、CPU/ディスクによる）。
時間切れ・中断・キャンセルを成功にしない。

旧NC mean/spatialはLinuxで4.440892e-16差により各setUpClassエラー・計8件未実行の履歴がある。
A-03では同じ元の厳密比較で8件成功した。今回の実件数は`validation.json`を参照し、
どちらの結果でも過去を改変しない。skip、許容差化、参照書換え、Windows分岐のLinux移植で通さない。
旧証拠の再監査成功はNC-01の数値受入を意味しない。
WindowsのNumPy longdoubleはLinuxの拡張精度と異なる場合がある。独立水理の厳しい残差を
満たせない場合はその生値・失敗を残し、基準や参照値を変更しない。今回のLinux成功をWindowsへ転用しない。

## 本体と独立参照を双方で新規求解

既存の受理範囲だけを使う。新しい出力ディレクトリを指定する。

```powershell
& $Py "$Study/study.py" --output runs/NEW-a04-PC-soil --steps 4 8 16
if ($LASTEXITCODE -ne 0) { throw '土の新規求解/比較/再開/改変のいずれかが未完了' }
& $Py "$Study/audit.py" --run runs/NEW-a04-PC-soil
if ($LASTEXITCODE -ne 0) { throw '全受理点・質量/荷重/凍結の保存場監査失敗' }
& $Py "$Study/convergence.py" --run runs/NEW-a04-PC-soil --steps 4 8 16
if ($LASTEXITCODE -ne 0) { throw '同一16物理時刻の比較に失敗' }
& $Py "$Study/joint_controls.py" --kind dry --steps 4 --output runs/NEW-a04-PC-dry
if ($LASTEXITCODE -ne 0) { throw '既存の乾燥JOINT対照が未完了' }
& $Py "$Study/joint_controls.py" --kind hydrostatic --steps 4 --output runs/NEW-a04-PC-hydrostatic
if ($LASTEXITCODE -ne 0) { throw '既存の静水圧JOINT対照が未完了' }
& $Py "$Study/refusals.py" --solver target/release/soil-fem-3d.exe --output runs/NEW-a04-PC-refusals
if ($LASTEXITCODE -ne 0) { throw '期待した拒否が確認できない。新しい受理にしない' }
```

土6対は272受理状態対。乾燥28対、静水圧16対で合計316対。
本体/独立の新規求解8対であり、別再開や保存場再監査をこれへ足さない。
土は5〜20分、乾燥/静水圧の観測・認証・復旧は各5〜20分程度をPC計画枠とする。
現行controllerのcheckpoint認証は受理済み経路も再演算するため、単なるJSON読込より重い。
実時間は個々のprocess.jsonとログを参照する。

|出力|内容|
|---|---|
|土`measurements.json`|全6経路、49場の最大差/相対L2・各受理点、time sensitivity、復旧/改変|
|`source-fingerprints.json` / `source-snapshot.zip`|新規求解開始時に存在したソース/SPECの指紋と内容|
|`bin/`, `binaries.json`|同一実行版の固定コピーとSHA。途中で入れ替えない|
|各`VG-nN` / `saturated-nN`の`native.json`|実CLIの全Run、初期平衡、全後続受理状態、台帳/監査/試行|
|各`independent-fields.jsonl`, `independent-result.json`|入力のみの独立全水理/3D力学求解。別PID|
|各`native-fields.jsonl`, `comparison.json`|対応場、活動ID・全点の差。保存場の派生観測を含む|
|各最細ケースの`restarts/`, `tamper/`|実checkpointと別PID全復旧、外側再hash後の同じ実行版による拒否|
|`supplemental-audit.json`|新規求解済み全状態の既存ゲート、支持反力投影、質量/凍結/解放時計の監査|
|`common-grid-convergence.json`|同じ16物理時刻・補間なしの絶対L2差/観測次数|
|JOINT対照の`case/`|全Run・全50/54場・静水圧準備点・全checkpoint・全後続・改変/世代/開閉|
|拒否群の`refusals.json` / 各`cli.log`|16入力の実CLI終了コードと固有ガード、結果ファイルなし|

`collect.py`は今回の保存フォルダー名を集計する配布証拠用で、独立ソルバーではない。
PCで個々の出力を読むだけなら実行不要。今回と同じ名前を親フォルダーの下へ指定すれば、
`--work 親フォルダー --output 空の集計フォルダー`で集計できる。

### 標準CLIで乾燥JOINTの追加対照

観測器と標準CLIは異なるバイナリなので、同一実行版の復旧一致とは区別する。

```powershell
$Dry = 'runs/NEW-a04-PC-dry/case'
& "$Dry/../bin/solver.exe" solve-stages "$Dry/mesh.json" "$Dry/ground.json" "$Dry/case.json" "$Dry/production-run.json"
if ($LASTEXITCODE -ne 0) { throw '乾燥JOINTの標準CLI対照失敗' }
```

静水圧対照はドライバ内で標準CLIを既に実行する。
標準CLIと観測器の全Run比較は時間とwrapper実行版SHAだけを除く。
checkpoint復旧では実行版SHAも含め、時間以外を一切除かない。

### 同一外水面の水頭整合に残る拒否漏れ

`boundary-consistency-gap.json`の先行入力は、同じfloorへ外水位2.2 mと水理H=1.8/2.1 mを
指定し、本体CLIが完走した。今回定義した整合契約では期待した拒否が未達である。
このケースは限定受入から除外し、16件の成功した入力拒否へ加算しない。
最終ドライバは両値を2.2/2.1 mに一致させる。独立参照だけには限定契約の明示検査があり、
本体のcontroller/schema/既定動作は変えていない。PCでも拒否漏れを合格にしない。

### 合否の扱い

- 全入力の定義済み受理点・全ID/場が揃う。MCの全活動GPが弾性枝である。
  各受理点の既存の変位/応力相対1e-7、全節点力1e-10 kN、幾何/拘束1e-10 mを維持する。
  本体水理残差基準1、独立側の厳しい.01、既存の真の線形基準を緩めない。
- 質量台帳は貯留＋除去水−初期水−正味流入。解放中の水理/時刻は厳密凍結。
  非活動Hを補完しない。保存状態の外側hashを再計算しても内容改変を拒否する。
- 実checkpoint14点、後続268状態、全Runがelapsed_secondsだけを除いて厳密一致する。
  生産CLIのSHAと別PIDを記録し、丸め差の許容化を行わない。
- 合計支持反力は節点成分の基準と別量。生差を記録し、新しい合計量の合格値は作らない。
  JOINTの全接触履歴・Richards全連成の包括許容値も新設しない。
- `numerical_acceptance=false`、A-04 `passed=false`を保持する。静水圧対照の流量/水量/物理秒のnullを
  0へ埋めて合格にしない。時間刻み差の減少を一般空間収束や湿潤JOINTへ外挿しない。

## PC専用の長い時間系列

クラウドではN=4/8/16の小型まで。次は同じ12 Tet4で時間系列だけを長くする。
4秒区間を細分化するもので、物理法則・境界値・力学解放4等分・空間格子は変更しない。

```powershell
$env:A04_PROCESS_BUDGET_SECONDS='21600'
try {
  & $Py "$Study/study.py" --output runs/NEW-a04-PC-fine --steps 64 128 256
  if ($LASTEXITCODE -ne 0) { throw '細分系列は未完了。途中・失敗ログを保持' }
  & $Py "$Study/audit.py" --run runs/NEW-a04-PC-fine
  if ($LASTEXITCODE -ne 0) { throw '細分系列の全点監査失敗' }
  & $Py "$Study/convergence.py" --run runs/NEW-a04-PC-fine --steps 64 128 256
  if ($LASTEXITCODE -ne 0) { throw '同一256時刻の細分比較失敗' }
} finally { Remove-Item Env:A04_PROCESS_BUDGET_SECONDS -ErrorAction SilentlyContinue }
```

各ケース/各実装の受理数は4N+8（264/520/1032）。全3水準・VG/飽和で3632状態対。
最細の保存点は1.5秒・2秒解放途中・3.5秒と同じ。出力の列・合否基準は上と同じ。
PCの未実測予算としてN64は10〜30分、N128は20〜60分、N256は40〜120分、
全体で1〜4時間と数GBの出力枠を準備する。6時間の子プロセス予算で停止しても成功にしない。
密な独立参照の32節点ガードを外して大メッシュを流す手順ではない。

大規模3D・実CAD・GUI/GPU・数時間の物理期間・別の施工履歴はこの契約外。
その入力と独立全場の受入仕様が揃うまでPCでもA-04を受入にしない。
とくにRichards湿潤JOINTや双方向の拒否解除・未定義水頭の補完・新しい材料を行わない。

## Linuxの同じ実行

```bash
source /workspace/soil-environment/activate.sh
cargo build --locked --release
cargo build --locked --release --examples
python interim/2026-11-29-richards-wet-joint-a04/validate.py --output runs/NEW-a04-linux-tests
python interim/2026-11-29-richards-wet-joint-a04/study.py --output runs/NEW-a04-linux-soil --steps 4 8 16
python interim/2026-11-29-richards-wet-joint-a04/audit.py --run runs/NEW-a04-linux-soil
python interim/2026-11-29-richards-wet-joint-a04/convergence.py --run runs/NEW-a04-linux-soil --steps 4 8 16
python interim/2026-11-29-richards-wet-joint-a04/joint_controls.py --kind dry --steps 4 --output runs/NEW-a04-linux-dry
python interim/2026-11-29-richards-wet-joint-a04/joint_controls.py --kind hydrostatic --steps 4 --output runs/NEW-a04-linux-hydrostatic
python interim/2026-11-29-richards-wet-joint-a04/refusals.py --solver target/release/soil-fem-3d --output runs/NEW-a04-linux-refusals
```

各コマンドの終了コードを個別に確認する。全件検証を中断したログを最終成功へ読み替えない。

## 配布された保存証拠のSHA監査

```powershell
& $Py "$Study/package_evidence.py" --verify "$Study/evidence-manifest.json" --scratch runs/NEW-a04-archive-verification
if ($LASTEXITCODE -ne 0) { throw '先行証拠part/ZIP/エントリーSHA不一致' }
& $Py "$Study/package_evidence.py" --verify "$Study/corrected/evidence-manifest.json" --scratch runs/NEW-a04-corrected-archive-verification
if ($LASTEXITCODE -ne 0) { throw '最終の整合境界系列part/ZIP/エントリーSHA不一致' }
```

元の分割証拠（有効なJOINT対照と先行の土系列）とcorrected（最終の整合境界土系列）の両方を監査する。
元のpart/hashは変更せず保持した。各part/ZIPのSHAと全エントリーのサイズ・CRC/SHAを再読込する。これは保存証拠の監査であり、
新規求解・物理受入・別PID復旧の代わりにはならない。Linux実行版はSHAを保存し、配布ZIPに含めない。
Windowsでは自身のビルドで新規求解し、その同一実行版でcheckpointを生成・復旧する。
