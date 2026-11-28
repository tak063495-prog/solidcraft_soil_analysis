# A 系の小型検証と NC-01 T9(c) の PC 照合手順

## Linux の小型検証

リポジトリのルートで、固定の Python 3.13、requirements-reference.lock.txt、rust-toolchain.toml を使う。出力は未使用のディレクトリを指定する。下記は bash 用。各手順はエラーの時点で止める。

```bash
set -eu
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
D=interim/2026-11-23-impl-root-causes-a-nc
O=runs/root-causes-a-reproduction
mkdir -p "$O"
cargo build --locked --release --bin soil-fem-3d --example t7_boundary
cli=target/release/soil-fem-3d
bridge=target/release/examples/t7_boundary
python "$D/a02.py" "$cli" "$O/a02"
python "$D/a04.py" "$cli" "$O/a04"
python "$D/a04.py" "$cli" "$O/a04-decimal40" --reference decimal40
for n in 18 144; do
  folder="$O/mesh-$n"
  python "$D/branch.py" prepare --elements "$n" "$cli" "$folder" "$bridge"
  for engine in native independent; do
    python "$D/branch.py" follow-same-contact "$engine" "$folder/$engine.gz" "$folder" "$bridge"
    part="$folder/$engine-4.gz"
    python "$D/branch.py" follow-same-contact "$engine" "$part" "$folder" "$bridge" --steps 4
    read -r digest filename < <(sha256sum "$part")
    python "$D/branch.py" follow-same-contact "$engine" "$folder/$engine-resumed.gz" "$folder" "$bridge" --resume "$part" "$digest"
  done
  python "$D/compare_branch.py" "$folder" "$bridge"
done
```

A-01 の prehistory.gz は受理状態、origins.gz は元入力・前史・版の束縛、各エンジンの gzip は open の試行・受理全場・仕事・再開用接線を含む。再開には同じソース・入力・評価実行版・数値環境と親セグメントの SHA-256 を使う。失敗した試行は親状態のまま。直接経路と分割経路の全受理点、試行列、終端チェックポイントの完全一致を compare_branch.py で確認する。スカラー表示が近い二点を同一座標とみなさない。

A-02 の出力 run.json と results.json に全節点力・摩擦円の照合結果がある。A-04 の native-fields.jsonl と independent-fields.jsonl には実時刻・施工進捗の付いた49場がある。全界面対照は T8 の water.py の affine_column を AFFINE_SLOPES × AFFINE_STEPS で実行し、assess.affine_water_passed で照合する。曲面は `cargo test --release --locked --lib curved_aperture_storage` で確認する。[T8 の照合実装](../2026-11-22-impl-t8-acceptance-criteria/assess.py)

## NC-01 PC-3　実行準備のみ

**今回の PC 実行は0件。** PC-3 は元の 5 kPa 履歴の **12 Tet10・54節点** で、入力から独立に前史を積分し、長い弧長経路と再開を調べる。96 Tet10 のクラウド／PC の線形停止は別の行に保存する。旧 PC-2 の独立0点、Linux の短い8点、別施工 NC-01b の完走は代用できない。

### 版と入力バイト

数値再生用の版は `e4e2dfd5c7bc189ccb471ab9cfe8621792e7b22f` に固定する。PR の最新 SHA が異なってもこの再生版を動かす。新しいディレクトリに detached worktree を作り、PC-2 の6本の native 監査があるディレクトリと Python 3.13 の実行版は絶対パスで渡す。PC-2 の元ファイルは読み取りのみ。依存はこの版の requirements-reference.lock.txt に従う。

```powershell
$pin = 'e4e2dfd5c7bc189ccb471ab9cfe8621792e7b22f'
git -c core.autocrlf=false worktree add --detach ../nc-pc3-pinned $pin
if ($LASTEXITCODE -ne 0) { throw 'worktree の作成失敗' }
Set-Location ../nc-pc3-pinned
if ((git rev-parse HEAD) -ne $pin -or (git status --porcelain --untracked-files=no)) { throw '固定版または追跡ファイルが不一致' }
$Py = (Resolve-Path (Read-Host '依存導入済み python.exe の絶対パス')).Path
$Pc2 = (Resolve-Path (Read-Host '既存 PC-2 native 監査の絶対パス')).Path
$Out = 'runs/new-pc3-pinned'
$Pre = 'runs/new-pc3-input-pin'
New-Item -ItemType Directory -Path $Pre -ErrorAction Stop | Out-Null
$env:PYTHONPATH = (Resolve-Path tools).Path
& $Py -c "from pathlib import Path; import sys; from probe_low_preload_nc import NAMES,expected_inputs,write; [write(Path(sys.argv[1])/(n+'.json'),v) for n,v in zip(NAMES,expected_inputs(),strict=True)]" $Pre
if ($LASTEXITCODE -ne 0) { throw '入力生成失敗' }
$expected = @{
    mesh   = '2e4190d4fda3a4304170ce04743b265e6de57a1081f8a2809b8db913c7ed8ebf'
    ground = 'e83e98e9481920b0a91484c9489ade6d5308493e228f1736439cce99193e6abf'
    stages = '14be5443f212cdc987599e86cf536f55a788a5086d6ffa11f0f175050208a400'
}
function Confirm-Inputs($folder) {
    foreach ($n in $expected.Keys) {
        $path = Join-Path $folder ($n + '.json')
        $actual = (Get-FileHash -LiteralPath $path -Algorithm SHA256 -ErrorAction Stop).Hash.ToLowerInvariant()
        if ($actual -ne $expected[$n]) { throw "入力バイト不一致 $path $actual" }
    }
}
Confirm-Inputs $Pre
```

上の SHA は固定ジェネレータが生成する UTF-8・BOM なしの JSON を Windows の書出し改行 CRLF で直列化して算出した期待値。PC の実測値は未取得。比較は各ファイルの生バイトで行う。不一致なら実ファイル・SHA・実行環境を保存し、入力や期待値を書き換えずに止める。

### PC で行う独立比較

以下は PC 実行時のコマンド。既存の runner は前史400/800（RK4）、各々の rate arc 128点、材料内4（ds=2e−5で64点、1e−5で128点）、材料内8（ds=2e−5で64点）を順次実行する。各前史で直接64点と別プロセス32＋32点を照合する。前史は各エンジン自身の入力から生成する。native の応力・保持荷重・変位を独立側に渡さない。

```powershell
& ./interim/2026-12-05-nc-independent-arc-residual/rerun-pc.ps1 -Python $Py -Output $Out -Pc2 $Pc2
if ($LASTEXITCODE -ne 0) { throw '途中成果と停止理由を保存して終了' }
Confirm-Inputs "$Out/reference400"
Confirm-Inputs "$Out/reference800"
Get-FileHash -Algorithm SHA256 -LiteralPath $Py,target/release/soil-fem-3d.exe
```

ヘッダ head.txt、ソース／入力／バイナリ SHA、Python・NumPy・SciPy・BLAS、全コマンドと rc、メモリと壁時計時間を添付する。途中停止時も生成済みの reference400/800 の入力 SHA を確認する。前史の最後が excavate_upper、rate arc の実際の受理128点、停止理由 null、元の拡大系真残差 ≤1e−10 を runner のゲートで確認してから局所経路を読む。全節点変位・保持荷重・反力、全GPの応力・pc・間隙比・接触履歴、仕事、左右弱モード、荷重方向・進捗、再開後の全受理行／実線形系／チェックポイントを同じ前史・座標で照合する。未設定の全場幅は保留する。[PC runner](../2026-12-05-nc-independent-arc-residual/rerun-pc.ps1)・[再開・比較の合格条件](../2026-11-17-nc-fold-stability/REPRODUCE.ja.md)

**資源見積もり。** 54節点×3並進で、162次の密行列は1枚約0.21 MB。履歴・Python・BLAS・OS を含む計画枠を8 GBと見積もる。PC実測ではない。前史5–120分、rate arc 最大180秒、各局所経路1–10分（既存上限1800秒）。1前史あたり局所5呼出しを全部上限まで使っても、求解予算の合計は約4.6時間。このほかビルド・監査・I/Oを含める。経路は直列、子プロセス合計を含めて RAM 40 GB・1経路24時間以内。上限前に保存・停止し、再開時も時間を合算する。768／6144要素は対象外。[既存の規模と時間見積もり](../2026-12-05-nc-independent-arc-residual/REPRODUCE.ja.md)

| 経路 | コミット・実行版SHA／入力3 SHA | 受理数 | 最終受理 λ | 次の目標 λ | 元系真残差 | 全場・再開・判定 |
|---|---|---:|---|---|---|---|
| 既報クラウド96、明示補正 | 過去PCとの同一性未確認 | 190 | 0.5647788047790527 | 0.5647790431976318 | 1.08826e−10 | 棄却 |
| 同日PC96、明示補正（ユーザー転記） | 過去クラウドとの同一性未確認 | 187 | 約0.56478 | 原票で照合待ち | 1.112e−10 | 棄却 |
| 新規PC-3 12、前史400・各分割 | 実行時に記入 | 未実行 | — | — | — | 未判定 |
| 新規PC-3 12、前史800・各分割 | 実行時に記入 | 未実行 | — | — | — | 未判定 |

元入力の受入には5 kPa 履歴で全撤去 λ=1、元系真残差 ≤1e−10 が必要。12 の局所 fold 約0.95140 と固定強度の下界0.95191・上側値0.95227は異なる問題。履歴依存の物理極限は判定保留とする。
