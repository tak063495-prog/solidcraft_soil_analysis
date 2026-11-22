# C-01・C-03 非線形制御の再実行

PR #13 の `codex/plastic-spatial-convergence` で、リポジトリのルートから実行する。既存の失敗出力は保存する。以下は新しい出力先を使う Windows PowerShell の例。PC の比較は未実施。

## 準備

```powershell
git rev-parse HEAD
$env:CARGO_BUILD_JOBS='2'
$env:RUST_TEST_THREADS='2'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
python -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE) { throw '依存関係の導入失敗' }
cargo build --locked --release --example wet_joint_path_probe --example c03_bulk_cause_probe
if ($LASTEXITCODE) { throw 'ビルド失敗' }
$S1='interim/2026-12-12-c01-switch-event-localization/study.py'
$S3='interim/2026-12-18-c03-plastic-bulk-cause-diagnosis/c03diag_pc.py'
```

## C-01 native / independent、512・1024

12 Tet10、54節点、材料324点、接触32点を使用する。`--event-proposal` は符号付きギャップによる再試行位置、基準格子までの次の刻み、接触構成則の最少4内分割を選ぶ。各時間水準を独立に初期化する。本体・独立とも元の入力から開始し、互いの受理列は渡さない。

```powershell
foreach ($n in 512,1024) {
  $out="work/NEW-PC-c01-event4-$n"
  python $S1 --mode init --event-proposal --solver target/release/examples/wet_joint_path_probe.exe --output $out
  if ($LASTEXITCODE) { throw "C-01 初期化失敗 $n" }
  python $S1 --mode native --event-proposal --pc-series --steps $n --timeout 21600 --output $out
  if ($LASTEXITCODE) { throw "C-01 native 未完了 $n。途中記録を保存" }
  python -c "import subprocess,sys; sys.exit(subprocess.run([sys.executable,*sys.argv[1:]],timeout=21600).returncode)" $S1 --mode independent --event-proposal --pc-series --steps $n --output $out
  if ($LASTEXITCODE) { throw "C-01 independent 未完了 $n。途中記録を保存" }
}
```

旧経路を比較する場合は別の空の出力先にして上記のすべての `--event-proposal` を外す。旧256は既存の計測を保持する。`--legacy` は別の経路で、今回の改善前対照には使わない。native と independent は別々に実行できる。一方が停止したらその失敗を保存し、もう一方は同じ水準・出力先の対応する未実行コマンドで検査できる。

`gN/native-audit.json` と `gN/independent-audit.json` の受理点数、全7工程の各 `i/N`、保存場を再組立した力残差（本体1e-8、独立1e-10）、粗い解と二分割解のすべり・散逸の差を確認する。失敗時は `native.log`、`native-raw/failed-run.json`、`independent-attempts.jsonl`、`independent-failure.json` の存在分を残す。打切り・時間切れ・全工程完了を区別する。

## C-03 stage 2 independent、g32-m0-c0.0001

既存の12 Tet10 native と対になる入力を作る。native の計算結果を独立求解へ与える操作はない。`--nonlinear-cutback` は独立側の線探索失敗または全体 Newton 上限でだけ再試行する。材料・接触の内分割上限を従来と同じ値に固定する。

```powershell
$P='work/NEW-PC-c03-cutback-g32'
$tag='g32-m0-c0.0001'
python $S3 --user-pc --mode init --refine 0 --solver target/release/examples/c03_bulk_cause_probe.exe --output $P
if ($LASTEXITCODE) { throw 'C-03 初期化失敗' }
python $S3 --user-pc --mode independent --nonlinear-cutback --output $P --folder "$tag-independent" --steps 32 --material 0.0001 --contact 0.0001 --timeout 14400
if ($LASTEXITCODE) { throw 'C-03 独立参照未完了。driver.log と途中場を保存' }
```

同一入力で完了した既存 native の `mesh.json`、`ground.json`、各実行フォルダの `stages.json` を新しいものと照合してから比較する。共通基準点224個の並びは `fields.jsonl`、受理した補助刻みは `substeps.jsonl`、棄却理由を含む試行列は `attempts.jsonl`、全受理点の材料・接触の分割数と残差は `result.json` に保存する。追加の平衡点は経路に影響し得る。固定32分割の切戻しなしという旧判定を補助刻みを使った経路へ引き継がない。比較コマンドは次のように既存の native フォルダを同じ出力の直下にコピーしてから実行する（コピー元は実際のパスに合わせる）。

```powershell
$Old='work/PC-C03-tet10'
foreach ($f in 'mesh.json','ground.json') {
  if ((Get-FileHash "$Old/$f").Hash -ne (Get-FileHash "$P/$f").Hash) { throw "入力相違 $f" }
}
if ((Get-FileHash "$Old/$tag-native/stages.json").Hash -ne (Get-FileHash "$P/$tag-independent/stages.json").Hash) { throw '工程入力相違' }
Copy-Item -Recurse "$Old/$tag-native" "$P/$tag-native"
python $S3 --user-pc --mode compare --output $P --folder "$tag-native" --other "$tag-independent"
if ($LASTEXITCODE) { throw 'C-03 比較未完了' }
```

## 資源見積もり

全コマンドを逐次実行する。対象は12要素のみ。C-01 は PC 実測の旧256・4,341受理・965秒を基準に、7工程×最大2,000試行、各試行の粗・中点・細解、増える接触内積分とファイル入出力を考慮し、各経路に6時間を割り当てた。native 子プロセスと independent 全プロセスはそれぞれ6時間で打ち切る。C-03 は既存12 Tet10の計画値0.5〜4時間を用い、独立ワーカーの上限は4時間。24時間以内の計算計画で、完走時間の保証はない。

メモリの計画値は C-01 16〜32 GB、C-03 32 GB 以下。C-01 は同時1経路、逐次保存、試行記録最大14,000件で見積もった。C-03 は162変位自由度で、自由度削減前の162×162の倍精度行列でも1枚約0.21 MB。材料324点と数値接線、全履歴の入出力を含め、旧12要素手順の32 GB枠を用いる。空きディスクは C-01 1水準30〜60 GB、C-03 20 GBが目安。メモリ実測が40 GBに近づいたら経路を停止して途中記録を残す。96、768、6144要素への細分化は行わない。

## クラウド比較の再現

Linux では拡張子なしの実行ファイルを使う。独立参照は既存と同じ Python 3.13・固定版パッケージ、BLAS は1スレッドで実行する。probe の出力名は未使用のものを用意する。

```sh
T=interim/2026-11-22-impl-c01-c03-nonlinear
O=work/NEW-cloud-c01-c03
mkdir -p "$O"
for n in 512 1024; do
  python "$T/probe_c01.py" --steps "$n" --output "$O/independent-$n-old.json" || exit 1
  python "$T/probe_c01.py" --steps "$n" --event-proposal --local-substeps 1 --output "$O/independent-$n-event.json" || exit 1
  python "$T/probe_c01.py" --steps "$n" --event-proposal --local-substeps 4 --output "$O/independent-$n-new.json" || exit 1
done
python "$T/probe_c03.py" --output "$O/c03-tet4-old.json"
python "$T/probe_c03.py" --cutback --output "$O/c03-tet4-new.json"
```

本体の観測用 `joint_path_window` は実 solver の Observer で指定位置の受理直後に停止する。既存の試行数を継承した単一の求解で、チェックポイントの書換えは行わない。同じ入力生成器で当該水準の工程 JSON を書き出す。以下は N512 の例で、N1024 は `open 0.04296875` に置き換える。

```sh
cargo build --locked --release --example joint_path_window
python "$T/probe_c01.py" --steps 512 --event-proposal --local-substeps 4 --native-input --output "$O/native-512-stages.json"
target/release/examples/joint_path_window examples/joint-water-column-tet10-mesh.json examples/ground-joint-water-mc.json "$O/native-512-stages.json" reload 0.921875 "$O/native-512.json"
```

`checkpoint_requested` のときは該当工程の最終受理座標を確認する。`stage attempt budget exceeded` は未達。C-03 の同じ12 Tet10 の参照を初期状態から検査する場合は、両方の probe コマンドに `--tet10` を加える。

C-01・C-03 は未達、passed=false。
