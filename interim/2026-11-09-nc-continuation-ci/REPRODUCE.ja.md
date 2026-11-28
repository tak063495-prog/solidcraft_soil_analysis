# 5 kPa NCの今回の継続・λ=0・再計測手順

PR #13の `codex/plastic-spatial-convergence` の本追補を含む先端を配置する。起点は42473341でありmainではない。元の形状、材料、荷重を変えない。100 kPaプレロードの合格は代用しない。
Python 3.13と `requirements-reference.lock.txt`、Rustは `rust-toolchain.toml` とlockfileを使う。
以下はPowerShell、リポジトリルートから実行する。各出力先は未存在にする。

## 材料を解かない全証拠再計測

```powershell
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
python interim/2026-11-09-nc-continuation-ci/stream_remeasure.py --output work/nc-remeasured
$env:NC_CONTINUATION_EVIDENCE='work/nc-remeasured/evidence'
$env:NC_CONTINUATION_MEASUREMENTS='work/nc-remeasured/recomputed/measurements.json'
python -m unittest discover -s interim/2026-11-09-nc-continuation-ci -p test_evidence.py -v
```

結合一時ファイルを作らない `stream_remeasure.py` で独立した別ファイルへ実際に展開する（hardlink不使用）。従来と同じ結合方式の `remeasure.py` も残す。分割アーカイブの全サイズ・SHA256、全エントリー・診断ソースを照合し、全受理点の応力から力を再組立する。実Newton/拡大系の残差、全GPの凍結、全対応点の刻み差、実checkpoint再開の全行一致、λ=0解を再計測する。
`remeasurement.json` の `exact_saved_measurements=true` を確認する。材料積分・釣合い求解は実行しない。差を手で合わせない。

## 元の実checkpointから同じ区間を求解

```powershell
python interim/2026-11-09-nc-continuation-ci/prepare_seed.py --output work/nc-seed
cargo build --locked --offline --release --manifest-path work/nc-seed/2026-11-09-nc-further-branches/observer-source/Cargo.toml --target-dir work/nc-observer --example material_diagnostic
```

prepareは過去3アーカイブの全5,067エントリーを検証し、必要な原入力・実checkpoint・観測69ファイルだけを取り出す。過去の数値履歴を再求解したとはしない。
通常のcloneには分割ファイルが含まれる。sparse checkoutを使う場合は対象3フォルダの `evidence.tar.xz.part*` も取得する。
依存がローカルにない初回は、lockfileを維持して通常の `cargo fetch --locked` で準備してからoffline buildする。

独立塑性3D写像の粗刻みは以下。Rustの場・応力・接線・追跡係数を独立場へ渡さない。

```powershell
python interim/2026-11-09-nc-beyond-fold/continue_path.py --map independent --reference work/nc-seed/2026-11-09-nc-beyond-fold/origin/reference --arc work/nc-seed/2026-11-09-nc-pseudo-arclength/independent-arc-history400 --resume work/nc-seed/2026-11-09-nc-further-branches/extended-independent-further/restart.json --output work/nc-run/extended-independent-next --step-size .0001 --count 4096 --budget-seconds 1800
```

nativeは `--map native --evidence work/nc-seed/2026-11-09-nc-beyond-fold/origin/prior --binary work/nc-observer/release/examples/material_diagnostic.exe` とし、native自身のcheckpointを使う。Linuxでは `.exe` を除く。
**観測バイナリSHA・モデル・数値設定・ソースの束縛が異なれば保存checkpointの再開は拒否する。束縛を書き換えない。**
今回のLinuxでは同じバイナリでそのまま再開できた。Windowsや異なるコンパイラで同じSHAを作れるとはしない。
不一致なら、以前の `2026-11-09-nc-beyond-fold/REPRODUCE.ja.md` と `reproduce.py` でそのバイナリ自身の0.1024 m checkpointを生成し、`2026-11-09-nc-further-branches/REPRODUCE.ja.md` の2区間を同じバイナリで走らせて0.9216 mを作り直す。そこから今回の区間へ進む。これらは診断経路であり、本番jobの履歴へ追加しない。

## 追加0.4096 m全区間の刻み半減

独立側は以下の2起動で8,192点。nativeも自身の引数・自身のseedで同じ2起動を行う。

```powershell
python interim/2026-11-09-nc-further-branches/branch_path.py --map independent --reference work/nc-seed/2026-11-09-nc-beyond-fold/origin/reference --arc work/nc-seed/2026-11-09-nc-pseudo-arclength/independent-arc-history400 --seed work/nc-seed/2026-11-09-nc-further-branches/extended-independent-further/restart.json --orientation 1 --output work/nc-run/extended-independent-half-a --step-size .00005 --count 4096 --budget-seconds 1800
python interim/2026-11-09-nc-further-branches/branch_path.py --map independent --reference work/nc-seed/2026-11-09-nc-beyond-fold/origin/reference --arc work/nc-seed/2026-11-09-nc-pseudo-arclength/independent-arc-history400 --seed work/nc-seed/2026-11-09-nc-further-branches/extended-independent-further/restart.json --resume work/nc-run/extended-independent-half-a/restart.json --orientation 1 --output work/nc-run/extended-independent-half-b --step-size .00005 --count 4096 --budget-seconds 1800
```

これは元checkpointを親として明示した別の数値刻み経路で、親の束縛や履歴は変更しない。fineの弧長カウンターは0から、粗刻みは継承した0.9216 mから始まる。比較は両者の追加弧長を合わせる。
最低1補正・最大20 Newton、力1e-8、真の線形残差1e-10、平面1e-12 m、ひずみL1上限0.001、RK4材料分割4・中央差分幅1e-8を維持する。
追加区間の刻み差が小さくても、それ以前の前史・材料分割・空間・全履歴の収束は認定しない。

## 残存域λ=0の非ゼロ初期推定からの釣合い

```powershell
python interim/2026-11-09-nc-continuation-ci/zero_equilibrium.py --map independent --reference work/nc-seed/2026-11-09-nc-beyond-fold/origin/reference --arc work/nc-seed/2026-11-09-nc-pseudo-arclength/independent-arc-history400 --output work/nc-run/zero-independent
```

nativeは前節と同じnative引数と別出力先を使う。切断前の各自の場・元の保持荷重を固定し、既存根探索の1e-4 m振幅で初期推定を摂動する。全履歴を固定した材料積分と実KによるNewtonでλ=0を解く。
力・真の線形残差・修正量1e-11 mは既存根探索と同じ。非収束なら未達として保存する。λ=1の全解放・全撤去の証明ではない。

## 途中再開

計算中に実際に出力された完全な `restart.json` を別名でコピーし、元バイトのSHA、payload SHA、step、時刻と「元実行がまだ終了していない」ことを記録する。書込み途中のJSONやpayload不一致は採らない。
同じ `continue_path.py`・同じ引数と `--resume <コピー> --count 128` を別プロセスで実行し、主経路の同じstep 128行と全値を照合する。
今回の実コピーは証拠中の `captured-native.json` と `captured-independent.json`、来歴は `checkpoint-captures.json`。
保存行からcheckpointを後で組み立てない。本番jobの全撤去完了までの再開は未確認。

## 判定・必要資源・CI

終了コード0は合格ではない。全上部撤去λ=1、全解放後の凍結、単調撤去の全履歴収束、独立3D全撤去、本番job完了までの再開、静的限界/安定性は未達・判定不能のまま。
時間切れ、Newton上限、点数上限を破壊や合格とみなさない。許容値変更・ピボット・カットバックの調整で通さない。
既存の小さい54節点・12 Tet10だけを対象とする。新しい形状・参照値・合格値はない。

全場JSONLとnative試行IOのため、準備・計算・再展開を同時に保持する場合は40 GB程度以上の空きを推奨する。上記の単一起動は最大1,800秒で停止し、遅いPCでは不完了になる。成功のために予算境界を受入条件へ読み替えない。
今回の長い経路は従来どおりCIの既定外。Windows CIの全標準試験・既存小規模検証は削除していない。[CI手順](CI.ja.md)を参照。
