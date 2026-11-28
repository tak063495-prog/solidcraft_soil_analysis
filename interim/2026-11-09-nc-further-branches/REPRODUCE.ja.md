# 既存5 kPa NCの再計測・有界延長

対象は `make_free_contact.py --mcc --pc 100 --preload-kpa 5 --force 1 .5 0` の既存小規模ケースだけ。
起点はPR #13、`codex/plastic-spatial-convergence` の `90d24e4181016a1ba8211719ee16b284d82b5eef`。
この追補を含む先端を配置する。PR #13はマージしない。

## 求解しない再計測

既存のPython 3.13環境と `requirements-reference.lock.txt` を使う。BLASのスレッドは1とする。
PowerShellでは以下をリポジトリのルートから実行する。出力先は毎回未存在のフォルダを指定する。

```powershell
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
python interim/2026-11-09-nc-further-branches/remeasure.py --output work/nc-further-remeasured
```

全分割アーカイブのサイズ・SHA256、全収録エントリー、診断ソースのSHAを検証する。
保存された応力から全受理点の力を再組立し、実拡大系・根探索系の真の線形残差、凍結、全場の刻み差、再開一致を再計測する。
材料積分・釣合い求解は行わない。`remeasurement.json` の `exact_saved_measurements=true` を確認する。
値を手で合わせない。ソースや記録が違う場合は両方のJSONを保存する。

## 同じ小規模計算の準備

Linuxでは既存venvを使い、上の環境変数を `export` する。Rustは `rust-toolchain.toml` と既存lockfileを使う。
このクラウドで禁止された大規模ケースを実行するコマンドは含めない。

```powershell
python interim/2026-11-09-nc-beyond-fold/remeasure.py --output work/nc-prepare/beyond-verified
python interim/2026-11-09-nc-further-branches/prepare_seed.py --output work/nc-prepare/prior-seed
python interim/2026-11-09-nc-pseudo-arclength/prepare_evidence.py --output work/nc-prepare/prepared --source-only
cargo build --locked --offline --release --manifest-path work/nc-prepare/prepared/isolated/Cargo.toml --target-dir work/nc-prepare/prepared/target --example material_diagnostic
```

最初のコマンドは前回fold後の202エントリーを検証・再計測し、2番目はその前の4,513エントリー全体を検証して独立起点だけを残す。
同じ69ファイルの隔離観測crateを作り、元の `src` や本番jobは変更しない。
Windowsではバイナリ名に `.exe` を付ける。

**既存Linux checkpointのnative再開には観測バイナリSHAの一致が必須。** 今回のLinux再ビルドでは一致した。
Windowsや異なるコンパイラのバイナリにcheckpointの束縛を付け替えてはならない。
一致しない場合、nativeは `interim/2026-11-09-nc-beyond-fold/reproduce.py` を用いて、そのバイナリ自身で前回の0.1024 m区間とcheckpointを生成し、その出力を以下の親起点に使う。
これはWindowsで今回の再現が確認済みという意味ではない。封印済みの本番jobは書き換えない。
独立側にはRustの場・接線・荷重係数を流用しない。

## 今回の2区間の延長

以下は独立写像の例。1区間目は数値弧長0.1024→0.512 m、2区間目は0.512→0.9216 m。

```powershell
python interim/2026-11-09-nc-beyond-fold/continue_path.py --map independent --reference work/nc-prepare/beyond-verified/evidence/origin/reference --arc work/nc-prepare/prior-seed/independent-arc-history400 --resume work/nc-prepare/beyond-verified/evidence/extended-independent-whole/restart.json --output work/nc-run/extended-independent-long --step-size .0001 --count 4096 --budget-seconds 1800
python interim/2026-11-09-nc-beyond-fold/continue_path.py --map independent --reference work/nc-prepare/beyond-verified/evidence/origin/reference --arc work/nc-prepare/prior-seed/independent-arc-history400 --resume work/nc-run/extended-independent-long/restart.json --output work/nc-run/extended-independent-further --step-size .0001 --count 4096 --budget-seconds 1800
```

nativeでは `--map native --evidence work/nc-prepare/beyond-verified/evidence/origin/prior --binary work/nc-prepare/prepared/target/release/examples/material_diagnostic.exe` とし、native自身のcheckpointを指定する。
Linuxでは `.exe` を除く。既存の力1e-8、真の線形残差1e-10、平面残差1e-12 m、最低1補正、最大20 Newton、ひずみL1上限0.001を変えない。
単一起動は最大4,096点・1,800秒。λの試行範囲は[0,1]。刻みを失敗回避のために調整しない。
各区間の実測・停止理由は `report.json` を読む。終了コード0は合格を意味しない。

## 他枝の局所診断

`branch_path.py` は束縛を照合した保存履歴を保ち、接線予測の向きだけを反転する。
例は独立側で、nativeでは前節の自身の引数・checkpointを使う。

```powershell
python interim/2026-11-09-nc-further-branches/branch_path.py --map independent --reference work/nc-prepare/beyond-verified/evidence/origin/reference --arc work/nc-prepare/prior-seed/independent-arc-history400 --seed work/nc-prepare/beyond-verified/evidence/extended-independent-whole/restart.json --orientation -1 --output work/nc-run/extended-independent-reverse --step-size .0001 --count 1024 --budget-seconds 900
python interim/2026-11-09-nc-further-branches/local_roots.py --map independent --reference work/nc-prepare/beyond-verified/evidence/origin/reference --arc work/nc-prepare/prior-seed/independent-arc-history400 --seed work/nc-prepare/beyond-verified/evidence/extended-independent-whole/restart.json --output work/nc-run/roots-independent-anchor1024 --budget-seconds 600
python interim/2026-11-09-nc-further-branches/local_roots.py --map independent --reference work/nc-prepare/beyond-verified/evidence/origin/reference --arc work/nc-prepare/prior-seed/independent-arc-history400 --seed work/nc-run/extended-independent-long/restart.json --output work/nc-run/roots-independent-anchor5120 --budget-seconds 600
```

各根探索は保存接線・対称化接線の最小2固有ベクトルの3方向、既存7振幅、保存λと1.0の2目標、計42試行。
固有ベクトルは数値的な初期方向だけで、安定性の受入条件ではない。
ひずみ上限・最大20 Newton・既存7段階line searchを維持する。同じ履歴からの有界増分根だけを探すため、見つからない枝は未探索として残る。
今回、非ゼロ増分の候補根は得られなかった。得られた場合も、本体の単調撤去履歴へ飛躍として追加しない。

### ゼロ増分接線からの追加経路

`unloading_path.py` は同じ保存点の実材料写像をゼロ増分で評価し、その実Kに対する `K v = L` を解いて `(v,1)` を正規化する。
Kを対称化して求解せず、真の線形残差1e-10と開始場の力/凍結を検査する。
同じ全材料履歴を保ったまま、正・負の予測方向から既存の移動平面補正を試す。既存接線反転の失敗を刻み縮小や許容値変更で回避する操作ではない。

```powershell
python interim/2026-11-09-nc-further-branches/unloading_path.py --map independent --reference work/nc-prepare/beyond-verified/evidence/origin/reference --arc work/nc-prepare/prior-seed/independent-arc-history400 --seed work/nc-prepare/beyond-verified/evidence/extended-independent-whole/restart.json --orientation -1 --output work/nc-run/extended-independent-unloading --step-size .0001 --count 1024 --budget-seconds 900
```

正方向は `--orientation 1` と別の出力先を使う。nativeでは自身の引数・保存点を指定する。
今回の除荷方向は各1,024点を追跡し、正方向は両写像とも既存ひずみ上限で0点受理・判定不能だった。
数値弧長はこの分岐候補の保存点から別に数える。元の下降経路と加算しない。
除荷側の経路も安定性を認定した代替静的経路ではなく、本体の単調撤去履歴へ追加しない。

## 刻み感度・再開

同じ保存点から `branch_path.py --orientation 1 --step-size .00005 --count 1024 --budget-seconds 900` を実行すると、0.1024→0.1536 mの局所刻み半減になる。
今回この全512対応点の全節点3D変位・324 GP応力/塑性ひずみ/pc/間隙比を比較した。0.9216 mまでの全区間の刻み収束へ拡張しない。
前回の細刻み0.1024 m終端を使う `continue_path.py --step-size .00005 --count 1024` の継続比較も別に保存した。
これは前回区間から累積した離散化差を含む比較なので、同一場起点の比較とは区別する。

nativeの再開試験は0.512 m以降の連続128点と、64点終了時の実checkpointを別プロセスで再開した64+64点を比較した。
独立側は進行中の主経路で実際に書かれたcheckpoint（step 2799）を完全なJSON・ペイロードSHA検証後にコピーし、別プロセスで64点を再開した。
元の連続4096点の対応64行と全保存場・監査行を照合した。終了後に保存行からcheckpointを捏造していない。
いずれも全撤去完了までの本番job再開を証明しない。

## 試験・JSON

```powershell
$env:NC_FURTHER_MEASUREMENTS='work/nc-further-remeasured/recomputed/measurements.json'
$env:NC_FURTHER_EVIDENCE='work/nc-further-remeasured/evidence'
python -m unittest discover -s interim/2026-11-09-nc-further-branches -p test_further.py -v
$env:NC_EXTENDED_MEASUREMENTS='interim/2026-11-09-nc-beyond-fold/measurements.json'
python -m unittest discover -s interim/2026-11-09-nc-beyond-fold -p 'test_*.py' -v
$env:NC_ARC_MEASUREMENTS='interim/2026-11-09-nc-pseudo-arclength/measurements.json'
python -m unittest discover -s interim/2026-11-09-nc-pseudo-arclength -p test_diagnostic.py -v
cargo test --locked --offline --release --test mcc_midpoint --test free_contact
```

`measurements.json` の `passed/history_converged/static_limit_confirmed`、本体確定履歴、各 `paths` の停止理由・力/線形残差・凍結、`local_root_searches` の有限試行範囲、`step_sensitivity` の限定区間、`restarts` の実全行一致を確認する。
全撤去未達、未探索区間・他枝、物理的安定性、全解放の凍結、完了までの再開、Windows新先端CIは未確認のまま。
時間切れ・未収束・カットバック枯渇・点数上限は破壊にも合格にもしない。

クラウド実測の時間は各数値報告と `numerical-processes.json` に残す。同程度CPUで主延長は数分〜十数分、準備・監査・アーカイブを含め十数分〜数十分、PC性能によってさらに長くなる。
再計測は材料求解より軽いが数GBの展開を伴う。今回と前回の全作業領域・生の材料要求/応答まで保持する場合は30 GB程度以上の空きが必要になり得る。
Windows実行は今回未実施。CIはこのinterimの新規Python試験を自動探索しないため、旧CIやLinux成功を新先端CI成功に転用しない。
