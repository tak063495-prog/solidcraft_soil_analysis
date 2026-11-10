# 再監査・新規計算・実checkpointからの続行

PR #13 はマージしない。指定起点 `b235824ea55a0e70b8788e586ee8562549f8baf5` を祖先に持つ同じブランチで作業する。未保存変更をresetしない。

```sh
git fetch origin codex/plastic-spatial-convergence
git switch codex/plastic-spatial-convergence
git pull --ff-only origin codex/plastic-spatial-convergence
git merge-base --is-ancestor b235824ea55a0e70b8788e586ee8562549f8baf5 HEAD
```

初回cloneでは `git -c core.autocrlf=false clone --branch codex/plastic-spatial-convergence https://github.com/tak063495-prog/solidcraft_soil_analysis.git` を使える。大きい既存証拠を含む。今回取得した履歴packは約13 GB、既存interim全体は約14 GBだった。展開とビルドにも空き容量を確保する。疎取得では必要なアーカイブを欠落させたまま試験を成功にしない。

Rust 1.99.1、Python 3.13、既存 `requirements-reference.lock.txt` を使用する。今回の環境は Python 3.13.17、NumPy 2.5.4、SciPy 1.18.1。

```sh
export CARGO_BUILD_JOBS=2 RUST_TEST_THREADS=2
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
export PYTHONUTF8=1 PYTHONIOENCODING=utf-8
```

PowerShellでは各変数を `$env:変数名='値'` で設定する。以下の出力先はすべて未存在のものを選ぶ。

## 保存した生証拠の再監査

新しいMCC釣合い求解は行わない。

```sh
python interim/2026-11-10-nc-postfold-branches/remeasure.py --output work/NEW-postfold-replay
```

全分割、結合全体、全エントリーのSHAを検査して安全に展開し、封印した `drivers/measure.py` を `-I -X utf8` の別プロセスで実行する。全受理材料場からの力再組立、元Hの成分境界、撤去済み材料・接触履歴、受理／棄却の区別、弧長経路の実拡大線形系を調べる。新しい釣合い解や独立解を保存場から生成するものではない。

`verification.json` の `archive_verified=true`、`measurements_byte_identical=true`、`numerical_acceptance=false` を確認する。派生浮動小数の厳密なバイト一致は記録環境と同じLinux/WSLを推奨する。WindowsのBLAS等で最下位桁が変われば差を隠さず失敗する。アーカイブのSHA一致と数値合格は別である。

旧入力の全照合結果は `prepared/prior/extraction.json` に入る。新規計算では次節のprepareが旧2,034,108,934バイト・122分割・1,349エントリーをすべて照合する。今回の再監査器は新アーカイブの生証拠を再組立するもので、旧アーカイブの再取得や旧前史の再求解を行ったとはしない。

## 同じ小さいケースを新規に計算する

```sh
cargo build --locked --release
python interim/2026-11-10-nc-postfold-branches/prepare.py --output work/NEW-postfold-prepared
python interim/2026-11-10-nc-postfold-branches/reproduce.py --evidence work/NEW-postfold-prepared/prior --binary work/NEW-postfold-prepared/observer/target/release/examples/material_diagnostic --output work/NEW-postfold-study
python interim/2026-11-10-nc-postfold-branches/extend.py --evidence work/NEW-postfold-prepared/prior --binary work/NEW-postfold-prepared/observer/target/release/examples/material_diagnostic --study work/NEW-postfold-study
python interim/2026-11-10-nc-postfold-branches/measure.py --source work/NEW-postfold-study --output work/NEW-postfold-measurements.json
```

Windowsの観測器には `.exe` を付ける。prepareは現在のRustソースを隔離コピーし、既存の読み取り専用観測コードを追加する。生産srcは変更しない。旧ヘルパーが保存する過去のコミットラベルと、今回の実際の呼出しコミットは新しい `preparation.json` で区別する。

入力は既存 `make_free_contact.py --mcc --pc 100 --preload-kpa 5 --force 1 .5 0` の54節点・12 Tet10・324材料GP・32接触点。pcの100 kPaと追加プレロード5 kPaを混同しない。本体は自身の元Hを使用し、独立側は入力のみから得た自身の400分割前史とHを使用する。Rust場を独立解の生成へ代入しない。

reproduceは最大2子プロセスで、長い弧長2経路、別起点1経路、弱いモードの正負6予測、Cartesian成分順序6経路、選んだ2順序の弧長、実middle checkpointからの部分再開を実行する。extendはz成分後の折返し後の追加長2,048増分と、最大s付近の正負2予測を各128増分実行する。部分再開を含む合計21実行を記録した。各経路の点数・時間予算は `processes.json` と `controls.json` に残る。PC速度で時間予算停止位置は変わる。終了コード0は診断処理の終了であり、`passed=true` ではない。途中の予算停止を破壊・静的限界・他経路の不存在にしない。

弧長では既存の接線超平面、1 mの荷重座標尺度、同一行列の反復改良を使う。力1e−8、真の線形残差1e−10、補正1e−11 m、拘束1e−12 m、増分適合ひずみL1≤0.001、Newton20回を維持する。元の保持荷重を超える成分、符号反転、元のゼロ成分への追加を拒否する。各経路は現在の区分線形荷重区間の外へ外挿しない。全解放端点を横切る候補には既存荷重制御による端点補正を用意したが、今回そこへ到達したとはしない。

sは各順序に固有の経路座標であり、異なる順序の同じsが同じ荷重を意味するとは限らない。s=1かつ保持荷重ゼロだけが全解放λ=1である。折返しではHの一部が増えるので、単調撤去の認定ではない。

## 実checkpointから別プロセスで再開

新規計算の本体長経路の実middleを使う例。

```sh
python interim/2026-11-10-nc-postfold-branches/trace.py --evidence work/NEW-postfold-prepared/prior --binary work/NEW-postfold-prepared/observer/target/release/examples/material_diagnostic --output work/NEW-postfold-resumed --map native --step-size .0002 --count 2048 --budget-seconds 600 --resume work/NEW-postfold-study/native-long/middle.json
```

driver、実行版、起点、数値設定とペイロードSHAの束縛を維持する。Linuxのcheckpointの実行版SHAをWindowsのものへ付け替えない。各PCで自身の観測器と新規経路を作る。これは今回の部分終端までの再開であり、未到達のλ=1までの再開検証ではない。

順序変更の予算停止後も、実restartを変更せず追加予算で続けられる。

```sh
python interim/2026-11-10-nc-postfold-branches/reorder.py --evidence work/NEW-postfold-prepared/prior --binary work/NEW-postfold-prepared/observer/target/release/examples/material_diagnostic --output work/NEW-x-first-continued --map native --order x-first --budget-seconds 1800 --resume work/NEW-postfold-study/native-x-first/restart.json
```

別の受理場からの新しい診断起点は `trace.py --anchor-samples .../samples.jsonl.gz --anchor-step N` で指定できる。弱いモードの予測変更には、その受理増分の実接線を保存した弧長samplesと `--probe -1`、`-.25`、`.25`、`1` を使う。新起点と実checkpointの再開は区別する。旧checkpointを編集・再封印しない。

## Linuxの全件検証

```sh
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo test --locked --release
cargo build --locked --release
python tools/run_ci_python_tests.py --group python-core-0 --output work/NEW-core-0
python tools/run_ci_python_tests.py --group python-core-1 --output work/NEW-core-1
python tools/run_ci_python_tests.py --group python-mean-evidence --output work/NEW-mean
python tools/run_ci_python_tests.py --group python-spatial-evidence --output work/NEW-spatial
python -m unittest discover -s tools -p 'test_*.py' -v
```

4群すべてを実行し、各inventoryの全発見IDが一致し、全割当IDの和集合に欠落・重複がないことを確認する。toolsの2件も4群の割当に含む。単独のtools実行と二重に合格件数を水増ししない。欠落・スキップ・失敗を成功にしない。

Windowsの起点CI失敗はPOSIX専用2試験のスキップによるものだった。POSIXでは従来の親だけへの外部シグナルを保持し、Windowsでは同じ実CLIの親メインスレッドへ実シグナルを届ける。子は別プロセスで動き、子終了確認、終了確認後のロック解放、安全な再実行を検査する。Windowsの強制終了をPOSIXの協調的SIGTERMと同一視しない。ランナーのスキップ拒否、assert、数値許容値を緩めていない。

## 合格へ不足する証拠とJOINT

全解放候補が得られても、撤去中の全材料・接触履歴、残存域λ=0の実Newton釣合い、全解放後の実ゼロ増分による全状態凍結、入力だけからの独立塑性3D完了全場、両実装の完了までの刻み感度、実checkpointから完了までの別プロセス再開をすべて揃えるまで合格にしない。今回の部分経路・部分再開で代用しない。

JOINTの部分亀裂先端は利用者の形状・先端位置・共有節点方針、接触BIRTHは既存履歴と設置参照・水理初期値、高次板／弱結合は補間・仮想仕事と独立3D証拠が必要。既存小ケースから新仕様を推測しない。今回のNC調査をJOINT全履歴収束、一般塑性・部材・水の独立3D統合へ転用しない。証拠不足は未検証・拒否またはユーザーPC手順として残す。P3b/P4各行、構造・水・JOINT全体、P6全体は未完である。
