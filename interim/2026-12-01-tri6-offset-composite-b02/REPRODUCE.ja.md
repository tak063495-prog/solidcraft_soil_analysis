# B-02の再現手順（Windows PC / Linux）

保存先の2026-12-01は依頼ラベル。実施日は2026-11-12（日本時間）。
PR #13と参照draft PRはマージしない。以下はリポジトリルートで実行する。
出力先は毎回未存在の名前にする。非ゼロ終了・未収束・時間切れ・取消は成功ではない。

## 固定版・環境・全件試験

```powershell
git fetch origin codex/plastic-spatial-convergence
git worktree add --detach D:\Work\b02-repro origin/codex/plastic-spatial-convergence
Set-Location D:\Work\b02-repro
git rev-parse HEAD
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
python -m venv .venv
.venv\Scripts\python -m pip install -r requirements-reference.lock.txt
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:CARGO_BUILD_JOBS='2'
$env:RUST_TEST_THREADS='2'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
$B='interim/2026-12-01-tri6-offset-composite-b02'
.venv\Scripts\python "$B/validate.py" --part rust --output runs/NEW-b02-validation-rust
.venv\Scripts\python "$B/validate.py" --part python --output runs/NEW-b02-validation-python
```

`validate.py --part rust`はfmt、locked Clippy全ターゲット・warnings-as-errors、
locked release build、Rust全ターゲット試験、doc試験を順に完了させる。
Pythonは既存の発見器を変えず、通常2群と旧NC mean/spatialの2群を全て実行する。
各コマンドのprocess JSONとログ、全ID inventoryが出力される。想定20〜120分、
初回コンパイル・ディスク・CPUで変わる未実測のPC目安。

正常の件数は本フォルダーの`validation.json`と照合する。実行前に既存の全試験証拠を
Gitの元バイトで配置する。疎checkoutでは対象ZIP/分割tarを省略しない。
通常の完全checkoutが前提で、必要なら`git sparse-checkout disable`を行う。
旧NCのLinux厳密検査で準備エラーとなった4+4件は実行済みに数えない。
`allow_platform_float_differences`、skip、参照・hash・許容値の変更で通さない。
Windowsの既存試験が持つOS分岐も変更しない。過去の8件未実行記録は保持する。

## 小〜中規模の新規求解・全受理点・再開

```powershell
.venv\Scripts\python "$B/build_observer.py" work/NEW-b02-observer-build
.venv\Scripts\python "$B/study.py" runs/NEW-b02-formal --binary target/release/b02-observer.exe
.venv\Scripts\python "$B/convergence.py" runs/NEW-b02-formal runs/NEW-b02-convergence.json
.venv\Scripts\python "$B/work_audit.py" runs/NEW-b02-formal runs/NEW-b02-work.json
.venv\Scripts\python "$B/scope.py" runs/NEW-b02-scope --cli target/release/soil-fem-3d.exe
.venv\Scripts\python "$B/mutations.py" runs/NEW-b02-formal runs/NEW-b02-mutations
```

診断実行版は`work/`の別Cargo manifestでビルドし、依存は生産Cargo.lockの全packageを
保持して診断rootだけを追加する。`--locked`で検査し、ルートmanifest/lockを変更しない。
生産ソースへ観測機能を追加せず、既存公開APIだけを呼ぶ。

`study.py`は10条件（板だけ／梁併用×正負ゼロの6条件、空間追加2、荷重刻み追加2）。
12/96/768 Tet10、45/225/1377節点、2/8/32枚の板、刻み4/8/16。
本体の観測経路10、同じ実行版の観測なし対照10、入力だけからの別Python PIDの独立10経路。
各経路の`native/accepted.jsonl`と`independent.jsonl`は新規求解で生成する。
`comparison.json`はその両方の全受理点・全ID・全成分の照合であり、保存JSON再監査の別名ではない。
想定15〜120分、RAM 4〜16 GB、出力空き5 GB以上を目安とする（PC未実測）。
大きいジャーナルは展開のまま保持するため、アーカイブだけの容量と混同しない。

`native/checkpoint-*.json`は途中受理点の実ファイル。
`resumed-*`は同じ固定コピー実行版・別PIDによる後続全状態と全Run。
`restarts.json`の16地点、後続184状態、Run全体の比較が全てtrueであることを確認する。
除外は`elapsed_seconds`だけで、入力/実行版hash・監査・試行・符号付き0を除外しない。
同一実行版の`plain/run.json`も一致することが必要。

既存ゲート：再組立相対1e-10、非線形1e-8、真の線形1e-10、
土変位/応力相対1e-7、構造並進/回転絶対1e-9、板の絶対1e-8/相対1e-7、
付着拘束1e-10 m。反力合計・荷重仕事・全場/細分化の未定義な受入値は追加しない。
`comparison.json`のゲート未定義nullや、ほぼ0の場の大きな相対差を削除しない。
空間系列は共通物理点・同じ受理progressで照合し、荷重刻み系列とは分ける。
弾性経路の丸め差の変化から漸近次数を認定しない。

`scope.py`は39件のCLI拒否と1件の低強度MC数値未完了を区別する。
後者の`complete=false`を専用拒否の成功へ加算しない。
`mutations.py`は再hash後の17内容改変、checksum未更新1件、型拒否4件を分ける。
22件とも生成元実行版で拒否されることが必要で、異なる実行版hashで止める代替は不可。
これらの拒否や小型ゲート成功をB-02/P6/NC-01の完了にしない。

集計：

```powershell
New-Item -ItemType Directory runs/NEW-b02-summary | Out-Null
.venv\Scripts\python "$B/collect.py" --formal runs/NEW-b02-formal --scope runs/NEW-b02-scope --mutations runs/NEW-b02-mutations --convergence runs/NEW-b02-convergence.json --work-audit runs/NEW-b02-work.json --output runs/NEW-b02-summary
```

`measurements.json.passed=false`、`b02_complete=false`、全場の包括受入nullを保持する。

## 大きいケース：PCだけの本体診断（クラウドでは未実行）

この節は同じ既存入力形式で6144 Tet10を準備する。本体単独の診断であり、
既存独立直接参照の6500土自由度上限を超えるため、独立全場受入の証拠にはならない。
既存の構造容量・幾何ガードに掛かったら拒否のまま保存し、容量を無断で広げない。

```powershell
.venv\Scripts\python "$B/prepare_pc.py" runs/NEW-b02-large --level 3 --steps 16 --offset 0.05
$L='runs/NEW-b02-large'
& target/release/soil-fem-3d.exe solve-stages "$L/mesh.json" "$L/ground.json" "$L/stages.json" "$L/run.json"
& target/release/soil-fem-3d.exe checkpoint-stages "$L/mesh.json" "$L/ground.json" "$L/stages.json" 40 "$L/checkpoint-40.json"
& target/release/soil-fem-3d.exe resume-stages "$L/mesh.json" "$L/ground.json" "$L/stages.json" "$L/checkpoint-40.json" "$L/resumed.json"
.venv\Scripts\python -c "import sys;sys.path.insert(0,'$B');from study import exact;from cases import read;assert exact(read('$L/run.json'),read('$L/resumed.json'))"
```

16〜32 GB RAM、空き20 GB以上、各求解1〜8時間を未実測の目安とする。
`scope.json`に節点・要素数と参照上限の判定、Runにcomplete/failure_reasonが残る。
40は初期平衡16＋設置16＋載荷途中8の受理数で、物理秒ではない。
打切り・OS中断の結果は成功へ変換しない。独立参照をこの大きさで動かす未検証オプションは
追加していない。中規模までの直接参照との一致を大規模へ外挿しない。

## 表示・実CAD・GUI/GPU（PC手順だけ）

小型またはPCの完成Runから既存VTU/CSVを出す：

```powershell
New-Item -ItemType Directory runs/NEW-b02-view | Out-Null
.venv\Scripts\python -c "import sys,json;from pathlib import Path;sys.path.insert(0,'tools');import plate_io;plate_io.export(json.loads(Path('runs/NEW-b02-formal/beam-positive/native/run.json').read_text()),Path('runs/NEW-b02-view'))"
```

`plates.pvd`、`plate-samples.pvd`、各工程VTU、CSV、`plate-fields.json`を確認する。
triangle6点と変位は付着面、neutral_reference_position/neutral_total_displacementは中立面、
生値の位置は中立面。厚さ.1 m、符号付きe=±.05/0、局所右手基底、撤去後空セルを再確認する。
外部ViewerでPVDを開く目視は10〜30分程度。Native GUI/GPUの受入とは区別する。
Native準備器の新方式拒否を削除しない。

実CADは既存`tools/prepare_plate_wall.py`でオフセットなしの別ジョブを作り、元CAD・
メッシュ・法線・境界/材料/荷重のsource SHAを保存する。例：

```powershell
$SolidCraftCli='D:\Apps\SolidCraft\solidcraft.exe'
.venv\Scripts\python tools/prepare_plate_wall.py runs/NEW-b02-CAD-baseline --cli $SolidCraftCli --mesh-sizes-m 4 3 --example plate-wall
```

実行ファイルの配置は利用者のPCに合わせる。準備5〜60分、計算は形状/格子依存。
現行のコンパクト/Native経路からオフセット板への自動変換は未対応である。
上のコマンドを新方式の成功に数えない。実CADでの直接6節点入力、同一付着面の所有・
平面/直線・e*n一致、剛体リンク、全場/仕事・別PID再開・GUI/GPUの工程追跡が
揃うまで実CAD/Nativeの受入は未達とする。今回未対応のJOINT・水・塑性土等を足して拒否を回避しない。

Linuxでは`.venv/bin/python`、`target/release/soil-fem-3d`/`b02-observer`、
`export NAME=value`に読み替える。新先端Windows CIは実際の当該commitの結果だけを
別途確認し、旧成功やLinux結果を転用しない。


## 生証拠の再読込（新規求解ではない）

```powershell
.venv\Scripts\python "$B/package_evidence.py" --verify "$B/evidence-manifest.json" --scratch work/NEW-b02-evidence-verify
```

全part、結合ZIP、全memberのサイズ/SHAが一致することが必要。出力は
`verification.json`と検証用の結合ZIP。数分・空きZIPの2倍程度を見込む。
この確認や`source-snapshot.zip`の監査は新規本体/独立求解に数えない。

`work-audit.json`は新規求解済み両履歴からの補助観測。
瞬間の`F・d`と、全受理点端点の台形積分による増分/累積仕事を分ける。
土RHSから二重計上となる構造付着荷重を引き、構造の力/偶力を明示的に加える。
新設節点の前値は各履歴自身の設置参照とする。施工保持荷重も含む。
設置/撤去の離散端点診断であり、連続物理経路の厳密仕事や新しいエネルギー合格基準とはしない。

Linuxの768 Tet10系列は、変位/応力/板の各場のゲートを満たす一方、
再組立相対残差が本体3.871041030426841e-10・独立1.0508229745493171e-10で既存1e-10を超えた。
この未達を保持する。PCでは元ゲートのまま実測を別出力に残し、Linux値を書き換えない。


封印したLinux実行版/checkpointをWindowsの別実行版へ直接resumeしない。
PCでは同じソースから新しく求解してPC自身のcheckpointを発行し、その固定実行版で復旧する。
