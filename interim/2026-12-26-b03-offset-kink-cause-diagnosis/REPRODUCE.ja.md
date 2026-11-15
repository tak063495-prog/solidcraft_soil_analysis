# B-03 原因診断の再現と小型PC手順

GPD WIN Max 2、使用可能メモリ約48GBを想定。**PCでは実行していない**。
以下の時間/メモリはPC未測定の目安で、Linux実測や合格条件とは別である。
旧証拠を上書きせず、新しいwork出力先を使う。終了コード・中断・未完了も保存する。
新たな768/6144要素、40GB超・1経路1日超の計算は計画しない。
上限を越える確認は**規模上限のため確認不可・未達**。許容値・支持・拘束は元のまま。

## 1. 環境と保存証拠の再解析（新規求解なし）

リポジトリ直下のPowerShell。Python3.13、requirements-reference.lock.txt、Rust1.99.1を使う。
Gitの旧封印を保持するためautocrlf=falseで取得する。新しい文書末尾試験は実CRLFにも対応するが、
数値JSON/アーカイブのSHAは生バイトで検証する。

```powershell
$B='interim/2026-12-26-b03-offset-kink-cause-diagnosis'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:CARGO_BUILD_JOBS='2'
$env:RUST_TEST_THREADS='2'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
python "$B/b03cause_evidence.py" --verify $B --extract work/b03cause-pc-unsealed
python "$B/b03cause_report.py" --patches work/b03cause-pc-unsealed/b03cause-final-patches --controls work/b03cause-pc-unsealed/b03cause-final-controls --output work/b03cause-pc-remeasured.json
python -m unittest discover -s tests -p test_b03_offset_kink_cause_diagnosis.py -v
```

想定1–3分、1GB未満。出力はwork/b03cause-pc-*。
見る値：全part/entry SHA、23＋23の保存全体解、力とモーメントを分けた残差、F*d−2U、
3面/異厚の宣言順、旧h系列の次数null、passed=false。
別OSで再集計の最終ビットが異なる場合は差と演算環境を別記し、許容値で吸収しない。

## 2. 同じ小型内部診断を新規実行

```powershell
cargo build --locked --release
python "$B/b03cause_build.py" --output work/b03cause-pc-build
python "$B/b03cause_driver.py" --probe work/b03cause-pc-build/b03cause_probe.exe --output work/b03cause-pc-patches
python "$B/b03cause_controls.py" --probe work/b03cause-pc-build/b03cause_probe.exe --solver target/release/soil-fem-3d.exe --output work/b03cause-pc-controls
python "$B/b03cause_report.py" --patches work/b03cause-pc-patches --controls work/b03cause-pc-controls --output work/b03cause-pc-observations.json
```

buildは10–30分/2–6GB、診断は1–5分/2GB未満の目安。
自由構造パッチは2 Tri6・最大12節点72DOF、23本体要素組立＋23入力のみ独立組立。
別PID再実行と頂点全順の単板評価は全体求解と別集計する。
Model::newだけの40要求は機械経路の完走ではない。
標準CLIは既存受理6＋拒否12、補足スクリプトは旧方式受理3＋拒否7。
コンパクト/Native準備はheadlessの拒否確認であり、表示は実行しない。

見る値：build-binding.jsonの無改変ソースSHA、production_scope_expanded=false、
剛体6/別ID12モードと全固有値、生の全自由度残差、独立剛性/解との差、
局所の腕偶力を落とした会計差、全回転/順序の結果。
同入力の再実行はresponse_bytes_exactを確認する。丸めで変えた入力を同一入力と呼ばない。
生産Modelの拒否を外さず、未定義継手をこの内部全体解で受理しない。

## 3. 既存許可域だけの12/96 Tet10・4/8/16刻み

旧B-03の生成器/独立参照/公開API観測器をそのまま使う。
これは既存の限定モデルの再現であり、追加の非ゼロ腕連続継手ベンチマークではない。

```powershell
$Old='interim/2026-12-02-tri6-offset-kink-b03'
$P='interim/2026-12-01-tri6-offset-composite-b02'
python "$P/build_observer.py" work/b03cause-pc-history-build
python "$Old/study_b03.py" work/b03cause-pc-limited-series --binary target/release/b02-observer.exe --kinds kink-zero --level-max 1
python "$Old/convergence_b03.py" work/b03cause-pc-limited-series work/b03cause-pc-limited-convergence.json
python "$Old/rigid_modes.py" work/b03cause-pc-limited-series work/b03cause-pc-limited-rigid.json
```

平面対照1＋選択幾何4、各経路は本体/独立/観測なし/実checkpoint復旧。
45/225土節点、12/96 Tet10、荷重4/8/16刻みのみ。
想定合計10–60分、4–16GB、ディスク5GB以上。出力work/b03cause-pc-limited-*。
同じコマンドで`--kinds normal-sign-pair`、`step-split`、`kink-split`、`height-step`を
**それぞれ新出力先に一つずつ**指定できる。別IDの結果を連続継手へ一般化しない。

見る値は元の全GP/節点/板点、支持反力、瞬間F*d、離散累積仕事、設置参照、別PID復旧。
元の再組立1e−10、非線形1e−8、真の線形1e−10、土変位/応力相対1e−7、
構造並進/回転絶対1e−9、板生値絶対1e−8/相対1e−7、付着拘束1e−10 mを保持。
反力合計/仕事/一般折れ目の新しい許容値は作らない。
同じ20受理座標・物理点・向き付き断面量で比較し、異なる格子の節点力を直接引かない。
この2水準だけから空間次数を出さない。4/8/16の仕事差と丸め領域の状態差を分ける。

## 4. 追加h系列・折板/交差殻ベンチマーク

**非ゼロ腕の連続継手は実装・仕様決定待ち。実行可能な生産コマンドは存在しない。**
上のcontrolsが現行拒否を確認するコマンドである。
MPC、director、回転解放、spine、荷重作用面、設置/保存を選んでいない状態で、
文献値に一致する支持・荷重・drilling係数を探す手順は書かない。

旧B-03の同じ土/板h生成器で2水準の次を作ると規模上限に該当する。
**元生成器の次水準による確認は規模上限のため確認不可・未達**とし、新しいPC実行コマンドを示さない。
上限内の第3格子を使うなら、元の物理領域/付着/支持/荷重・全体基底のN/M/Q・
共通物理点と反力の転置制限を保持する非二進格子生成器/比較器を先に具体化する必要がある。
現時点はその実装・比較仕様決定待ちで、架空のオプションや収束率を補わない。

将来の継手仕様が決まった後も、最初の候補は2/8/32 Tri6以下、全構造DOF1000以下、
固定の小さい土12/96 Tet10以内での剛体6モード/余剰モード・仕事共役荷重・独立全体解とする。
目標予算は1経路30分/8GB以下だが、未実装なので実測時間/メモリ・実コマンドは未確定。
出力先計画名はwork/b03cause-pc-junction-{mesh}。**実装・仕様決定待ち**である。
8要素以下で実測して上限内と見込める場合にだけ次のコマンドを具体化する。
文献のねじり解放/soft支持/交差部3Dを、本件の既存問題へ黙って代入しない。
上限内で参照モデル/元基準の充足を証明できなければ未達。外挿/Richardson補正で合格にしない。

## 5. 外部表示・実CAD・Native/GPU（PCだけ、未実行）

既存の限定直接入力で作ったRunに対する小型の外部表示準備だけは実行可能。

```powershell
New-Item -ItemType Directory work/b03cause-pc-view | Out-Null
python -c "import sys,json;from pathlib import Path;sys.path.insert(0,'tools');import plate_io;plate_io.export(json.loads(Path('work/b03cause-pc-limited-series/kink-zero-h0-s4/native/run.json').read_text()),Path('work/b03cause-pc-view'))"
```

想定準備1分/1GB未満、外部Viewer確認10–30分/2GB以内。
見る値：plates.pvd/plate-samples.pvd、VTU/CSV、付着面と中立面、符号付きe、
局所軸、生値位置、設置基準、撤去後の空セル。別IDをmerge/weld/法線平均で統合しない。
この外部表示をNative表示/GPUの受入へ数えない。

コンパクト/Nativeのoffset準備はcontrolsの補足コマンドで拒否される。
その先のNative表示は**実装・仕様決定待ち**。
実CADも元CAD・法線/付着/中立面・接続・同じ支持/荷重を確定し、
12/96 Tet10以内の入力へ保持できることが前提。現状、対象CADとその変換仕様が未確定なので
準備/求解コマンド・メモリ/時間の確定は**実装・仕様決定待ち**。
GPUを選択する実装済みソルバーオプションはない。架空の起動フラグを示さない。
GUI/GPU・実CAD・Native表示をクラウドで実行した証拠はない。

## 6. 全件試験とP6台帳監査

```powershell
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release --all-targets
foreach ($g in @('rust-core','rust-integration-0','rust-integration-1')) {
  python tools/run_ci_rust_tests.py --group $g --output "work/b03cause-pc-validation/$g"
}
cargo test --locked --release --examples
foreach ($g in @('python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence','python-c01-switch','python-c01-path','python-c01-cause')) {
  python tools/run_ci_python_evidence_shards.py --group $g --output "work/b03cause-pc-validation/$g"
}
python interim/2026-12-06-p6-01-status-reconciliation/audit_status.py --root . > work/b03cause-pc-p601-audit.json
```

想定1–4時間/2–16GB。既存のソフトウェア全件で、1経路1日超の新しい物理計算ではない。
必要な旧アーカイブと一時展開用のディスク空きを確保し、欠落をskip/除外で回避しない。
見る値：全発見/実行ID、成功/失敗/error/未実行/skip、SOURCE_NOT_INVENTORIED=0、
integrity_ok、全体未達の型付きfalse/null。旧失敗記録と今回の実行を分ける。
push後Windows CIは未確定。Linuxの成功を新先端Windowsの成功へ転用しない。
