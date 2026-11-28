# 再現手順（ユーザーPC）

リポジトリルートで実行する。Rust 1.99.1、Python 3.13と固定依存を使う。各コマンドの終了コードを確認し、非ゼロ・未収束・途中停止を合格にしない。出力先は新しい名前にする。

```powershell
git fetch origin codex/plastic-spatial-convergence
git switch codex/plastic-spatial-convergence
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
python -m venv .venv
.venv\Scripts\python -m pip install -r requirements-reference.lock.txt
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release
cargo test --locked --release
.venv\Scripts\python -m unittest discover -s tests -p 'test_*.py' -v
.venv\Scripts\python -m unittest discover -s tools -p 'test_*.py' -v
$env:OPENBLAS_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
.venv\Scripts\python tools/record_plate6_curved_offset.py runs/plate6-offset-user-pc
```

Linuxでは `.venv/bin/python`、環境変数は `export NAME=1` に読み替える。生成する `acceptance.json` は、曲面未達のため全体 `passed=false`、平面オフセットの全ゲート成功時だけ `thickness_offset.passed=true` となる。全体falseをエラーなくtrueへ変更しない。途中で例外が出た場合、完成した受入JSONがないことを成功にしない。

`examples/plate6-offset-{mesh,ground,stages}.json` は検証用の同一入力。直接部材で方式 `mindlin_tri6_offset`、`attachment_offset_m` を明示する。付着面の節点座標を中立面へ移してから入力すると別の幾何になるので、元の地盤節点座標を保持する。既存の `match_soil_node` / 一節点集合を使う。省略eは0。コンパクト面宣言とNative準備器への新方式追加は未検証のため拒否のまま。

## 大規模・実CAD・GUI/GPU（今回のクラウドでは未実行）

1. 従来モデルと4.740/3.466 μmを含む元結果を保存し、別ジョブを作る。平面・直線の対象地盤Tri6面と向き付き法線を確認する。曲面の平坦化、曲面を多数の平板へ分けたものを今回の曲面Tri6合格と扱うこと、法線の自動平均で折れ目を隠すことはしない。
2. 平面だけの直接部材入力へ明示方式とeを指定する。有限 `|e|<=t/2`、共有節点で同じe n、乾燥弾性MC、板と梁の線形小変形、設置/載除荷/撤去の範囲を保持する。JOINT・水・土BIRTH/DEATH/MATSET・塑性土等の拒否を回避しない。梁を共有する場合、その節点座標・直接荷重は付着面基準であることを確認する。
3. 既存の構造 `storage_budget` で容量を明示する。大規模疎構造系の性能は今回未達なので、容量拒否を切り捨て・自由度省略で回避しない。実モデルの粗密・厚さ・drilling・境界感度は別系列として測定し、小型直接解からの一般収束を推定しない。
4. 受理点での保存と同一バイナリの別CLI再開、棄却時の確定状態、全節点釣合い、入力からの独立解、エネルギー・力/モーメント、付着面内点/辺内点の並進一致を元の許容値で再検査する。`tools/plate6_scope_rejections.py` の拒否記録は、力学ベンチマークの合格にはならない。
5. `tools/plate_io.py` の出力は付着面triangle6、中立面の位置/並進と生値を明示する。VTUを再読込し、厚さ、eの符号、局所軸、設置基準と撤去時の空フレームを確認する。現在のNative構造フィールド抽出は新方式を拒否するため、公式Viewerが表示できたとは報告しない。Nativeの独立監査と表示側にオフセット意味論を実装・検証してから、実GUI/GPUの全工程、断面・生値・失効/Undo・別GPUを検査する。
6. 曲面を進める場合はREADMEの未実装定義を実装し、6剛体モード/余剰モード、曲面適合パッチ、平面退化のround-off一致、二次地盤面の内点一致を先に検証する。Scordelis-Lo屋根・ピンチ円筒・半球殻のうち可能な参照問題の原典・単位・支持・荷重・測定点を固定し、二次曲面幾何を維持したh系列の参照差を保存する。差が減らなければ不合格として残す。曲率上限、全域Jacobianの正性と法線連続性の受理基準も、定式化の検証前に緩めない。現行本体は曲面を拒否するので、ここは今後の開発手順である。

低プレロードNC全撤去や湿潤JOINT塑性すべり全履歴収束を、この手順の成功で完了へ変更しない。Windows CIは対象pushの実際の結果を確認して別途記録する。

## 指定先端との従来機能回帰を再計算する場合

元ソースは別worktree/別targetでビルドし、現行ソースや実行中のバイナリを置換しない。原点はmainではなく指定先端cecbb39である。

```powershell
git worktree add --detach runs/plate-baseline-src cecbb39edf7359f02ca3968998532ebd8612d544
cargo build --locked --release --manifest-path runs/plate-baseline-src/Cargo.toml --target-dir runs/plate-baseline-target
.venv\Scripts\python interim/2026-11-14-plate6-curved-offset/legacy_compare.py --baseline runs/plate-baseline-target/release/soil-fem-3d.exe --output runs/plate-regression-user-pc
```

Linuxでは実行ファイル名の `.exe` を除く。従来CST/DKT・平面Mindlin・乾湿JOINT・先端JOINT・範囲拡張(a)の11完成Runについて、計測時間と実行バイナリSHAのみを除き厳密比較する。これは全過去問題の再解析を意味しない。元の70ビルドソースのGit blobとSHAは `baseline-source-proof.json` に保持する。

符号付き0を含めた厳密JSON監査も実施する。通常の数値 `==` で+0と−0を同一視せず、計測時間を除く正規化JSONの一致を確認する。

```powershell
.venv\Scripts\python interim/2026-11-14-plate6-curved-offset/audit_exact_histories.py runs/plate6-offset-user-pc runs/plate6-offset-user-pc-exact.json
```

今回の `exact-history-audit.json` は正負eの10再開・棄却なし経路・キャッシュ・e=0について、`legacy-exact-history-audit.json` は従来11ケースについて符号付き0を含む一致を保存する。棄却なし経路で区別するのはREADMEに明記した棄却操作/初回刻みとその入力ハッシュだけであり、物理場や受理履歴は除外しない。

疎チェックアウトの場合、過去の失敗JSONや封印済み証拠ZIP/分割tarも省略しない。全件試験はそれらのハッシュと不合格判定の保持を検査する。必要に応じて `git sparse-checkout disable` で元の追跡ファイルを展開する。不足をテスト除外や参照値変更で回避しない。

## 今回のLinux完走結果（再接続後に確認）

|検査|結果|件数|
|---|---|---:|
|cargo fmt --all --check|成功|書式検査|
|cargo clippy --locked --all-targets -- -D warnings|成功|全ターゲット|
|cargo build --locked --release|成功|release版|
|cargo test --locked --release|成功・失敗/無視/除外0|451件（69集計）|
|Python tests 自動発見・無フィルター全件|成功・失敗/エラー/スキップ0|497件|
|Python tools 自動発見・全件|成功・失敗/エラー/スキップ0|2件|

Python合計499件。testsの実測時間は1689.599秒、toolsは0.017秒。再実行・集中試験を重複加算しない。全件ログは logs/rust-complete-tests.log、logs/python-clean-complete.log、logs/python-tools-clean.log、機械可読の集計と各ログSHAは validation.json。

前回の接続断の間も全件プロセスは継続し、今回その完了ログを確認した。環境の再初期化や差分再適用は不要だった。検証済み828ソースとのSHA一致と現在の実行バイナリSHAを確認し、fmt/clippyも再接続後に追確認した。曲面Tri6は未実装・明示拒否で acceptance.json.passed=false のまま、平面・直線辺Tri6のオフセットだけ限定受入である。Windows CIは未確定。
