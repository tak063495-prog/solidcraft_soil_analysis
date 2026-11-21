# T4 の再現手順

PR #13 の履歴を含むチェックアウトを使う。実行場所はリポジトリ直下。T4 は Rust 1.99.1、Python 3.13 に固定した。Windows ではチェックアウト前に `git config --global core.autocrlf false` を設定する。既存の診断と台帳は過去の Git blob を参照するため、単一コミットの ZIP 展開では検証できない。

## 乾燥パッチと既存ケース

以下は PowerShell 7 の例。出力先には未使用のディレクトリ名を指定する。

```powershell
$env:CARGO_BUILD_JOBS = '2'
$env:RUST_TEST_THREADS = '2'
$env:OMP_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$env:OPENBLAS_NUM_THREADS = '1'
$env:PYTHONUTF8 = '1'
$env:PYTHONIOENCODING = 'utf-8'
$T = 'interim/2026-11-21-impl-t4-offset-plate-composite'
$O = 'runs/t4-reproduction'
python -m pip install -r requirements-verification-win-cp313.lock.txt
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release --bins --examples
cargo test --locked --release --test t4_offset_composite
python "$T/study.py" "$O/patch"
python interim/2026-12-01-tri6-offset-composite-b02/build_observer.py "$O/observer-build"
python "$T/existing_cases.py" B-02 "$O/b02"
python "$T/existing_cases.py" B-03 "$O/b03"
python "$T/collect.py" "$O" "$O/summary"
```

`study.py` は既存12 Tet10 メッシュの上下層を分け、上層を 0.1 m 移して適合する一対の面を作る。入力は `mesh.json` と `case.json` に保存する。自由度順は、全土節点の世界座標並進 (m)、次に付着面の6節点それぞれの世界座標回転 (rad)。付着側の並進は土 A の節点と同じ自由度を使う。露出側の節点順は板の頂点3・辺中点01、12、20に対応させる。法線の向きは板節点順、実接触面の選択は `exposed_side` で指定する。入力の `plate_traction_kpa` は中央面の追加死荷重、板自重は断面の単位体積重量と厚さから別に1回加える。`nodal_load` は付着参照点系の外力と偶力で、板面荷重を二重に入力しない。

`install` と `remove` では名前を与えた保持機構がその瞬間の配置を固定する。拘束反力はその配置で再計算する。次の `load` で自由にする成分、規定配置と荷重を指定する。各荷重行が一つの受理増分である。未収束の増分は結果を出さない。`transfer_force` は供給者名、全16 GP の有効 traction、付着節点の並進・回転36成分と相手節点の並進18成分を要求する。面上の一定圧縮・せん断を既知の力として移す場合も、各 GP の非引張・摩擦条件と全54成分の照合を通す。

新しい CLI を直接使う場合の例を示す。

```powershell
$cli = 'target/release/soil-fem-3d.exe'
& $cli offset-composite "$O/patch/mesh.json" "$O/patch/case.json" "$O/frame9.json" --until 9
& $cli offset-composite "$O/patch/mesh.json" "$O/patch/case.json" "$O/resume9.json" --resume "$O/frame9.json"
```

再開には作成時と同じ実行版を使う。実行版・入力の SHA-256 を確認し、確定済みイベント列を再計算して保存内容全体を照合した後に次の増分から継続する。同じ受理点の再加算や新しい世代の作成は行わない。独立系は既存の多項式 Tri6 と GL125 土剛性、GL25 面積分を使う。自由な A 側土並進と板付着並進を含む18個の tie 成分を乗数で組み立てる。独立解を先に求め、本体の受理変位、行列、反力を初期値に使わない。GP 比較では今回の固着・開口経路での二次 traction 場を本体の求積位置で評価する。すべり境界を横断する一般場への適用は未確認。

`existing_cases.py` は B-02・B-03 の元の入力生成器を呼び、3入力の SHA-256 を各旧 measurements.json と照合する。本体の観測あり・なし、独立解、別プロセス再開を経路ごとに直列実行する。各出力先に全比較量の comparison.json、入力、ログを保存する。経路の照合終了後に一時全場ストリームを解放する。旧拒否入力も元の probe / scope で再実行し、診断、終了コード、結果ファイルの有無を記録する。既存の12要素入力だけを再求解し、96・768・6144要素への細分は行わない。T4 新規経路の上限は24 Tet10、1経路8 GB・2時間。

## 全件試験と台帳

```powershell
cargo test --locked --release
cargo test --locked --release --examples
foreach ($g in 'python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence','python-c01-switch','python-c01-path','python-c01-cause') {
    python tools/run_ci_python_evidence_shards.py --group $g --output "$O/tests/$g"
    Write-Output "$g exit=$LASTEXITCODE"
}
python interim/2026-12-06-p6-01-status-reconciliation/audit_status.py
python -m unittest discover -s interim/2026-12-06-p6-01-status-reconciliation -p test_status_audit.py
```

Linux は `python -m pip install -r requirements-reference.lock.txt pytest==9.2.0` を使い、上の7群に加えて `python -m pytest -q tests/linux` を実行する。PowerShell に依存する既存試験には pwsh 7.5.4 を使う。書き込み可能な領域へ XDG_CONFIG_HOME、XDG_CACHE_HOME、XDG_DATA_HOME を設定する。pytest、Rust examples、P6-01 専用57件は、それぞれの実行範囲を区別して集計する。Python の各群は出力先の python-test-inventory.json に全割当・実行数・失敗・エラー・スキップ・終了コードを記録する。

固定証拠は inventory.toml の T4 登録にある40桁コミットの measurements.json。register_sources.py が元バイト、SHA-256、Git blob ID と JSON の状態値を検査し、audit_status.py が台帳全体の指紋を検査する。再現時に新しい実行ファイルや時刻に合わせて既存の登録を更新する必要はない。

新規再現では経路ごとに3入力の SHA-256、全場の比較量、再開後の確定履歴、拒否診断を照合する。`evidence.zip` にはプロセス ID・所要時間・実行ファイルの SHA-256 も入るため、別環境で作った ZIP のハッシュは変わり得る。登録済みアーカイブ自体の検証には `measurements.json` の `evidence_sha256` を用いる。
