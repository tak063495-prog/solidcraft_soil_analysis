# A-01・A-02・A-04 の小型比較を再実行する

リポジトリのルートで実行する。クラウドで下記の4経路を実行済み。ユーザー PC での実行は任意である。出力には新しいディレクトリ名を使う。各ドライバは入力、受理履歴、棄却理由、観測量をその下に保存する。

## 準備

CI と同じ Python 3.13、Rust（`rust-toolchain.toml`）、依存固定ファイルを使う。以下は PowerShell 用。Windows では履歴文書のバイト検査のため、チェックアウト時に Git の自動 CRLF 変換を無効にする。Linux はパス区切りをそのまま使い、実行ファイルの `.exe` を外す。

```powershell
python -m pip install -r requirements-reference.lock.txt
cargo build --release --locked --bin soil-fem-3d --example joint_a02_projection_probe
$t = 'interim/2026-11-22-impl-missing-features-a01-a02-a04'
$w = 'runs/missing-features-small'
$cli = 'target/release/soil-fem-3d.exe'
$probe = 'target/release/examples/joint_a02_projection_probe.exe'
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
```

## 小型経路

```powershell
python "$t/a01.py" "$w/a01" $cli
python "$t/a01.py" "$w/a01-common" $cli --open-fraction 0.02
python "$t/a02.py" "$w/a02" $cli $probe
python "$t/a04.py" "$w/a04" $cli
```

|出力先|規模と読み方|クラウド経過時間|子プロセス最大 RSS|
|---|---|---:|---:|
|a01|18・144 Tet4、元入力。BIRTH、shear、open の最終受理点・停止理由を読む|2.82 秒|71.4 MiB|
|a01-common|同じ格子で open の目標を元幅の 0.02 とした別入力。世代・向き・枝が一致する点を読む|3.28 秒|72.3 MiB|
|a02|24 Tet10。元の S1-A、定応力軸圧縮、節点力投影と拒否対照を読む|2.95 秒|107.6 MiB|
|a04|12 Tet4。既定と貯留ありを各2時刻進め、独立全連成解、面の体積差・流量、水収支を読む|0.48 秒|107.6 MiB|

時間は release ビルド後のドライバ全体の実測。RSS は Linux の子プロセス最大値で、プロセス全体の同時合計は未測定。PC ではビルド後、各行につきメモリ 2 GB、10 分を目安とする（推定）。順次実行する。この手順の最大格子は144要素。A-01 の元入力の open 棄却、A-02 の不適合入力の拒否も結果に含む。閾値や増分試行数を変えて再試行しない。

各経路の `results.json` と同梱の `observations.zip` 内の同名ファイルで比較する。アーカイブは生成入力、ネイティブ履歴、独立解との差と SHA-256 表を含む。新規実行時は、ソルバのタイミング・浮動小数点結果に環境差があり得る。配布時のバイトは次で照合できる。

```powershell
python -c "import hashlib,json,pathlib,sys,zipfile; d=pathlib.Path(sys.argv[1]); p=d/'observations.zip'; assert hashlib.sha256(p.read_bytes()).hexdigest()==json.loads((d/'acceptance.json').read_bytes())['observations_zip_sha256']; z=zipfile.ZipFile(p); m=json.loads(z.read('sha256.json')); assert set(m)==set(z.namelist())-{'sha256.json'}; assert all(hashlib.sha256(z.read(k)).hexdigest()==v for k,v in m.items()); print(len(m),'files verified')" $t
```

## 全試験と台帳

証跡を含むチェックアウトで実行する。リポジトリの CI が定める7群は、`tests` と `tools` の unittest をモジュール単位で重複なく分割する。追加の Linux 専用 pytest は Linux 環境で実行する。下記の全試験は上表の小型ドライバの資源見積もりとは別である。

```powershell
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --release --locked --bins --examples
cargo test --release --locked --all-targets --no-fail-fast
cargo test --release --locked --doc
$groups = 'python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence','python-c01-switch','python-c01-path','python-c01-cause'
foreach ($g in $groups) {
    python tools/run_ci_python_evidence_shards.py --group $g --output "$w/tests/$g"
}
python interim/2026-12-06-p6-01-status-reconciliation/audit_status.py
```

Linux 専用の追加群は `python -m pip install pytest==9.2.0` の後、`python -m pytest -q tests/linux` で実行する。P6-01 の監査はソルバを起動せず、台帳に固定したコミット・Git blob・SHA-256 と保存状態値を照合する。
