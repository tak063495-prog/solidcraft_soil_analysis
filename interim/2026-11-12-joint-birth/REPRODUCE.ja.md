# ユーザーPCでの再現（案2 接触BIRTH）

同じPRブランチのpush後先端を取得し、`git rev-parse HEAD`、OS、CPU、メモリー、実行版SHAを保存する。PRはマージしない。出力先は毎回新しいフォルダ名にする。過去の不合格・未収束・キャンセル記録を上書きしない。

## 小型受入と全試験

Rust 1.99.1、Python 3.13、固定依存を使用する。全Python試験は過去の証拠ファイルを参照するため、`.github/workflows/fem-ci.yml` と同じ取得範囲を使う。分割された大容量 `interim/2026-11-10-nc-postfold-branches/evidence.tar.gz.part*` は全試験の取得から除外できる。他の必要な過去証拠を省略したまま、試験失敗を仕様の拒否や合格として扱わない。

Windows PowerShell:

```powershell
git fetch origin codex/plastic-spatial-convergence
git switch codex/plastic-spatial-convergence
git pull --ff-only
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements-reference.lock.txt
$env:OPENBLAS_NUM_THREADS='1'; $env:MKL_NUM_THREADS='1'; $env:OMP_NUM_THREADS='1'
$env:CARGO_BUILD_JOBS='2'; $env:RUST_TEST_THREADS='2'
$env:PYTHONUTF8='1'; $env:PYTHONIOENCODING='utf-8'
cargo fmt --all --check
if ($LASTEXITCODE -ne 0) { throw 'fmt failed' }
cargo clippy --locked --all-targets -- -D warnings
if ($LASTEXITCODE -ne 0) { throw 'clippy failed' }
cargo test --locked --release
if ($LASTEXITCODE -ne 0) { throw 'Rust tests failed' }
cargo build --locked --release
if ($LASTEXITCODE -ne 0) { throw 'build failed' }
cargo build --locked --release --example joint_birth_points
if ($LASTEXITCODE -ne 0) { throw 'point driver failed' }
python -m unittest discover -s tests -p 'test_*.py' -v
if ($LASTEXITCODE -ne 0) { throw 'Python tests failed' }
python -m unittest discover -s tools -p 'test_*.py' -v
if ($LASTEXITCODE -ne 0) { throw 'Python tools tests failed' }
python tools/verify_joint_birth.py --solver target/release/soil-fem-3d.exe --point-driver target/release/examples/joint_birth_points.exe --output runs/joint-birth-new
if ($LASTEXITCODE -ne 0) { throw 'BIRTH acceptance failed' }
```

Linuxは `source .venv/bin/activate`、実行版名は拡張子なしに読み替える。固定版のRust/Python依存を使用する。終了コード0に加え、新しい `acceptance.json` の `passed: true` と各範囲を確認する。実測時間を除く比較以外は数値許容値を緩めない。

独立参照は、入力から全体3D弾性問題を組み立てる経路、非線形接触を点積分する経路、保存材料応力から水圧・接触を再組立する経路を別々に記録する。塑性3D解を独立弾性解の合格として扱わない。

## 自分の二ボディへの設定

1. 対応するTri3/Tri6面と独立節点IDを用意し、既存JOINTの両側面・法線・所有要素を確認する。両側ボディが活動すると誕生前は厳密な等変位拘束になる。先端宣言、非適合面、部分再誕生は今回の対象へ混ぜない。
2. JOINT定義へ `"birth_policy":"transfer_force"` または `"stress_free"` を明示する。`joint_birth` 工程の `joints` にIDを列挙する。非活動土の新設は先に乾燥の土 `birth` を行い、その後に接合の `joint_birth` を置く。材料履歴はゼロから作り、旧材料をコピーしない。
3. Aでは既知の圧縮・摩擦円内せん断を誕生前の拘束で伝える。誕生瞬間の変位と全節点伝達力、GP初期力、`U₀`、すべり・散逸0を確認する。引張・摩擦円外・投影不能の拒否を回避しない。GP回収は2000 GP・6000物理節点成分を上限とする。
4. Bでは `transition_audits[].prepared_state` の接触力が厳密に0であり、`birth.jump_m` がその時点の相対変位であることを確認する。工程終了時には外荷重との釣合いで新しい閉合が生じ得る。誕生瞬間と工程終了時を混同しない。
5. 水は乾燥の次に、静水圧・CONTINUOUS・SOIL・PRESCRIBEDを個別に検証する。体積間隙水圧、接合面、露出面を各1回だけ組み立て、開離後にも水圧荷重が残るか確認する。新設土の水理状態を省略した湿潤BIRTHとRichards湿潤JOINTは使わない。
6. ボディDEATHでは除去済み履歴が凍結することを確認し、乾燥再BIRTHとJointBirthで世代が増え、旧すべり・散逸・接触基準が戻らないことを確認する。接合だけの `joint_death` はボディを消さず、接触世代を凍結して活動面を誕生前の等変位へ戻す。

## 別CLIでの途中再開

```powershell
# 上の独立ドライバが作成したケースを利用する例
$dir='runs/joint-birth-new/wet-rebirth-open'
# n は初期平衡を含む受理増分数。attempts と geostatic.attempts から数える。
target/release/soil-fem-3d.exe checkpoint-stages "$dir/mesh.json" "$dir/ground.json" "$dir/case.json" n runs/birth-checkpoint-new.json
target/release/soil-fem-3d.exe resume-stages "$dir/mesh.json" "$dir/ground.json" "$dir/case.json" runs/birth-checkpoint-new.json runs/birth-resumed-new.json
```

`n` を実数値へ置き換える。誕生直前・工程途中・すべり後・撤去後・再誕生途中で行う。ドライバはB経路の11地点を自動検査する。Aは `dry-inherited-shear` に切り替え、初期平衡の受理数を g、contact の受理数を b として n=g、g+1、g+b+1 の3地点を追加検査する。A途中の `last_committed.applied_nodal_load_kn` は初期平衡フレームと完全一致し、`continuation.holding_load` は全成分0であることを確認する。直接計算と再開で同じ実行版を固定する。BIRTH付きチェックポイントは保存地点まで入力から再計算して検証するため、その分の時間とメモリーを確保する。封印の外側ハッシュだけを書き換えた方針・世代・基準・初期力・水圧面節点力の改変は拒否される。

## 大規模・実CAD・GUI/GPU

今回の小型受入を通した後、ユーザーPCで実CADの適合面を作る。体積・面積、節点対応、法線、支持条件、材料・初期応力、荷重・水頭を固定して粗密比較する。Aの回収上限を超える場合は小型証拠の合格を流用せず、疎な節点拘束反力回収・投影の実装と独立検証を別途行う。単に上限や許容値を緩めない。

要素・GP・自由度数、ピークメモリー、計算時間、非線形/線形残差、拘束誤差、全節点力、初期エネルギー、塑性・数値散逸、棄却刻みを保存する。粗密比較や大規模の未収束を受入の完了としない。既存の低プレロードNC全撤去・湿潤塑性すべり全履歴収束の未達は保持する。

GUI/GPUでは `tools/export_joint_birth.py mesh.json run.json 新出力フォルダ` のVTU/PVDと全GP生値を使い、誕生前の節点対、誕生時の基準、世代、開閉・すべり、水圧面荷重、撤去された旧世代を表示する。`joints.pvd` で `joint_generation`、`joint_prebirth_tied`、`has_birth_reference`、誕生基準/初期力/初期エネルギーと水圧を選択し、`joint-generations-NNN.json` の凍結履歴も確認する。表示値とJSON/CSVの節点・GP値を照合する。GPU表示上の一致を数値受入の代わりにしない。実CAD・GUI/GPU・大規模計算の今回の実測合格は宣言しない。
