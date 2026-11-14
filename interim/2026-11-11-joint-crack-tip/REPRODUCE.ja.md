# ユーザーPCでの再現

同じPRブランチの今回のpush先端を取得し、作業前に `git rev-parse HEAD` を保存する。PRはマージしない。過去の `interim` 記録を上書きしない。大容量の分割証拠を必要としない全試験の取得範囲は `.github/workflows/fem-ci.yml` のsparse設定と同じにできる。

## 小型の自動受入（Windows PowerShell）

Rust 1.99.1、Python 3.13を使用し、固定依存をインストールする。出力先は毎回新しい名前にする。既存出力への上書きは行わない。

```powershell
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements-reference.lock.txt
$env:OPENBLAS_NUM_THREADS='1'; $env:MKL_NUM_THREADS='1'; $env:OMP_NUM_THREADS='1'
$env:CARGO_BUILD_JOBS='2'; $env:RUST_TEST_THREADS='2'
$env:PYTHONUTF8='1'; $env:PYTHONIOENCODING='utf-8'
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo test --locked --release
cargo build --locked --release
python -m unittest discover -s tests -p 'test_*.py' -v
python -m unittest discover -s tools -p 'test_*.py' -v
python tools/verify_joint_crack_tip.py --solver target/release/soil-fem-3d.exe --output runs/crack-tip-new
foreach ($mode in @('continuous','soil','prescribed')) {
  python tools/verify_joint_crack_tip_water.py --solver target/release/soil-fem-3d.exe --output "runs/crack-tip-water-new-$mode" --mode $mode
  if ($LASTEXITCODE -ne 0) { throw "wet crack-tip failed: $mode" }
}
python tools/verify_joint_crack_tip_refinement.py --solver target/release/soil-fem-3d.exe --output runs/crack-tip-local-h-new
```

各コマンドの終了コード0と `measurements.json` の `passed: true` の双方を確認する。時間切れ・キャンセル・不完全Run・未収束は合格にしない。Linuxでは仮想環境の有効化を `source .venv/bin/activate`、実行版を `target/release/soil-fem-3d` に読み替える。最終Linux版の実行SHA・全試験件数は同梱JSONと照合する。

チェックポイントCLIは `checkpoint-stages mesh.json ground.json stages.json 受理増分数 checkpoint.json`、再開は `resume-stages mesh.json ground.json stages.json checkpoint.json resumed.json`。受理増分数には初期平衡も含む。自動ドライバが各工程途中で別プロセスを起動し、実測時間以外の全Run一致を検査する。先端宣言・面集合・所有側・使用ソルバー/入力の改変を封印の付け替えで通さない。

## 実CADの先端形状

1. SolidCraftで亀裂を体積分割面として準備し、亀裂の終端線を**実在のTet10辺列**へ一致させて共形メッシュを作る。要素内部で終端する形状は入力しない。分離面と未分離接着面を明確に分ける。同一ボディ内部でも両側所有要素を列挙する。
2. `faces` に分離前の頂点面と所有要素対を、`tip_edges` に両端・中間節点を列挙する。`split_joint_mesh.py --prepared 元準備フォルダ --spec split.json --output 新準備フォルダ` を実行する。
3. 出力の `mesh.provenance.json` 内 `joint_split.crack_tip_declarations` をJOINT定義へ設定する。面集合全体・所有側を省略しない。座標・材料・領域体積を元メッシュと比較し、未宣言の接着面は同じIDであることを確認する。
4. 先に乾燥で、次に完全水没の静水圧/定常Darcyで計算する。濡れた先端の一般塑性・すべり全履歴・部分水没の追加受入は別途必要。RichardsやBIRTHへ変更して既存拒否を回避しない。
5. `tools/export_stages.py mesh.json result.json 新出力フォルダ` でVTU/PVD/CSVを書き出し、既存 `tools/verify_stage_package.py` のCLI（`--help`）に従って入力・実行SHAと全場を照合する。表示では先端3節点の共有、内部ジャンプ、未分離接着、片側撤去後の露出面と接触力消失を確認する。

## 大規模計算・GUI/GPU

実CAD先端、大規模の局所h系列、GUI/GPU表示はユーザーPCで実施する。疎ソルバーと十分なメモリーを用意し、要素数・GP数・自由度数、実行SHA、ピークメモリー、経過時間、各増分の真の線形残差/非線形残差、棄却理由を保存する。許容値はREADMEの値を維持する。

局所h細分では幾何・材料・支持・荷重を固定し、同じ観測点の変位、同一定義のエネルギー、先端から一定距離以上離れた**同一領域**の応力を比較する。特異先端応力の点値一致を求めない。開口ヒンジが生じる拘束では、接触支持があるものとして解析を継続しない。途中保存から別CLIで再開し、実測時間だけを除いて直接経路と比較する。

GUI/GPUでは亀裂面の両側IDと全GP生値を保持し、変形倍率・法線・開閉/すべり状態・撤去済み面を確認する。表示上の一致を数値受入の代わりにしない。別PC/GPUでの成功を今回のLinux証拠へ追加する場合も、実行版と入力・出力・ログを別フォルダに残す。低プレロードNC全撤去と湿潤すべり全履歴の既存未達は、今回の限定合格で更新しない。

WindowsでビルドしたexeのSHAはLinux版と異なる。同じPCの直接計算と再開では同じexeを使い、Linuxのチェックポイントを別exeへ封印し直して移すことはしない。実行版を変える場合は保存した入力から新規計算として再現する。

撤去保持力と解放率0の独立監査は、上記ドライバの出力に対して次で再現できる。湿潤では `--run` を各出力内の `continuous` / `soil` / `prescribed` フォルダへ変更する。

```powershell
python interim/2026-11-11-joint-crack-tip/audit_holding.py --run runs/crack-tip-new/construction --output runs/crack-tip-new/holding.json
```
