# Windows CIの45分打切りと検証範囲を維持した分割

起点は指定されたPR #13先端 `42473341dcfaa3c6c969a28f66fe0f5e0d312037`。
元のActions run [42615618270](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/42615618270)、job `129154626423` の実ログを取得した。計測抜粋は `ci-timeout-observed.json`。
チェックアウトされたPR merge SHAは `f775ab0d21c4778086f9eeaa5452e7184a705ddb`（42473341をmainへ仮統合したGitHubの試験用コミット）。PRをマージしたという意味ではない。

- release全ターゲットビルド：22分11秒。
- Rust試験：12:43:37〜12:59:23 UTC、約15分47秒。全Rust試験は成功。
- Python試験開始：12:59:24 UTC。ジョブ開始12:16:29から約42分55秒経過。
- 13:01:28 UTCにジョブ全体の45分上限でキャンセル。Pythonのアサーション失敗は記録されていない。後続の湿潤構造4種・NC診断は未実行。

これらは元ログの時刻である。日本時間ではそれぞれ+9時間。
原因はビルドと全検証の直列累積時間で、キャンセル時に表示されていたPython試験だけの問題とは認定しない。

## 対処

45分の上限を維持し、Windows 2025上の5つの必須matrix群に分けた。
Rust lib/bin/doctestと例のreleaseビルド、Rust integration 2群、全Python tests/tools、既存の湿潤構造4種と有界NC計測の5群である。
各群の失敗で他群をキャンセルせず、旧チェック名 `Rust and Python on Windows` の集約ジョブは**全群成功時だけ成功**する。

`tools/run_ci_rust_tests.py` はCargo metadataの実integrationターゲットを名前順の交互割当で全列挙する。固定された部分リスト、名前除外、ignore追加、許容値変更、release最適化変更はない。将来追加したターゲットも自動的に入る。全群にソース指紋・試験一覧・終了コード・ログを保存する。
Pythonの全discover 2本、既存の湿潤4種、NCの20/60秒native診断と120秒独立参照をそのまま維持した。
重い計算をCI既定から外してはいない。interimの長い擬似弧長追跡は以前からCIの自動探索外であり、本追補のユーザーPC手順で別に再現する。

CI成功は5 kPa NC全撤去の合格を意味しない。時間上限・未収束・点数上限を破壊や合格に変換しない。
新構成のWindows実測結果は別の記録に保存する。Linux試験結果や旧先端CIを新構成の成功として扱わない。

## ユーザーPCでの同一検証

Rust toolchain/lockfileとPython 3.13、`requirements-reference.lock.txt` を使う。BLASは1スレッド。
リポジトリルートのPowerShellで、各行の終了コードが0であることを確認する。

```powershell
$env:CARGO_BUILD_JOBS='2'; $env:RUST_TEST_THREADS='2'
$env:OPENBLAS_NUM_THREADS='1'; $env:MKL_NUM_THREADS='1'; $env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'; $env:PYTHONIOENCODING='utf-8'
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
python tools/run_ci_rust_tests.py --group rust-core --output work/ci-core
python tools/run_ci_rust_tests.py --group rust-integration-0 --output work/ci-integration-0
python tools/run_ci_rust_tests.py --group rust-integration-1 --output work/ci-integration-1
cargo build --locked --release
python -m unittest discover -s tests -p 'test_*.py' -v
python -m unittest discover -s tools -p 'test_*.py' -v
```

追加の湿潤4種・有界NCの同一コマンド、同一予算、ログ出力先は `.github/workflows/fem-ci.yml` の `verification` 群を参照。出力先は毎回未存在のフォルダにする。再実行のために封印済みjobを上書きしない。
