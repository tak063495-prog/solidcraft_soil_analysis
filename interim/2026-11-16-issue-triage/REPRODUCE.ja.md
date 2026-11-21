# 再現手順とユーザーPCで残す作業

保存先ラベルは2026-11-16、実施日は2026-11-11（日本時間）。起点は`d8616cefb526e70f422ec8c164afd2353c8f8149`、PR #13・`codex/plastic-spatial-convergence`。修正コミットは`issues.json`の`repair_commit`を参照する。PRはマージしない。入力・旧archive・過去の失敗記録・hashを更新しない。各コマンドの終了コードを確認し、失敗・時間切れ・キャンセルを後続の成功と混ぜない。

## 取得とLinux全件チェック

既存の実studyから分離した検証checkoutで、同ブランチを取得する。

```bash
git fetch origin codex/plastic-spatial-convergence
git switch codex/plastic-spatial-convergence
git pull --ff-only origin codex/plastic-spatial-convergence
git rev-parse HEAD
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
python3.13 -m venv .venv
.venv/bin/python -m pip install -r requirements-reference.lock.txt
export PATH="$PWD/.venv/bin:$PATH"
export CARGO_BUILD_JOBS=2 RUST_TEST_THREADS=2
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
export PYTHONUTF8=1 PYTHONIOENCODING=utf-8
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release
cargo test --locked --release
```

Pythonは既存の無欠落分割runnerで全4群を実行する。出力先`work/triage-recheck`は新規にする。runnerはtestsとtoolsの全IDを発見し、skip/expected failure/unexpected success・欠落・重複を成功にしない。

```bash
for group in python-core-0 python-core-1 python-mean-evidence python-spatial-evidence; do
  .venv/bin/python tools/run_ci_python_tests.py --group "$group" --output "work/triage-recheck/$group" || exit
done
```

4つの`python-test-inventory.json`で`discovered_test_ids`が等しく、実行IDの和集合が発見IDと一致し、重複なし・全件実行・失敗/エラー/skip等0を確認する。今回の結果と照合するファイルは`validation.json`と`python-coverage.json`。集中回帰や再実行を全件数へ足さない。

今回Rustは、取得待ちの間にGit blobから隔離コピーしたbuild checkoutで実行し、144 Rust/config入力が起点・最終ブランチとも一致することを確認した。Python全4群は同ブランチの通常checkout上で実行する。隔離コピーと実行場所は`rust-build-source-binding.json`、各コマンド/終了コード/時間は`validation/*-commands.json`に残す。初回の隔離コピーでは過去の斜面不合格JSONが不足し、lib 59成功/1失敗で停止した。元Git bytesを復元して全件再実行し、初回ログは`retained-failures/initial-rust-*`へ保持する。試験や不合格値を書き換えず、初回59件を最終成功数へ足さない。

## #5/#14の最小反例と正常対照

新規回帰5メソッドを単独実行する。テストは通常実行の内部からnormal/`-O`/`PYTHONOPTIMIZE=1`の別Pythonを起動する。試験自体をskipしていない。

```bash
.venv/bin/python -m unittest discover -s tests -p 'test_issue_triage_followups.py' -v
.venv/bin/python interim/2026-11-16-issue-triage/probe_baseline_followups.py --output work/NEW-triage-baseline-comparison.json
```

probeは`git show d8616ce:...`から元checkerの実bytesを読み、現行と同じ解析式fixture・依存を使って比較する。元のGitファイル、旧archive、実studyは変更しない。#5は元スクリプトから本人確認blockだけをAST抽出し、偽のWindows APIだけでnormal/最適化の違いを観測する。CIM・実保持handle・signal・実process終了・成功receipt生成は行わない。

#14は土の実外面12面（head2/no_flow10）を使う。末尾no_flowの`gradient_outflow_m3_s=NaN`、`relative_nodal_balance_error=-Infinity`、同`-1`は起点でAPI受理・全3モードexit0/passed=true、修正後はAPI拒否・全3モードexit1/report未生成となる。正常reportの修正前後一致も出力する。既存GPのNaN/±Infinity・形状不正は`test_p2_soil_source.py`に残る。

修正前の入口テストは、最適化ガードがないため失敗する。`baseline-followup-tests.log/json`はその期待した失敗記録で、最終成功数には含めない。正常モードのhelp、snapshot親子判定と#1の所有子停止回帰、他の既存P2は次で実行できる。

```bash
.venv/bin/python -m unittest discover -s tests -p 'test_p2_*.py' -v
cargo test --locked --release --test epoch_srm_gate
```

## #15：小型契約試験と本物の旧証拠再計測

```bash
.venv/bin/python -m unittest discover -s tests -p 'test_nc_replay_contract.py' -v
.venv/bin/python -m unittest discover -s tests -p 'test_historical_utf8_isolation.py' -v
.venv/bin/python -m unittest discover -s tests -p 'test_spatial_replay_isolation.py' -v
.venv/bin/python -m unittest discover -s tests -p 'test_nc_mean_stress_evidence.py' -v
.venv/bin/python -m unittest discover -s tests -p 'test_nc_spatial_evidence.py' -v
```

小型契約試験2メソッド/16条件は合成archiveの隔離と破損拒否であり、本物の数値結果の検証を代用しない。本物はmean 328 part・5,490,818,201 bytes・1,314 entries、spatial 55 part・920,331,265 bytes・411 entriesを読み、各4試験を実行する。新しいMCC求解や大規模計算は起動しないが、旧証拠の読込・復元・再計測は行う。

既存`remeasure.py`の手動実行も可能だが、出力先は存在しない新規ディレクトリにする。旧source/hashに現在のbytesを上書きして合格にしない。試験はWindowsに限り既存の派生float差記録を使い、Linuxはstrictな計測JSONバイト一致を要求する。今回この条件・例外を拡張していない。両者とも`numerical_acceptance=false`、`current_checkout_numerical_certification=false`を維持する。

```bash
.venv/bin/python interim/2026-11-10-nc-mean-stress-release/remeasure.py --output work/NEW-mean-replay
.venv/bin/python interim/2026-11-10-nc-spatial-release/remeasure.py --output work/NEW-spatial-replay
```

元の封印からの不変性は次で照合できる。差分0を期待する対象はpart/manifest/指紋/旧計測であり、後続修正を含む現行sourceと旧sourceの差そのものを0にする手順ではない。

```bash
git diff --exit-code f19ce642647207a87b9dda6cc31b006b85b5936c HEAD -- \
  'interim/2026-11-10-nc-mean-stress-release/evidence.tar.gz.part*' \
  interim/2026-11-10-nc-mean-stress-release/archive-sha256.json \
  interim/2026-11-10-nc-mean-stress-release/measurements.json \
  'interim/2026-11-10-nc-spatial-release/evidence.tar.gz.part*' \
  interim/2026-11-10-nc-spatial-release/archive-sha256.json \
  interim/2026-11-10-nc-spatial-release/source-fingerprints.json \
  interim/2026-11-10-nc-spatial-release/measurements.json
```

元実行バイナリは両archiveに同梱されない。元の実行版SHA/実行記録が不変であることと、旧バイナリを再ビルドして再実行することを混同しない。現行sourceは元epoch sourceと異なるままでよい。archive内source/依存は元SHAへ厳密に束縛し、現行sourceを旧実行に付け替えない。

## Windows CIとユーザーPCでの確認

Windowsの通常Pythonと固定依存で同じCargoチェックを実行し、`tools/run_ci_python_tests.py`の全4群を実行する。PowerShellでは各コマンド直後の`$LASTEXITCODE`を確認する。`test_p2_resume_ownership`のWindowsはCRT main-thread signalと実所有子を使う経路であり、外部POSIX signalやCIM停止全体の受入ではない。最適化safe-pauseは常に拒否されるべきで、実studyへ`--execute`して試さない。

新先端Windows CIでは、head SHAとmerge tree、全job、各NC4試験名の`... ok`、`Ran 4 tests`、最終集計とskipなしを確認する。非選択matrix stepのskipと試験のskipを区別する。今回の旧d8616ce run 42951479845の成功を新先端へ転用せず、タイムアウト/キャンセルは成功にしない。

以下は**今回クラウドで実行しないユーザーPC作業**。必要な元入力・固定参照・仕様を揃え、既存の以下の手順と`todo.ja.md`の受入不足を対応させる。入力がない実CADや亀裂形状・基準を創作しない。

- (a) `../2026-11-13-joint-ext-combos/REPRODUCE.ja.md`：湿潤先端BIRTH/複数JOINT/構造、定常流BIRTH、Richards湿潤JOINTは拒否継続。既存の小型限定合格で代用しない。
- (b) `../2026-11-14-plate6-curved-offset/REPRODUCE.ja.md`：曲面幾何/director/曲率/MITC・パッチ・古典殻h系列、オフセットの未検証複合、実CAD/GUI/GPUは未達。現行曲面拒否を外して計算しない。
- (c) `../2026-11-15-wet-joint-convergence/REPRODUCE.ja.md`：256増分以上・768要素以上等の系列、共通物理点の全履歴、一般塑性、完了までの実checkpoint別プロセス再開を検証。親面平均や7端点だけで全履歴を合格にしない。
- `../../docs/LOW_PRELOAD_NC.ja.md`、`../../docs/EXACT_FS_REAL_SECTION_SRM.ja.md`、`../../docs/GRIFFITHS_2D_COMPARISON.ja.md`、残件表：NC全撤去、斜面h系列5%、正確なFs/実断面SRM、大規模/長時間Richards、実CAD・Native UI・別PC/GPU・完成配布は未達。却下済みNC調整は行わない。
- 実Windows安全停止を別途確認する場合は隔離環境で自分で起動した使い捨て親子だけを使い、snapshot本人情報・PID作成時刻・保持handleのimage/作成時刻、確認不能時の拒否を観測する。無関係processや実studyを終了しない。

## issue操作

`github/issue-snapshot.json`は操作前の全件取得記録。`issues.json`は起点/修正後判定と推奨操作、`comments/`は全issueの日本語文案。push後、権限がある場合は#5/#14へ追記の修正根拠、#15へ解決根拠を投稿し#15だけcloseする。既存closedは再closeやreopenを要しない。拒否された場合は文案のみ残す。PR #13のマージ操作は含めない。
