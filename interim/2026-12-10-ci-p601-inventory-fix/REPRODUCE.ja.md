# 再現・台帳追補手順

リポジトリのルートで実行する。Python 3.13、固定依存、Rust toolchain、
元の Git 履歴が必要。checkout では `core.autocrlf=false`、`fetch-depth=0` を用い、
既存の疎取得の除外を保持する。証拠の改行を変換しない。

## Linux の検証

```bash
source /workspace/soil-environment/activate.sh
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
python -m pip install -r requirements-reference.lock.txt
cargo build --locked --release
python interim/2026-12-06-p6-01-status-reconciliation/test_status_audit.py -v
python interim/2026-12-07-p6-03-cloud-reconciliation/test_audit.py -v
python tools/run_ci_python_tests.py --group python-core-0 --output work/p601-python-core-0
python tools/run_ci_python_tests.py --group python-core-1 --output work/p601-python-core-1
python tools/run_ci_python_tests.py --group python-mean-evidence --output work/p601-python-mean-evidence
python tools/run_ci_python_tests.py --group python-spatial-evidence --output work/p601-python-spatial-evidence
```

全件は CI と同じ `tests` と `tools` の unittest discovery であり、4 群に全 ID を重複なく割り当てる。
各群の `python-test-inventory.json` で discovered ID と割当の和集合、実行件数、failure/error/skip、
expected failure、終了コードを確認する。群どうしや該当試験の件数を二重加算しない。
PowerShell を使う試験も実行するため、実物の `pwsh` を PATH に入れる。
未導入による skip は全件成功に数えない。

本作業では全追跡 JSON とソースを配置したうえで、必要なアーカイブ 413 ファイル（初回漏れの cut-volume 1 ファイルを含む）を起点 `5e3d191` の raw URL から取得し、
各ファイルの Git blob SHA-1 を固定 tree と照合してから配置した。試験内の既存 manifest/entry SHA 検査も実行する。
不要な大規模アーカイブは疎配置で省略し、検証したとは主張しない。公式 PowerShell 7.5.2 を作業ディレクトリへ配置し、
`XDG_CACHE_HOME`・`XDG_CONFIG_HOME`・`XDG_DATA_HOME` を書込可能な作業ディレクトリへ設定した。
証拠 JSON、旧 hash、受入基準を変更する準備処理は行っていない。
全件試験前に 3,919 ファイルの指紋を固定し、終了後も全て同一であることを確認した。

## 台帳への追加

```bash
git add -- interim/NEW/validation.json
git commit -m 'Preserve new evidence snapshot before inventory registration'
snapshot_sha=$(git rev-parse HEAD)
python interim/2026-12-10-ci-p601-inventory-fix/register_sources.py \
  --commit "$snapshot_sha" --id N005 interim/NEW/validation.json > work/new-registration.toml
cat work/new-registration.toml >> interim/2026-12-10-ci-p601-inventory-fix/inventory.toml
sha256sum interim/2026-12-10-ci-p601-inventory-fix/inventory.toml
```

ID は現在の台帳にないものを選ぶ。最後の SHA を
`interim/2026-12-06-p6-01-status-reconciliation/audit_status.py` の
`SUPPLEMENTAL_SHA256` へ反映し、既存登録が変更されていないことをレビューする。
追加する JSON ごとに登録する。ヘルパーは不変の 40 桁コミット・その中の実バイトとの一致を必須にし、
boolean/null の JSON Pointer と型付き値も生成する。
追補台帳は TOML であり、自己参照する JSON 証拠を兼ねない。

## Windows CI

修正を `codex/plastic-spatial-convergence` へ push して PR #13 の FEM run を確認する。
HEAD SHA を照合し、8 検証群と `Rust and Python on Windows` 集約の完了・success を確認する。
PC rate-arc gate compatibility も別 workflow として確認する。
キャンセル・時間切れ・実行中・旧先端の成功は、この先端の成功に数えない。

検証済み先端は `08b58fef4727659ad91e2fafaa087aef53e7ee52`。
[FEM run 43330893373](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/43330893373) は完了・success、
8 検証群と集約の全 9 job が success。
[PC compatibility 43330893601](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/43330893601) も
Windows PowerShell / PowerShell 7 の両 job が完了・success。
Windows Python は 321 + 347 + 4 + 4 = 676 件、両固定監査の `... ok` をログで照合した。

結果記録コミット自身で起動する新しい CI は PR の最新 Checks で確認する。
上記成功の SHA と後続の文書コミットを区別し、最終確認は作業完了時の返信にも記録する。

Windows PowerShell 5.1 で台帳を追補する場合は `>` / `Out-File` の既定 UTF-16 保存を避け、
Python の `subprocess.check_output` が返す bytes を `Path.write_bytes` / `open(..., "ab")` で保存する。
台帳の hash は保存された実バイトから計算し、改行や文字コードの差を許容値で吸収しない。

今回の Linux 最終実測は core-0 321、core-1 347、mean 4、spatial 4、合計 676 固有 ID 成功。
初回 core-0 の旧アーカイブ配置漏れ 1 error を保持し、同群 321 件を再実行した。
全群と再実行の method 試行 997 件と、最終 676 固有 ID を分けて確認する。
