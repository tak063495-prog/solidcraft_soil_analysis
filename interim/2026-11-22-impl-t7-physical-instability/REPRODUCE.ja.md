# T7 の再現

PR #13 の履歴を持つ作業ツリーのルートで実行する。入力と証拠の元バイトを保つため、取得前に Git の core.autocrlf を false にする。以下は bash 表記。Windows では Python 3.13 と同じ Rust ツールチェーンを使い、実行版の末尾に .exe を付ける。O は未使用のディレクトリにする。計算は経路ごとに逐次実行する。

```sh
git log -4 --oneline --decorate origin/codex/plastic-spatial-convergence
git merge-base --is-ancestor c3e26502ba5e05f94982734bd5b5fd77f080beca HEAD
python -m pip install -r requirements-reference.lock.txt pytest==9.2.0
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
export CARGO_BUILD_JOBS=2 RUST_TEST_THREADS=2 PYTHONUTF8=1 PYTHONIOENCODING=utf-8
T=interim/2026-11-22-impl-t7-physical-instability
O=runs/t7-reproduction
mkdir -p "$O"
```

## 過去診断と台帳

[T6 の再現手順の冒頭](../2026-11-21-impl-t6-offset-kink-joint/REPRODUCE.ja.md)にある試験前の照合を実行する。過去の斜面失敗5ファイルは各 SHA-256 と 99959691 での Git blob の両方が一致すること。CI 全群が参照する ZIP・分割ファイルは HEAD の Git blob に照合する。今回の証拠も含めて展開する。完全取得の代わりに疎取得を使う場合、旧 nc-postfold-branches/evidence.tar.gz.part* だけを除外できる。Linux の PowerShell 依存試験には pwsh 7.5.4 を使う。P6 の固定コミット検査には Git の履歴と対象 blob が必要である。

```sh
python interim/2027-01-04-a01-high-precision-nullspace/a01hp_evidence.py --verify "$T" --extract "$O/reopened"
python "$T/trace_archive.py" --restore "$O/reopened/science/mesh-18"
python "$T/trace_archive.py" --restore "$O/reopened/science/mesh-144"
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release --examples --bins
python -m unittest discover -s tests -p test_t7_physical_instability.py -v
```

## 参照解と平行2面

```sh
cli=target/release/soil-fem-3d
bridge=target/release/examples/t7_boundary
python "$T/references.py" "$O/references"
python "$T/arch_elements.py" "$O/beam"
python "$T/parallel.py" "$O/mesh-18" "$cli" "$bridge" --elements 18
```

references.py は完全柱の解析解、公刊入力の深円弧の射撃法、元入力から独立に解く頂点変位制御の選点法を実行する。arch_elements.py は同じ深円弧の全弧を15・30・60梁要素で解く。15要素では集中荷重を頂点の節点に置くため、固定側を7要素、ヒンジ側を8要素とする。

各メッシュの run.json は通常の荷重制御の結果。数値停止時には failure_reason と、停止前の last_committed が残る。independent-last.json.gz は独立計算で得た受理履歴の終端と次の段階の定義を含む。BIRTH の前後は別の状態として扱う。独立計算で初期保持力を解放する場合は新世代側で進捗ゼロの平衡を再検査する。

ラベルの順序は計算系・無次元弧長・接触小刻み倍率。最初の64受理増分で観測区間の終端に達しなかった経路は step budget と記録される。最後の受理点から別プロセスで続ける。次の例は18 Tet4 の細かい弧長である。固定された入力、コード、実行版、親履歴の SHA が一致する同一実行環境で再開する。

```sh
dir="$O/mesh-18"
label=native-0.01-1
src="$dir/$label.json.gz"
digest=$(python -c 'import hashlib,sys;print(hashlib.file_digest(open(sys.argv[1],"rb"),"sha256").hexdigest())' "$src")
python "$T/resume_interval.py" "$dir" "$label" "$bridge" "$src" "$digest" "$dir/$label-segment-1.json.gz"
python "$T/parallel.py" "$O/mesh-144" "$cli" "$bridge" --elements 144
```

144 Tet4 は同じ面対、支持、元の18四面体を保つ8分割。144の場合は弧長0.02・接触倍率1を固定する。summary.json にある step budget の各経路について、前コマンドと同様に新しいセグメントを生成する。各セグメントの標準出力の stop が step budget なら、そのセグメントを新しい親として、末尾の番号を1つ増やす。数値棄却で終わった経路はその理由と最後の受理点が結果になる。封印から開いたファイルの実行版 SHA が手元の実行版と異なる場合は元入力から再計算する。Windows/Linux の実行時間や全浮動小数点表記のバイト一致は受入基準に含めていない。

```sh
python "$T/compare.py" "$O" "$bridge"
python "$T/scan_accepted.py" "$O" "$bridge" "$O/accepted-spectra.json"
```

compare.py は同じ open 進捗で、同じ活動枝の単調区間内に実際の平衡を求める。分岐、開閉や世代をまたぐ補間は行わない。支持・反力・応力・接触すべりの単位、残差の分母、臨界区間の補正行列と右辺は対応する経路の JSON に含まれる。再開・刻み依存の確認では parent、result_checkpoint、generation と進行方向を照合する。規模上限を超える768／6144要素の手順は設けていない。

## 全件

```sh
cargo test --locked --release --all-targets --no-fail-fast
cargo test --locked --release --doc
all_groups=0
for group in python-core-0 python-core-1 python-mean-evidence python-spatial-evidence python-c01-switch python-c01-path python-c01-cause; do
  rc=0
  python tools/run_ci_python_evidence_shards.py --group "$group" --output "$O/$group" || rc=$?
  echo "$group $rc" >> "$O/python-exit-codes.txt"
  if [ "$rc" -ne 0 ]; then all_groups=1; fi
done
python -m pytest -q tests/linux
python interim/2026-12-06-p6-01-status-reconciliation/audit_status.py --root .
test "$all_groups" -eq 0
```

各群の python-test-inventory.json で発見件数・割当・実行・成功・失敗・エラーを確認する。P6-01 は出典 JSON の固定コミット、SHA-256、Git blob、状態キー、未登録資料件数を確認する。Windows CI は push した PR 先端に紐付けて確認する。

A-01 (d) は未達、結論 D。
