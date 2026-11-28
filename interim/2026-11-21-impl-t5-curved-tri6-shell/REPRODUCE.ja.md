# T5 曲面 Tri6 の再現手順

リポジトリ直下で、PR #13 の履歴を含むチェックアウトを使う。Rust は 1.99.1、Python は 3.13。以下の Linux 例では、各出力に未使用のディレクトリ名を指定する。Windows では実行ファイル名に .exe を追加する。Python の試験は過去診断の展開後に起動する。

```sh
git log -4 --oneline --decorate origin/codex/plastic-spatial-convergence
python -m pip install -r requirements-reference.lock.txt pytest==9.2.0
export CARGO_BUILD_JOBS=2 RUST_TEST_THREADS=2
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
export PYTHONUTF8=1 PYTHONIOENCODING=utf-8
T=interim/2026-11-21-impl-t5-curved-tri6-shell
B=interim/2026-12-03-tri6-curved-b01
O=runs/t5-reproduction
mkdir -p "$O"
python - <<'PY'
from pathlib import Path
import re,subprocess
src=Path('src/plane_strain_2d.rs').read_text().split('fn archived_slope_failure_files_are_unchanged()')[1]
paths=re.findall(r'"(interim/[^"]+-FAILED.json)"',src)
assert len(paths)==5
for p in paths:
    sha=subprocess.check_output(['git','rev-parse','99959691:'+p]).strip()
    data=Path(p).read_bytes()
    assert subprocess.check_output(['git','hash-object','--stdin'],input=data).strip()==sha,p
PY
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release
cargo build --locked --release --bins --examples
cargo test --locked --release --test t5_curved_shell -- --nocapture
python -m unittest discover -s tests -p test_t5_curved_shell.py
python "$T/modes.py" target/release/examples/curved_shell_operators "$O/modal"
python "$T/manufactured.py" target/release/soil-fem-3d "$O/mms"
python "$T/standard.py" target/release/soil-fem-3d "$O/standard"
python "$T/section_work.py" "$O/section-work.tsv" "$O/mms" "$O/standard"
python "$T/attached.py" target/release/soil-fem-3d "$O/attached"
```

## 入力、自由度、設置、再開

`curved-shell CASE NEW_RUN` は乾燥・等方線形弾性、一定厚さの曲面シェルを扱う。座標は付着面、単位は m、kN、kPa を用いる。標準殻は各公刊問題の一貫単位系で数値を入力する。Tri6 の順序は頂点3個、辺01・12・20の中間節点。各構造節点は全体並進3成分、全体回転3成分をこの順で持つ。`group` は滑らかな面群を指定する。`directors` の明示値は面群・節点 ID を組にして与える。乾燥土がある場合、全付着面を所有セルと6節点で指定する。土の付着節点の並進は構造と共有し、残りの土並進を土節点 ID 順に末尾へ置く。

`steps` は T4 の install・load・remove を使う。load の nodal_load は付着参照点の力と偶力、plate_traction_kpa は全要素の中央面の死荷重。自重を別途体積積分する。face_loads を指定する場合は、1始まりの load イベント番号・面 ID・中央面からの z により面を選ぶ。pressure_kpa は選んだ外面への内向きを正とする。同じ層への重複指定を検査する。soil.owners がある場合には exposed_side の符号と実付着の宣言を照合する。

別プロセス再開の例を示す。保存点は受理済みイベント番号。設置時にその時点の配置を無応力参照として固定し、撤去世代を凍結する。同じ入力・実行版で受理履歴を再演算し、保存内容を照合して次の受理点から進む。

```sh
cli=target/release/soil-fem-3d
"$cli" curved-shell "$O/attached/case.json" "$O/cut.json" --until 4
"$cli" curved-shell "$O/attached/case.json" "$O/resumed.json" --resume "$O/cut.json"
cmp "$O/attached/full.json" "$O/resumed.json"
```

## B-01 の旧受入ケース、全件試験と台帳

```sh
python "$B/probe.py" "$O/legacy-scope" --cli "$cli"
python "$B/spherical_probes.py" "$O/legacy-sphere" --cli "$cli"
python interim/2026-12-01-tri6-offset-composite-b02/build_observer.py "$O/observer-build"
python "$B/study_b01.py" "$O/legacy-fields" --binary target/release/b02-observer
python "$B/checkpoint_mutations.py" "$O/legacy-fields" "$O/legacy-mutations"
python "$B/rigid_modes.py" "$O/legacy-fields" "$O/legacy-rigid.json"
python "$B/curved_comparison.py" "$O/legacy-fields" "$O/legacy-comparison"
python "$B/work_audit_b01.py" "$O/legacy-fields" "$O/legacy-work.json"
cargo test --locked --release --no-fail-fast
cargo test --locked --release --examples
for group in python-core-0 python-core-1 python-mean-evidence python-spatial-evidence python-c01-switch python-c01-path python-c01-cause; do
  python tools/run_ci_python_evidence_shards.py --group "$group" --output "$O/$group"
done
python -m pytest -q tests/linux
python interim/2026-12-06-p6-01-status-reconciliation/audit_status.py
```

各 Python 群が参照する ZIP と分割アーカイブは .github/workflows/fem-ci.yml の当該群と同じ条件で展開する。過去診断5ファイルの欠落、異なる SHA、あるいはアーカイブ欠落があれば、試験前に該当 Git blob を復元する。PowerShell 依存試験は T4 と同じ pwsh 7.5.4 で実行する。Windows の物理結果は同じ Windows 実行版で再開して照合する。

B-01 は未達、passed=false。
