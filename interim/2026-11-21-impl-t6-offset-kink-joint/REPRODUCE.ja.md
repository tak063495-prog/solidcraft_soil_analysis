# T6 オフセット折れ継手の再現手順

PR #13 の履歴を取得した作業ツリーで実行する。科学計算の起点は `d390fb7327f11cdbf38739be46187068a0997601`。パスはリポジトリのルートを基準とする。以下は bash の例。Windows では Python 3.13、同じ Rust ツールチェーン、拡張子 `.exe` の実行版を使用する。出力先には未使用の名前を指定する。

## 依存関係と過去診断の先行照合

```sh
git log -4 --oneline --decorate origin/codex/plastic-spatial-convergence
python -m pip install -r requirements-reference.lock.txt pytest==9.2.0
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
export CARGO_BUILD_JOBS=2 RUST_TEST_THREADS=2
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
export PYTHONUTF8=1 PYTHONIOENCODING=utf-8
T=interim/2026-11-21-impl-t6-offset-kink-joint
B=interim/2026-12-02-tri6-offset-kink-b03
P=interim/2026-12-01-tri6-offset-composite-b02
O=runs/t6-reproduction
mkdir -p "$O"
```

試験を起動する前に失敗診断5件と各群で使う分割ファイルを作業ツリーに展開する。完全チェックアウトを使える場合も下記の照合を行う。CI の各専用群は結合済みアーカイブや数値結果の代用品を使わず、復元した原本から一時領域に展開する。疎チェックアウトの場合、CI 全群の和集合から外れる `nc-postfold-branches/evidence.tar.gz.part*` は不要。PowerShell を起動する Python 試験には実物の `pwsh` が必要で、Linux の検証版は 7.5.4。

```sh
python - <<'PY'
from pathlib import Path
import hashlib, json, subprocess
def git(*args):
    return subprocess.check_output(['git',*args])
manifest=Path('interim/2026-11-26-wet-soil-birth-darcy-a03/validation/restored-original-slope-failures.json')
reports=json.loads(manifest.read_bytes())
assert len(reports)==5
for row in reports:
    p=Path(row['path'])
    assert hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256'],p
    assert git('hash-object',str(p))==git('rev-parse','99959691:'+p.as_posix()),p
for entry in git('ls-tree','-rz','HEAD').split(b'\0'):
    if not entry: continue
    meta,raw=entry.split(b'\t',1)
    p=Path(raw.decode())
    if p.as_posix().startswith('interim/2026-11-10-nc-postfold-branches/evidence.tar.gz.part'):continue
    if not any(ext in p.name for ext in ('.zip','.tar','.gz','.npz')):continue
    assert git('hash-object',str(p)).strip()==meta.split()[2],p
PY
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release --bin soil-fem-3d
cli=target/release/soil-fem-3d
python -m unittest discover -s tests -p test_t6_offset_kink.py -v
python "$T/t6_study.py" "$cli" "$O/structural"
python "$T/t6_convergence.py" "$O/structural"
```

## 入力と出力

`offset-kink-joint CASE NEW_RUN [--until INDEX | --resume RUN]` は T5 の `shell` と明示した `joint` を含む JSON を読む。生成入力は `structural/n2-distort0/case.json`。各翼に別々の構造節点 ID を与え、同じ土節点に競合する付着条件を設定すると診断する。有向直線上の `line_m` と各翼の `edge_nodes` は端点、中間点、次の端点の順に並べる。一つの接続に剛な翼群を列挙できる。ヒンジ・線ばねの対象は明示した2翼。翼の `section` 省略時は T5 の `shell.section` を使用する。切り詰め後の各翼本体を Tri6 で入力し、翼対ごとの `separators` で材料領域の分離を宣言する。判定は板厚を含む。局部の接合材料は入力に含めていない。

長さ・力・応力の単位は m、kN、kPa。回転は全体 XYZ の微小回転 rad。支持と節点荷重は各付着点の物理自由度順で与える。未拘束一覧から継手条件を消去する際、独立な反力を回収できない支持・重複条件は節点と階数を含む診断となる。線荷重はイベント、荷重 ID、作用する翼、有向弧長の区間、力／長さ、偶力／長さを指定する。線ばねの剛性は kN·m/(rad·m)。力と偶力は符号付きの面別集約値と、整合線質量で射影した密度を記録する。断面量の軸は配列に順に並ぶ3基底ベクトル。圧力・自重は T5 の面／体積求積を使う。

受理フレームごとに `active.joint.stations` と `active.joint.cuts` を読む。設置・撤去の保持イベントで一意な継手反力を回収していない欄は null。膜・曲げ・せん断の原場、物理／drilling エネルギーは `active.elements`。撤去済み世代は `archived`。再開時は同じ入力・同じ実行版で受理列を再演算する。異なる OS からの再開結果のバイト一致は未確認。

```sh
"$cli" offset-kink-joint "$O/structural/n2-distort0/case.json" "$O/cut.json" --until 1
"$cli" offset-kink-joint "$O/structural/n2-distort0/case.json" "$O/end.json" --resume "$O/cut.json"
cmp "$O/end.json" "$O/structural/n2-distort0/run.json"
```

## B-03 の旧ケースと全件試験

```sh
python "$P/build_observer.py" "$O/observer-build"
python "$B/probe.py" "$O/b03-probe" --cli "$cli"
python "$B/supplemental.py" "$O/b03-supplemental" --cli "$cli"
python "$B/geometry_audit.py" "$O/b03-probe" "$O/b03-geometry.json"
python "$B/orientation_audit.py" "$O/b03-orientation.json"
python "$B/study_b03.py" "$O/b03-formal" --binary target/release/b02-observer
python "$B/convergence_b03.py" "$O/b03-formal" "$O/b03-convergence.json"
python "$B/rigid_modes.py" "$O/b03-formal" "$O/b03-rigid.json"
python "$B/checkpoint_mutations.py" "$O/b03-formal" "$O/b03-mutations"
python "$B/verify_contract.py" "$O/b03-formal" "$O/b03-probe/scope.json" "$O/b03-supplemental/supplemental.json" "$O/b03-mutations/mutations.json" "$O/b03-contract.json"
python "$B/collect_b03.py" "$O/b03-formal" "$O/b03-probe/scope.json" "$O/b03-supplemental/supplemental.json" "$O/b03-mutations/mutations.json" "$O/b03-summary"
cargo test --locked --release --all-targets --no-fail-fast
cargo test --locked --release --doc
for group in python-core-0 python-core-1 python-mean-evidence python-spatial-evidence python-c01-switch python-c01-path python-c01-cause; do
  python tools/run_ci_python_evidence_shards.py --group "$group" --output "$O/$group"
  echo "$group $?" >> "$O/python-exit-codes.txt"
done
python -m pytest -q tests/linux
python interim/2026-12-06-p6-01-status-reconciliation/audit_status.py
```

各群の `python-test-inventory.json` の割当・実行・失敗・エラー・未実行を照合する。Linux の NC 準備エラーは元の厳密検査の結果を記録する。三水準の自己収束は同一の観測量、等比細分、同方向・減少の隣接差がそろうときだけ算出する。768／6144 要素メッシュは規模上限のため確認不可・未達。

B-03 は未達、passed=false。
