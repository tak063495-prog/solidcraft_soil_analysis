# T4 入力検査と T7 分岐走査の再現手順

## 修正前の作業ツリー

リポジトリのルートから、PR #13 の修正前の固定コミットを一時作業ツリーに展開する。

```sh
base=7457a1c96e320d4f3c0d1ac8e2fc0432f2c46137
baseline_dir="$(mktemp -d)/baseline"
git worktree add --no-checkout --detach "$baseline_dir" "$base"
git -C "$baseline_dir" sparse-checkout set --no-cone \
  /Cargo.toml /Cargo.lock /rust-toolchain.toml /src/ /tests/t4_offset_composite.rs \
  '/interim/2026-11-22-impl-t7-physical-instability/*.py'
git -C "$baseline_dir" reset --hard "$base"
cd "$baseline_dir"
```

## #39 修正前の Rust 実測

一時作業ツリーの `tests/t4_offset_composite.rs` の末尾に、既存 fixture を利用する次の診断を追加する。反力差分は同じ実行で 0.01 m と 0.010001 m の荷重フレームから求める。

```rust
fn duplicate_exposed_contact(mesh: &mut Mesh, case: &mut Case) {
    let mut contact = case.plate.contacts[0].clone();
    let source = contact.soil_b.cell;
    let old = mesh.tetrahedra10[source];
    let new = std::array::from_fn(|i| mesh.nodes.len() + i);
    let coordinates = old.map(|node| mesh.nodes[node]);
    mesh.nodes.extend(coordinates);
    let remap = |node| new[old.iter().position(|&n| n == node).unwrap()];
    let mut region = mesh.regions[1].clone();
    region.id = 3;
    region.body = "soil3".into();
    contact.id = "coincident-B".into();
    contact.soil_b.cell = mesh.tetrahedra.len();
    contact.soil_b.nodes = contact.soil_b.nodes.map(remap);
    mesh.tetrahedra.push(mesh.tetrahedra[source].map(remap));
    mesh.tetrahedra10.push(new);
    mesh.material_ids.push(mesh.material_ids[source]);
    mesh.region_ids.push(region.id);
    mesh.regions.push(region);
    case.plate.contacts.push(contact);
}

fn normal_compression(mesh: &Mesh, case: &mut Case, depths: &[f64]) {
    case.steps = vec![Step::Install {
        contacts: case.plate.contacts.iter().map(|contact| {
            (contact.id.clone(), Placement::StressFree)
        }).collect(),
        holder: "displacement supports".into(),
    }];
    for &depth in depths {
        let mut d = vec![0.; mesh.nodes.len() * 3 + 18];
        for contact in &case.plate.contacts {
            for n in mesh.tetrahedra10[contact.soil_b.cell] {
                d[3 * n + 2] = -depth;
            }
        }
        case.steps.push(load(d));
    }
}

fn plate_normal_reaction(case: &Case, frame: &path::Frame) -> f64 {
    case.plate.soil_a.nodes.iter().map(|n| frame.reactions[3 * n + 2]).sum()
}

#[test]
fn reproduce_coincident_exposed_contacts() -> Result<()> {
    let (mut mesh, mut case) = fixture()?;
    duplicate_exposed_contact(&mut mesh, &mut case);
    for (label, contacts) in [
        ("single", vec![case.plate.contacts[0].clone()]),
        ("duplicate", case.plate.contacts.clone()),
        ("reversed", case.plate.contacts.iter().rev().cloned().collect()),
    ] {
        let mut local = case.clone();
        local.plate.contacts = contacts;
        normal_compression(&mesh, &mut local, &[0.01, 0.010001]);
        let run = path::run(&mesh, &local, 3, None)?;
        let force = plate_normal_reaction(&local, &run.frames[2]);
        let perturbed = plate_normal_reaction(&local, &run.frames[3]);
        let slope = (perturbed - force) / 0.000001;
        println!("{label}: force={force:.12} kN; reaction slope={slope:.6} kN/m");
        let count = local.plate.contacts.len() as f64;
        close(force, 5. * count, 1e-10);
        close(slope, 500. * count, 1e-6);
    }
    Ok(())
}
```

```sh
cargo test --locked --release --test t4_offset_composite \
  reproduce_coincident_exposed_contacts -- --nocapture
```

```text
single: force=5.000000000000 kN; reaction slope=500.000000 kN/m
duplicate: force=10.000000000000 kN; reaction slope=1000.000000 kN/m
reversed: force=10.000000000000 kN; reaction slope=1000.000000 kN/m
```

再現用の試験は 1 passed / 0 failed / 0 ignored、既存 6 件は名前指定でフィルタされた。

## #40 修正前の Python 実測

NumPy と SciPy は `requirements-reference.lock.txt` の固定版を使う。受理点はいずれも `x=0`、相対力許容値は `1e-11`。

```sh
PYTHONPATH=interim/2026-11-22-impl-t7-physical-instability python - <<'PY'
import numpy as np
from continuation import Evaluation
from scan_accepted import scan

for sign in (1, -1):
    for loads in ((-.25, 0., .25), (-.25, .25), (0., .25)):
        def op(parent, z):
            x, load = z
            return Evaluation(
                np.array([x * (x*x + sign*load)]),
                np.array([[3*x*x + sign*load]]), np.array([sign*x]), {}, 1.,
            )
        points = [dict(z=[0., load], state={}) for load in loads]
        candidates = scan(op, points, 1e-11, lambda state: state)["branch_candidates"]
        print(sign, loads, [
            (c["accepted_indices"], c["kind"], c.get("reason"),
             [[float(z[-1]) for z in branch] for branch in c.get("branches", [])])
            for c in candidates
        ])
PY
```

| 残差 | 受理荷重 | 候補の順序と結果 |
| --- | --- | --- |
| `x(x²+λ)` | `[-0.25, 0, 0.25]` | 低荷重側 `D`、高荷重側 `D` |
| `x(x²+λ)` | `[-0.25, 0.25]` | `bifurcation`。両枝で負の荷重を 3 点ずつ取得 |
| `x(x²+λ)` | `[0, 0.25]` | `D` |
| `x(x²−λ)` | `[-0.25, 0, 0.25]` | 低荷重側 `D`、高荷重側 `bifurcation` |
| `x(x²−λ)` | `[-0.25, 0.25]` | `bifurcation`。両枝で正の荷重を 3 点ずつ取得 |
| `x(x²−λ)` | `[0, 0.25]` | `bifurcation`。両枝で正の荷重を 3 点ずつ取得 |

表の `D` の理由はすべて `outside the accepted parent interval`。3 点履歴で低荷重側の実在枝を見落とすことを確認した。

## 修正後の試験

修正後のルートで実行する。Python ランナーは `tests` と `tools` の unittest を 7 群へ割り当てる。Linux 専用試験には pytest 9.2.0 を使う。
疎チェックアウトで証拠ファイルが欠けている場合は `git fetch --refetch --no-filter` で Git オブジェクトを取得し、必要なファイルを展開して元の blob と照合する。

```sh
python -m pip install -r requirements-reference.lock.txt pytest==9.2.0
python -m pip check
python -m pytest --version
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo test --locked --release --no-fail-fast
cargo test --locked --release --examples --no-fail-fast
cargo build --locked --release

out="$(mktemp -d)"
for group in python-core-0 python-core-1 python-mean-evidence python-spatial-evidence \
    python-c01-switch python-c01-path python-c01-cause; do
  python tools/run_ci_python_evidence_shards.py --group "$group" --output "$out/$group"
  printf '%s exit=%s\n' "$group" "$?"
done
python -m pytest -q tests/linux
python interim/2026-12-06-p6-01-status-reconciliation/audit_status.py
```

### Linux の実測件数

Rust 1.99.1、Python 3.13.17、pytest 9.2.0。固定版の依存関係を照合し、`pip check` は終了 0。フォーマット、Clippy 全ターゲット、release ビルドはすべて終了 0。

Rust は通常ターゲット、サンプル内の試験、doctest を合わせて **512 passed / 0 failed / 0 ignored**。フィルター除外は 0 件。doctest は 0 件。T4 の 8 件、Python の T7 モジュール 14 件は、それぞれの全件集計に含む。

| Python の群 | 割当 | 実行 | 成功 | 失敗 | エラー | スキップ | 未実行 | 終了コード |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `python-core-0` | 516 | 516 | 515 | 0 | 0 | 1 | 0 | 1 |
| `python-core-1` | 492 | 492 | 491 | 0 | 0 | 1 | 0 | 1 |
| `python-mean-evidence` | 4 | 4 | 4 | 0 | 0 | 0 | 0 | 0 |
| `python-spatial-evidence` | 4 | 4 | 4 | 0 | 0 | 0 | 0 | 0 |
| `python-c01-switch` | 4 | 4 | 4 | 0 | 0 | 0 | 0 | 0 |
| `python-c01-path` | 4 | 4 | 4 | 0 | 0 | 0 | 0 | 0 |
| `python-c01-cause` | 8 | 8 | 8 | 0 | 0 | 0 | 0 | 0 |
| `tests/linux`（pytest） | 7 | 7 | 7 | 0 | 0 | 0 | 0 | 0 |
| **計** | **1,039** | **1,039** | **1,037** | **0** | **0** | **2** | **0** | |

unittest の 1,032 固有 ID を 7 群へ重複・欠落なく割り当てた。2 件のスキップはいずれも Linux 環境の `pwsh` 不在による。CI ランナーはスキップを含む群で終了コード 1 を返す。mean と spatial の既知の `setUpClass` エラーは今回の実行では 0 件だった。

初回の `python-c01-path` は疎チェックアウトの指定漏れで `retained-failures/prototypes.zip` が欠け、3 件成功、失敗 0、エラー 1、スキップ 0、終了 1 だった。元の Git blob `0b38e567edb10b3a19a8997ca862c53885a26814` を展開・照合して群全体を再実行した。上表は再実行の件数を示す。

### 証拠の照合

`git fetch --refetch --no-filter` で起点の到達可能な 15,561 オブジェクトを取得し、欠落が 0 件であることを確認した。過去の斜面失敗 5 ファイルの SHA-256 と `99959691` の Git blob を照合した。資料の末尾は UTF-8 の LF で固定し、試験では実際に CRLF へ変換したコピーも同じ末尾であることを確認した。

P6-01 は `validation.json` を固定コミット `23e8ea6418b683215141236ea9d3a76486aeb3c1`、SHA-256 `c602fd1bc747ee3a1efb640344db759d3ebf348726a043db9063282ecb269437`、Git blob `abf9f40b6782f393183df9a62c35aa25c36afa22` と型付き状態で追補台帳に登録した。監査は 409 ソースで `integrity_ok=true` だった。

A-01 (d) は D（判定不能）、B-02 / B-03 は未達。
FEM の全項目は未達、passed=false。
