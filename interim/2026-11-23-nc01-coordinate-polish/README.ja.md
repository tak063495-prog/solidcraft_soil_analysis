# NC-01・96 Tet10 の倍精度座標仕上げ

明示設定は `settings.dense_residual_method = "coordinate_polish"`。dense LU に限って選択できる。元の NC-01 は未達・分類 D。旧7系の通過と保存新終端の通過は、それぞれ単独の線形系に対する限定受入である。

## 求解と停止条件

初めに通常の部分ピボット LU を使い、元の判定を超えた場合だけ同じ因子で既存の補償残差を1回補正する。それでも超過する場合に、元系の最小二乗目的関数 `‖b−Ax‖₂²` を座標ごとに最小化する。全行列で列順は先頭から末尾までの固定順であり、行列・荷重・要素数に固有の座標は指定しない。

各掃引の先頭で、積の FMA 誤差と2段の TwoSum を使って補正用の `r=b−Ax` を求め直す。列 `aⱼ` ごとの更新幅は `α=(aⱼᵀr)/(aⱼᵀaⱼ)`。列最大絶対値で割って内積を計算し、その倍率を更新幅に戻す。解を `y=fl(xⱼ+α)` とし、残差には実際の変化 `(y−xⱼ)aⱼ` を反映する。探索方向・ステップとも方程式から定まり、解のビットや隣接浮動小数点数を操作する処理はない。

掃引後の合否には nalgebra の元の binary64 行列積を使う。有限解かつ `‖Ax−b‖₂/max(‖b‖₂,1e-30) ≤ 1e-10` が通過条件。全成分が不変、非有限残差、または64掃引で仕上げを打ち切り、超過時は棄却する。追加の計算量は1掃引あたり O(n²)、LU 以外の追加作業領域は O(n)。報告の `iterations` は初回 LU、補償残差の補正、各座標掃引をそれぞれ1と数える。絶対残差は元系の2ノルムで、保存された受理報告では `true_absolute_residual` に記録する。

## 348自由度の保存系

初期値を保存解から借りず、各系を通常 LU から2回独立に解き直した。行列・右辺を含む保存行の SHA-256 を過去の原票と照合した。絶対残差、相対残差、反復数が2回とも一致した。絶対残差の単位は kN。

|保存系|元の絶対残差|元の相対残差|1e-10 判定|補償補正後の座標掃引数|
|---|---:|---:|---|---:|
|旧1|1.051084802e-14|3.144092806e-11|通過|0|
|旧2|6.161279517e-15|3.704163771e-11|通過|0|
|旧3|2.364982820e-15|2.824676385e-11|通過|0|
|旧4|1.480316330e-15|3.567710009e-11|通過|0|
|旧5|8.008133879e-16|3.867753637e-11|通過|0|
|旧6|3.069341408e-16|2.947433169e-11|通過|0|
|旧7|1.597241298e-16|3.046017178e-11|通過|0|
|保存新終端|5.073706993e-16|9.618474879e-11|通過|1|

旧7系は既存の1回補正と同じ値で通過した。保存新終端での既存1回補正は、相対残差 1.0882613663e-10 で棄却されている。同じ系に一般則を適用した座標仕上げが、1掃引で通過した。既存6経路×8保存系の再生でも、以前の原票と判定・元残差・反復数または棄却理由が48件すべて一致した。

## 元入力からの96 Tet10 本計算

同じメッシュ・地盤・工程からソルバ選択だけを加え、初期状態から Linux の release 実行版を1経路実行した。実測版のコミットは `2bf5740fe8e45e06d7dde7e5ab8178a0fbea3b16`（evidence.zip 内の measured-source.bundle）。実行ファイル SHA-256 は `4d997416f02e87ffaf6e8736068f86c2b2d9ffb413444a2ea89bf8bfef5b3a57`。実時間は 916.20 秒、観測ピークメモリは 32.8 MiB、プロセス終了値は 0。

|項目|結果|
|---|---|
|最終受理 λ|0.5647809505462646|
|最終試行 λ|0.5647811889648438|
|最終工程|excavate_upper|
|工程・初期平衡を合わせた受理／棄却増分|188／18|
|完走フラグ|`false`|
|終端の絶対／相対残差|5.5534860799e-16 kN ／ 1.0541749207e-10（棄却）|
|最終受理増分の最後の線形求解の絶対／相対残差|6.0444134732e-18 kN ／ 2.5985141223e-12（通過）|
|受理増分に記録された最大相対残差と、その同じ求解の絶対残差|9.9954859332e-11 ／ 3.0365453893e-14 kN|

現行1回補正の最終受理は 0.5647788047790527、その次に棄却された目標は 0.5647790431976318。新経路でその目標以上となった最初の受理は 0.5647792816162109 だった。受理増分の線形報告 1951 件は全て既存判定内。全棄却で確定状態の保持を確認した。

終了理由は `excavate_upper: cutback exhausted: numerical solve failed: true linear residual 1.0541749207184199e-10 (coordinate_polish, absolute 5.553486079915115e-16, initial 3.608449325051563e-10, steps 21, unscaled min abs U pivot 1.638118908745895e0)`。

並行作業を取り込む前後の2実行を照合した。取り込み前の最終受理 λ は 0.5647809505462646。全工程試行の開始・目標・受理／棄却・理由と、線形報告の反復数・絶対／相対残差は前後で一致した。両方の Run を保存した。

元の NC-01 全撤去は未達・分類 D。保存系の通過や途中の進捗を全撤去の達成として扱わない。

Linux では Rust ライブラリ 78 件、関連 Python 4件、`cargo fmt --all --check`、`cargo clippy --locked --all-targets -- -D warnings` が成功した。

## 保存系とクラウド本計算を再実行する

リポジトリのルートで、次を Linux シェルから実行する。実測時点を厳密に再現する場合は、新しい作業コピーで原票 ZIP から `measured-source.bundle` を取り出し、`git fetch /path/to/measured-source.bundle refs/heads/nc01-coordinate-measured`、`git switch --detach FETCH_HEAD` の順に実行する。前提コミットは `bd0a479cda71a16004b245d19142d5e7d50ebd3a`。保存8系の出力は絶対残差、相対残差、右辺ノルム、反復数を含む。元入力のメッシュ・地盤・工程の SHA-256 は測定 JSON に収録した。

```bash
cargo test --locked --release --lib coordinate_polish_replays -- --nocapture
out="$(mktemp -d)/nc96"
python interim/2026-11-23-impl-nc01-residual-methods/prepare.py "$out" --cloud-original
python - "$out" <<'PY'
import json, sys
from pathlib import Path
folder = Path(sys.argv[1])
case = json.loads((folder / "stages.json").read_bytes())
case["settings"]["dense_residual_method"] = "coordinate_polish"
with (folder / "stages-coordinate_polish.json").open("x", encoding="utf-8") as f:
    json.dump(case, f, indent=2)
    f.write("\n")
PY
cargo run --locked --release --example nc_history_trace -- direct \
  "$out/mesh.json" "$out/ground.json" "$out/stages-coordinate_polish.json" "$out/run"
```

[測定値と SHA-256](measurements.json)・[Run、実行記録、保存系再生と検証のログ](evidence.zip)・[従来の6経路](../2026-11-23-impl-nc01-residual-methods/README.ja.md)
