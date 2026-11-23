# C-01 方法別・組合せ別の再現

PR の同一コミットでリポジトリ固定の Python 依存を導入する。N は荷重進捗の基準分割数。全経路は同じ12 Tet10、32接触 GP。工程上限2,000、基準点への着地、力残差、すべり・散逸の粗解／二分割誤差を全方式で維持する。再試行はアーカイブと別の新しい出力先を選ぶ。

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
F=interim/2026-11-23-c01-method-improvements
W=work/new-c01-methods
cargo build --locked --release --example c01_methods
for side in native independent; do
  for n in 256 512; do
    python "$F/run_case.py" "$W/$side-$n-base" --solver "$side" --steps "$n" --save-seed || exit 1
  done
  python "$F/run_case.py" "$W/$side-1024-base" --solver "$side" --steps 1024 || exit 1
  for n in 512 1024; do
    for mode in w s a ws wa sa wsa; do
      opts=()
      [[ $mode == *w* ]] && opts+=(--seed "$W/$side-$((n/2))-base")
      [[ $mode == *s* ]] && opts+=(--active)
      [[ $mode == *a* ]] && opts+=(--arc)
      python "$F/run_case.py" "$W/$side-$n-$mode" --solver "$side" --steps "$n" "${opts[@]}" || exit 1
    done
  done
done
python "$F/report.py" "$W" "$W/results"
```

これは Bash の手順。Windows で個別実行する場合は同じ Python 引数と `--binary target/release/examples/c01_methods.exe` を使う。PC 実行は行っていない。規模は旧 C-01 と同じで、PC 計画上限は1経路6時間・32 GB、直列。容量・時間の上限に達する前に中止し、途中記録を保存する。768/6144は実行対象外。

`--save-seed` は受理座標の変位を各工程順に保存する。identity.json はメッシュ・材料・ケース・粗経路全記録と seed の SHA-256、実装と N を結び付ける。全7工程を完走した同一実装・同一入力の直前水準のみ使用する。N256 は N512 の初期推定専用で、N512 側の対照・N1024 の粗解取得には event4 の N512 を共用する。基準計算の試行費用は routes.json の coarse_seeds に別記する。Newton の途中の状態や棄却はデータ源に含まれない。

`--active` は全体 Newton の接触接線を局所内分割ごとの開口・固着・すべりの一般化微分に替える。再閉鎖時には閉鎖している部分長の微分も伝播する。残差は従来と同じ接触則で組み立てる。`--arc` は初期地盤平衡後の6工程で平衡解に沿って球面弧長を解く。変位は全節点全成分の RMS を用い、長さ尺度は入力の上面指定変位の最大成分0.02 m、球半径は1/N。各試行は力・経路誤差を再検査し、球面の補正でやり直した荷重点も工程内の1試行に数える。

```sh
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo test --locked --release --example c01_methods
cargo test --locked --release --lib open_sticking_sliding_and_reclosing
python -m unittest discover -s tests -p test_c01_method_choices.py -v
```
