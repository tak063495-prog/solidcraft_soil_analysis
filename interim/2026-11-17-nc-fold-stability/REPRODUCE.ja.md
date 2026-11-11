# 再監査・Windows PC での追加計算

保存先の日付は依頼ラベル。実施日は2026-11-11（日本時間）。すべてリポジトリのルートから実行する。出力先は毎回未存在の名前を使う。PR #13をマージしない。指定起点 `5f564621c7396c2a605edbfe5451cd800cf318d3` の子孫であることを確認する。

以下の長い前史・刻み系列は **Windows PC 専用**。クラウドで実行したのは旧証拠の読取再監査、解析例、局所材料写像の再生と6点の短経路だけである。終了コード0、点数上限、時間切れ、最小刻み枯渇は数値受入やfold確定を意味しない。旧 `passed=false` を書き換えない。

## 1. 環境と全件試験

既存のRust toolchain/lockfileとPythonの固定依存を使う。PowerShell例：

```powershell
git branch --show-current
git merge-base --is-ancestor 5f564621c7396c2a605edbfe5451cd800cf318d3 HEAD
if ($LASTEXITCODE -ne 0) { throw '起点が異なる' }
$TaskRun = 'runs/NEW-nc-fold-PC'
$Py = '.venv/Scripts/python.exe'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo test --locked --release
cargo build --locked --release
foreach ($group in @('python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence')) {
  & $Py tools/run_ci_python_tests.py --group $group --output "$TaskRun/validation/$group"
  if ($LASTEXITCODE -ne 0) { throw "試験失敗: $group" }
}
```

想定：既ビルドで10〜45分程度、初回ビルド・証拠展開はさらに時間とディスクを要する。全群の発見IDが一致し、和集合が欠落・重複なしで全件となること、失敗・error・skip・ignored・filter除外が0であることを確認する。CIの取消・45分打切りを成功にしない。

## 2. 求解しない旧fold近傍の再監査

```powershell
& $Py interim/2026-11-17-nc-fold-stability/extract.py --output work/NEW-fold-replay
& $Py tools/analyze_nc_fold.py --evidence work/NEW-fold-replay/evidence --output work/NEW-fold-analysis
& $Py tools/analyze_nc_fold_failure.py --evidence work/NEW-fold-replay/evidence --output work/NEW-failed-Newton
& $Py tools/nc_fold_analytic.py --output work/NEW-analytic-fold
& $Py -m unittest discover -s tests -p test_nc_fold_stability.py -v
```

想定1〜5分。`extract.py` は旧12分割493,384,680バイトを結合し、分割・全体SHA256・全4,513メンバーを検証する。必要な場とfold近傍のnative材料応答だけを新規展開する。圧縮アーカイブ約0.5 GBと選択展開分の空きが必要で、全3.34 GBはストリームで検証し、全量展開しない。

出力は `fold-results.json`、`low-stress.json`、`native-newton-matrices.json`、`failed-low-stress.json`、`analytic-fold.json`。保存版との比較ではLAPACK/丸め差も記録し、既存値を新しいPC値で上書きしない。解析例は既知の `(u,y,λ)=(0,0,0.95)` のfold、正則区間の非検出、実checkpointの6点＋6点再開を検査する。

## 3. 新しい比例経路の初期1/128・1/256（PC専用）

```powershell
& $Py tools/prepare_nc_fold_series.py --output "$TaskRun/series"
cargo build --locked --release --example nc_history_trace
foreach ($denom in @(128,256)) {
  $p = "$TaskRun/series/$denom"
  & target/release/examples/nc_history_trace.exe trace "$p/mesh.json" "$p/ground.json" "$p/stages.json" "$p/history"
  if ($LASTEXITCODE -ne 0) { throw "履歴保存の異常終了: $denom" }
  # 元入力のコピーであり、checkpointを編集する操作ではない。
  Copy-Item "$p/mesh.json" "$p/history/mesh.json" -ErrorAction Stop
}
& $Py tools/prepare_nc_fold_observer.py --output "$TaskRun/observer"
$Probe = "$TaskRun/observer/target/release/examples/material_diagnostic.exe"
foreach ($denom in @(128,256)) {
  $p = "$TaskRun/series/$denom"
  & $Py tools/nc_fold_path.py --map native --fresh-native "$p/history" --binary $Probe --step-size 0.00002 --count 128 --budget-seconds 1800 --output "$TaskRun/native$denom-coarse"
  & $Py tools/nc_fold_path.py --map native --fresh-native "$p/history" --binary $Probe --step-size 0.00001 --count 256 --budget-seconds 1800 --output "$TaskRun/native$denom-fine"
  & $Py tools/analyze_nc_fold_path.py --path "$TaskRun/native$denom-coarse" --output "$TaskRun/native$denom-coarse-audit"
  & $Py tools/analyze_nc_fold_path.py --path "$TaskRun/native$denom-fine" --output "$TaskRun/native$denom-fine-audit"
}
```

想定：全前史の各経路5〜60分、観測器の初回ビルド2〜10分、局所128/256点は各1〜10分程度。実機差が大きく、この時間は成功保証ではない。時間打切り・前史の未完・近傍の保存点不足なら未確定を残す。`FreshNative` は最高受理点が撤去工程の0.94〜0.96にない入力、物理条件差、初期刻み以外の変更を拒否する。

`series.json` の2条件の差は `stages[2].initial_increment` だけ。材料中点法・最小刻み・差分幅・Newton上限・力/線形許容値・幾何・接触・荷重は既存値を保持する。前史の「初期・載荷・除荷」が2条件で一致することを、`accepted.jsonl` と `run.json` の全GP・接触履歴で確認する。既存の非比例mean/spatial系列をこの比較へ代入しない。

`nc_history_trace` は各受理点を実際に `Checkpoint::seal` → JSON → `resume` しており、入力・実行版の封印を維持する。`history/accepted.jsonl` に全受理場、`history/run.json` に全試行、各工程の途中checkpointに元の生産保存点が残る。`initial-replay.json` には最後の増分を幾何から再構成した全活動GPとの差を記録する。保存された確定場自体は置き換えない。

観測器は現在の `src` の隔離コピーへ既存の読取専用観測関数だけを追加し、`observer.json` に元ソースとバイナリのSHAを保存する。本番の `src`、生産CLI、構成則は変更しない。`nc_fold_path` は既存 `continue_path.py` の移動接線平面法を呼び出す。成分別解放・追加荷重・追加拘束はない。

## 4. 入力だけからの独立前史400/800（PC専用）

```powershell
foreach ($steps in @(400,800)) {
  & $Py tools/low_preload_nc_reference.py --output "$TaskRun/reference$steps" --steps $steps --method rk4 --budget-seconds 7200
  & $Py interim/2026-11-09-nc-pseudo-arclength/independent_arc.py --reference "$TaskRun/reference$steps" --output "$TaskRun/rate-arc$steps" --step-size 0.00002 --count 128
  & $Py tools/nc_fold_path.py --map independent --reference "$TaskRun/reference$steps" --arc "$TaskRun/rate-arc$steps" --step-size 0.00002 --count 64 --material-steps 4 --budget-seconds 1800 --output "$TaskRun/independent$steps-coarse"
  & $Py tools/nc_fold_path.py --map independent --reference "$TaskRun/reference$steps" --arc "$TaskRun/rate-arc$steps" --step-size 0.00001 --count 128 --material-steps 4 --budget-seconds 1800 --output "$TaskRun/independent$steps-fine"
  & $Py tools/analyze_nc_fold_path.py --path "$TaskRun/independent$steps-coarse" --output "$TaskRun/independent$steps-coarse-audit"
  & $Py tools/analyze_nc_fold_path.py --path "$TaskRun/independent$steps-fine" --output "$TaskRun/independent$steps-fine-audit"
}
# 独立の材料内4→8は前史と弧長刻みを固定して別出力へ。
& $Py tools/nc_fold_path.py --map independent --reference "$TaskRun/reference400" --arc "$TaskRun/rate-arc400" --step-size 0.00002 --count 64 --material-steps 8 --budget-seconds 1800 --output "$TaskRun/independent400-material8"
& $Py tools/analyze_nc_fold_path.py --path "$TaskRun/independent400-material8" --output "$TaskRun/independent400-material8-audit"
```

想定：前史各5〜120分、局所各1〜10分。独立側にはnativeの変位・応力・pc/e・保持荷重・受理係数を渡さない。閉鎖固着枝の範囲外、元の降伏ドリフト・線形残差検査による停止ならそのまま拒否として残す。`independent_arc.py` の既存180秒制限による停止も成功へ変えない。空の `samples.json` や不受理起点をnative保存場で補完しない。

## 5. 新規診断checkpointからの実再開

3または4と**同じソース・バイナリ・ds・材料内分割・入力**で、例えば直接64点と32＋32点を比較する。次は独立側の例。

```powershell
& $Py tools/nc_fold_path.py --map independent --reference "$TaskRun/reference400" --arc "$TaskRun/rate-arc400" --step-size 0.00002 --count 32 --material-steps 4 --budget-seconds 1800 --output "$TaskRun/restart-first"
& $Py tools/nc_fold_path.py --map independent --reference "$TaskRun/reference400" --arc "$TaskRun/rate-arc400" --step-size 0.00002 --count 32 --material-steps 4 --budget-seconds 1800 --resume "$TaskRun/restart-first/restart.json" --output "$TaskRun/restart-second"
& $Py -c "import json,pathlib; r=lambda p:[json.loads(s) for s in pathlib.Path(p).read_text().splitlines()]; root=pathlib.Path('$TaskRun'); a=r(root/'independent400-coarse/samples.jsonl'); b=r(root/'restart-first/samples.jsonl')+r(root/'restart-second/samples.jsonl'); assert a==b; print('all fields exact:',len(a))"
```

`restart.json` は実行中の受理点で作られる。制御方式 `common_lambda_moving_tangent_plane/1`、λ、枝方向、弧長尺度、全GP・接触履歴、全荷重成分、モデル・ソース・実行版・dsを束縛する。旧保存行を編集して新checkpointを作らない。Linuxの封印へWindows版exeを付け替えない。`samples.jsonl`、`gp-tangents.jsonl`、`linear-systems.jsonl`、`geometry.json`、`original-load.json`、`evidence-sha256.json`を同じフォルダに保存する。再開では元のcheckpointバイトと受理行、これらのsidecarのSHAも照合する。checkpoint単独のコピーを再封印して代用しない。投影/B行列・重みのhash、Python/NumPy/SciPy/プラットフォームも束縛する。再開の保存行全体と実線形系は、試行/出力パス等の来歴を数値内容と区別して比較し、省略した場がある場合は明記する。

## 6. 比較表・判定基準

```powershell
& $Py tools/compare_nc_fold_reports.py --reports "$TaskRun/native128-coarse-audit/fold-results.json" "$TaskRun/native128-fine-audit/fold-results.json" "$TaskRun/native256-coarse-audit/fold-results.json" "$TaskRun/native256-fine-audit/fold-results.json" "$TaskRun/independent400-coarse-audit/fold-results.json" "$TaskRun/independent800-coarse-audit/fold-results.json" --output "$TaskRun/comparison.json"
```

比較器は数値表を作り、欠けた系列を補わず、自動でNC-01を完了にしない。新しいλピーク許容値を作らない。以下をすべて照合してから、人が結論を追記する。

- 元の入力・前史・制御方式を固定。各実装自身の切断前場から `Ffinal+(1−λ)H` を別構成し、拘束/非活動DOFを含む全成分を比較。nativeは「撤去前の外力−撤去土内力−消失接触力−最終外力」と非活動節点ゼロの元規則を保つ。拘束点の反力を保持荷重に混入しない。
- 力残差1e−8、線形残差1e−10。実際の拡大系残差と力ブロック残差を別々に合格させる。接線の同次方程式は行列・ベクトル尺度による後退誤差も保存する。弧長拘束は既存1e−12 m。
- **受理した連続枝**の実非対称Kで左右特異方程式、σ最小と次の特異値の分離、`wᵀH`、向きを継承した `dλ/ds` の正→負反転、拡大系の正則性を確認。SVDの算術的rank閾値は `n*eps*σmax` であり、平衡許容値の緩和ではない。LU失敗や最小刻み到達だけではfoldとしない。
- 本体、入力由来独立側、元の比例1/128・1/256、局所弧長の刻み半減、独立前史/材料内分割を別表で比較。ピークだけでなく全GP・接触履歴、モード、荷重の差を残す。差が未収束なら未確定。
- 全GPのp′・`K=(1+e)p′/κ` と位置、弱モードとの相関を保存。診断しきい値は報告用で、材料点へ下限を入れない。低応力化が観測されない場合もそのまま記録する。
- 受理点列が不足、力/線形監査失敗、枝の非滑らかな切替、入力/前史不一致なら未確定。局所foldが確定しても、別の非接続解の不存在・動的/熱力学的安定性・全撤去の受入へ拡張しない。NC-01/P6は未完のまま。

長い空間系列や長い下降枝は今回の判定に代入しない。必要なら既存の入力生成・メッシュ系列を別PC計算として保存し、元の比例荷重と各実装自身の前史を維持する。成分別解放、追加拘束、正則化、粘性、p′下限、許容値緩和を試す手順は本書に含めない。
