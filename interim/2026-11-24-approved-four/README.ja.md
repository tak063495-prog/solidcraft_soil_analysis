# C-01、T7、Kirsch、排水 MCC と C-03 の検証

2026-11-24 10:59 JST 承認分。C-01 の最終比較は issue #44 を含む PR #13 先端 `bb2ef41c5dc5315ca1d17be209b399567db499bb` を起点とする。

## C-01 の接触事象

`--mode friction-events` は既存の `local-active` に事象位置推定を加える。二つの半増分で、閉じた固着点の弾性試行せん断力が Coulomb 限界に達する位置を推定する。開閉には法線ギャップの零点、滑りから固着への遷移には状態の異なる半増分の中点を使う。座標は再試行の提案であり、各再試行で全体を釣り合わせる。端点からの分解能は元の刻みの 1/64 と最小刻みの大きい方。12 Tet10 の各工程の試行上限は 2,065、固定予算の既定は 2,000。

独立系で最初の初期釣合い増分が停止した N512・N1024 のログを保存した。境界変位関数 `offset()` の呼び出しと同じスコープにあった局所変数を `event_offset` とした。[N512 の停止記録](retained-failures/independent-512.log)・[N1024 の停止記録](retained-failures/independent-1024.log)

|実装・モード|3水準の共通座標|符号付き塑性すべりの次数／余弦|塑性散逸の次数／余弦|工程内最多試行 N256／512／1024|
|---|---:|---:|---:|---:|
|本体・従来比較|1,792|2.1102／-0.3866|0.9791／+0.9995|576／993／1,906|
|本体・事象位置推定|1,792|1.9061／-0.2791|0.9808／+0.9996|614／1,058／1,951|
|独立・従来比較|1,792|1.9438／-0.4991|0.9853／+0.9996|568／994／1,861|
|独立・事象位置推定|1,792|2.0044／-0.3178|0.9842／+0.9997|620／1,071／1,920|

両系とも連続20場のうち13場で次数・余弦がともに正、6場で両差がゼロだった。符号付き塑性すべりの余弦は両系で負のままで、全場の刻み収束は**未達**。[同一先端の比較対照](c01-control.json)・[事象位置推定後の全場、工程別試行、ストリーム照合](c01-fields.json)・[範囲別の判定](scope-rows.json)

## T7 固定ヒンジ深いアーチ

`arch.py --extrapolate` は直近最大4個の受理平衡点から Newton 初期値を外挿し、現在の弧長補正面へ戻す。受理判定と単純折返しの識別条件は従来値。荷重区間は公開図からの読取り 885〜905。[公開例](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEEXARefMap/simaexa-c-snapbuckling.htm)

|梁要素|受理／棄却|折返しと判定|
|---:|---:|---|
|12|66／0|単純折返し 935.41250。区間外・未達。|
|60|99／15|単純折返し 898.40341。区間と荷重反転のみ限定合格。|
|96|145／81|受理荷重 896.91921 → 897.36420 → 896.95378。反転を観測、単純折返しは分類 D・未達。|

全曲線の定量受入は未達。[荷重・変位曲線](arch.json)・[全試行](arch-trials.json.gz)

## Kirsch 円孔の段階掘削

半径 1 m、外周半径 10 m。初期体とコア撤去後の解放率 0、0.5、1 で既存と同じ境界荷重・材料を使用した。孔の近傍へ等比に半径方向の節点を集め、周方向にも分割した。平面ひずみ、厚さ方向は1層。外挿には各軸から内側2半径区間・2角度区間の24個の Tet4 積分点を使う。各点の周方向応力を局所極座標で求め、体積重みの最小二乗で定数・半径の1次と2次・軸からの角度の2次の項を当てはめた。孔壁かつ軸上での定数項を評価した。解析応力は回帰に使用していない。[公開閉形式](https://docs.itascacg.com/itasca950/3dec/docproject/source/examples/CylindricalHoleInAnIEM.html)

|周方向 × 半径方向|残存 Tet4|従来の軸隣接 GP 平均、水平／鉛直 MPa|孔壁外挿、水平／鉛直 MPa|
|---|---:|---:|---:|
|4 × 2|48|35.596／22.456|38.578／19.548|
|8 × 4|192|46.344／22.342|58.564／21.451|
|16 × 8|768|57.779／18.555|72.315／16.656|
|24 × 16|2,304|65.841／16.360|75.642／14.661|

従来の厚さ2層・96 Tet4 は 35.618／22.556 MPa。理論値 75／15 MPa に対し、最細外挿の符号付き誤差は +0.642／−0.339 MPa（+0.86%／−2.26%）。初期体を含む全16解の自由力・線形相対残差は最大 9.94×10⁻¹⁴ 未満。釣り合いは限定合格。公刊の場誤差許容幅は未確定で、場精度の定量受入は未達。[各解放段階・回帰の条件数](kirsch-wall.json)

## 一定 G の排水 MCC

公開例の G = 250 kPa、初期平均有効応力 5 kPa、初期先行圧密圧力 8／40 kPa を用いた。材料内刻み 10⁻⁴／5×10⁻⁵、12／96 Tet4 と OCR 1.6／8 の8ケースで、各512せん断増分を記録した。軸ひずみ 0.5 の p′、q、v は全ケースで既存の 2% 基準内。最大相対誤差はそれぞれ 0.0842%、0.2475%、0.0160% 未満。**終端値のみ限定合格**。軸ひずみに対する全履歴の点別受入幅は公刊資料から確定できず、**参照なし（全文待ち）・未達**。[公開例・表1・終端理論式](https://docs.itascacg.com/flac3d700/flac3d/zone/test3d/VerificationProblems/TriaxialCompressionTest/TriaxialCompressionTest.html)・[各増分と終端値](mcc-constant-g.json)

## C-03 の局所検証

局所点の外部参照との定量検証を独立した未達行に加えた。材料点には Sysala・Čermák の MC リターンマッピング、接触点には 3DEC の閉じた傾斜き裂の載荷・除荷問題が候補。6親入力での GP273 の分類一致は既存の限定結果。湿潤経路全23場と塑性仕事の定量受入は **参照なし（全文待ち）・未達**。[範囲別台帳](../2026-12-06-p6-01-status-reconciliation/ledger-ac-20261123.ja.md)・[局所原票と公開参照](c03-local.json)

## 再現

リポジトリのルートから、固定版の Rust・Python 依存関係を使用する。BLAS のスレッド数を1にし、出力には新規ディレクトリを指定する。

```bash
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
cargo build --locked --release --bin soil-fem-3d --example c01_methods --example t9d_mcc
out=interim/2026-11-24-approved-four
raw=runs/approved-four
mkdir -p "$raw"
for side in native independent; do
  for n in 256 512 1024; do
    for mode in local-active friction-events; do
      python interim/2026-11-24-c01-event4-retry/route.py "$raw" --solver "$side" --steps "$n" --mode "$mode"
    done
  done
done
python "$out/c01_fields.py" "$raw" "$raw/c01-fields.json"
python "$out/c01_fields.py" "$raw" "$raw/c01-control.json" --mode local-active
python interim/2026-11-24-t9d-benchmarks/arch.py "$raw/arch.json" --extrapolate
python "$out/kirsch_wall.py" target/release/soil-fem-3d "$raw/kirsch"
target/release/examples/t9d_mcc --constant-g > "$raw/mcc-constant-g.json"
```

## Linux の試験結果

実装コミット `8efe778e96178ac5cc410919710978a026999c81` で fmt、Clippy（全ターゲット、警告をエラー扱い）、Rust 全ターゲット 541 件が合格。doctest は対象0件、終了コード0。CI の Python 7群は最終実行で計1,112件が合格し、Linux 専用7件も合格した。mean-evidence と spatial-evidence はそれぞれ4件すべて合格した。[件数・終了コード・生ログの SHA-256](tests/linux.json)

初回は PowerShell 不在で各 core 群1件がスキップされた。Microsoft の署名元リポジトリで掲示されたパッケージの SHA-256 を照合し、PowerShell 7.6.9 を作業環境に展開した。再試行では core-1 の SIGINT 試験で生成途中の空の PID ファイルを1回読み、`ValueError` が発生した。同じソースで当該試験単体と core-1 全548件が合格した。[初回 core-0](tests/python-core-0-initial.log.gz)・[初回 core-1](tests/python-core-1-initial.log.gz)・[PID 読み取り失敗](tests/python-core-1-powershell-first.log.gz)・[最終 core-1](tests/python-core-1.log.gz)
