# B-03クラウド追補の小型再現

GPD WIN Max 2、使えるメモリ約48GBを想定。**PCでは未実行**。
時間/メモリはPC未測定の目安で、合格条件ではない。
旧証拠を上書きせず、出力先は存在しない新しいworkディレクトリとする。
新たな768/6144要素、40GB超、1経路1日超の計算は計画しない。
元の次の土h水準は**規模上限のため確認不可・未達**。許容値を変えない。

## 保存証拠だけの検証

リポジトリ直下のPowerShell、Python3.13とrequirements-reference.lock.txt、Rust1.99.1。
旧封印の生バイトを保つためgit core.autocrlf=falseで取得する。
新しい文書末尾検査は実CRLFもLFに正規化する。数値/封印のSHAには正規化をしない。

```powershell
$B='interim/2026-12-27-b03-offset-kink-cloud-followup'
$Old='interim/2026-12-26-b03-offset-kink-cause-diagnosis'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
$env:CARGO_BUILD_JOBS='2'
$env:RUST_TEST_THREADS='2'
python "$B/b03follow_evidence.py" --verify $B --extract work/b03follow-pc-sealed
python "$Old/b03cause_evidence.py" --verify $Old --extract work/b03follow-pc-prior
python "$B/b03follow_report.py" --rotation work/b03follow-pc-sealed/b03follow-rotation --boundary work/b03follow-pc-sealed/b03follow-boundary --work work/b03follow-pc-sealed/b03follow-work-final --precision work/b03follow-pc-sealed/b03follow-precision --output work/b03follow-pc-reported.json
python -m unittest discover -s tests -p test_b03_offset_kink_cloud_followup.py -v
```

想定1–3分、1GB未満（展開ディスクは余裕を持って2GB）。
見る値：全part/entry SHA、45＋45解、全零モード6/余剰0、元の回転絶対量を超えた13件、
Decimalとの違い、5,859ガード構築と境界表、46自由会計/30土・板状態×2演算子、全体passed=false。
同じ環境の固定順集計ではmeasurements.jsonをバイト再現する。
別OSで最終ビットが変わったら差・演算環境を別記し、新しい許容値で吸収しない。

## 新規に同じ小型診断を実行

```powershell
cargo build --locked --release
python "$B/b03follow_build.py" --output work/b03follow-pc-build
python "$B/b03follow_rotation.py" --probe work/b03follow-pc-build/b03follow_probe.exe --output work/b03follow-pc-rotation
python "$B/b03follow_boundary.py" --probe work/b03follow-pc-build/b03follow_probe.exe --output work/b03follow-pc-boundary
python "$B/b03follow_precision.py" --rotation work/b03follow-pc-rotation --output work/b03follow-pc-precision
python "$B/b03follow_work.py" --prior-extracted work/b03follow-pc-prior --probe work/b03follow-pc-build/b03follow_probe.exe --output work/b03follow-pc-work
python "$B/b03follow_report.py" --rotation work/b03follow-pc-rotation --boundary work/b03follow-pc-boundary --work work/b03follow-pc-work --precision work/b03follow-pc-precision --output work/b03follow-pc-observations.json
```

build想定10–30分/2–6GB、診断全体1–10分/2GB未満、出力ディスク2GBを目安にする。
最大18 Tri6・49節点294DOF。45本体＋45独立の自由全体系、8回のDecimal代数解。
ガードは構築のみであり施工完走件数には数えない。
仕事は前回の同じ保存値を使い、施工経路を解き直さない。

見る値：同じQで移した荷重、剛体6基底・全固有値・残差・板生値・応力テンソル・仕事、
条件数と50桁/80桁の一致、異厚の比較で使う「現在」の板厚、境界の最後の受理/最初の拒否、
作用面と設置参照項。誤った腕/偶力省略を本体に実装する操作はない。
元の回転1e−9等の量を超えた観測を記録したままにし、値を丸めて合格へ変えない。

## 全件ソフトウェア試験と台帳監査

```powershell
python "$B/b03follow_validation.py" rust --output work/b03follow-pc-full/rust
python "$B/b03follow_validation.py" python --output work/b03follow-pc-full/python
python "$B/b03follow_validation.py" collect --output work/b03follow-pc-full
python interim/2026-12-06-p6-01-status-reconciliation/audit_status.py --root . > work/b03follow-pc-p601.json
```

rustにはfmt、probe rustfmt、clippy --all-targets、build --all-targets、全ターゲット試験、doc試験を含む。
pythonにはtoolsを含む全7群を含む。PowerShell自身を必要とする既存試験も実行する。
既存の必要アーカイブを欠落させず、skipで回避しない。
想定1–4時間/2–16GB。既存証拠の取得・展開用ディスクは70GB以上の空きを推奨する。
これは既存ソフトウェア全件の再生を含む目安で、新しい大規模FEM経路の計画ではない。
発見ID・実行ID、成功/失敗/error/未実行/skip、SOURCE_NOT_INVENTORIED=0、integrity_okを確認する。
ソフトウェア試験成功とB-03数値受入falseを区別する。

## 未実装の継手・Native表示

今回の3点をPCへ移送する必要はない。小型再現は任意。
非ゼロ腕の共有ID折れ目/段差、片側e=0、e段差、法線だけ/eだけの反転、範囲外e、曲がった中間節点、
コンパクト/Native準備の現行拒否を解除するコマンドは存在しない。
ユーザー判断を要する接続・回転・作用面・設置参照/所有・認証・表示仕様を選ぶ前に、
文献値に合わせて支持や物性を変える手順は書かない。
将来の実装後の確認も、最初は今回以内の小型パッチから始める。
上限内で元の基準を満たす証拠を出せなければ未達のまま。
実CAD・Native表示は対象入力と表示仕様の確定待ちで、今回クラウド/PCの実行は0件。
push後のWindows CIは未確定。
