# A-02クラウド追補の小型再現

GPD WIN Max 2、使用可能約48GBの上限を守る。**PC未実行**。
同じ81節点・24 Tet10・648材料GP・32接触GPだけを使い、N=2/4/8、最小接触2/4/8。
768/6144要素、40GB超、1経路1日超の計画は作らない。
上限超の一般化は**規模上限のため確認不可・未達**。許容値・既存保存値を変えない。

## 取得・前提

PR #13の同じブランチを取得し、マージしない。新しいwork出力名を毎回使う。
Python3.13、Rust1.99.1、PowerShell、requirements-reference.lock.txtを使う。

```powershell
git switch codex/plastic-spatial-convergence
git pull --ff-only
git config core.autocrlf false
git sparse-checkout set --no-cone '/*' '!/interim/2026-11-10-nc-postfold-branches/evidence.tar.gz.part*'
python -m pip install -r requirements-reference.lock.txt
$env:CARGO_BUILD_JOBS='2'
$env:RUST_TEST_THREADS='2'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
$D='interim/2026-12-29-a02-joint-birth-cloud-followup'
```

旧postfoldの既定除外を除いて、全件試験に必要な旧証拠を取得する。
partial cloneは固定コミット/出典blobも必要。生ファイル取得だけでGitオブジェクト取得済みと思い込まない。
Gitの一括履歴取得失敗は試験の数値失敗と分け、必要OIDを限定取得・照合する。
[b383aa0の準備失敗記録](../2026-12-27-b03-offset-kink-cloud-followup/validation/preparation-events.json)も参照。
全件の旧証拠保管・ビルド・展開には大きなディスク余裕が要る（前回目安70GB以上）。
これは新しい大規模FEMの計画ではない。欠落をskipや封印の再生成で回避しない。

## 封印の検査・再集計

```powershell
python "$D/a02follow_evidence.py" --verify $D --extract work/NEW-a02follow-sealed
python "$D/a02follow_report.py" --run work/NEW-a02follow-sealed/formal --output work/NEW-a02follow-reported
python "$D/a02follow_contact_report.py" --run work/NEW-a02follow-sealed/retained/a02follow-contact-exchange --main-run work/NEW-a02follow-sealed/formal --main-report work/NEW-a02follow-reported --output work/NEW-a02follow-contact-reported
python "$D/a02follow_cone_report.py" --cone work/NEW-a02follow-sealed/retained/a02follow-cone-v2 --rounding work/NEW-a02follow-sealed/retained/a02follow-rounding-v2 --output work/NEW-a02follow-cone-summary.json
python -m unittest discover -s tests -p test_a02_joint_birth_cloud_followup.py -v
```

想定1〜5分/1〜4GB。part/全entry SHA、全86主経路（78＋接触活動対照8）、枝位置・差ベクトル・coneを確認する。
再集計は新しい全体求解ではない。別OSで最終ビットが違ったら差と環境を記録し、新許容値を作らない。
文書末尾検査はuniversal newlineで読み、CRLF/autocrlfでも同じ論理末尾を検査する。

## 同じ小型入力の新規求解

```powershell
python "$D/a02follow_build.py" --output work/NEW-a02follow-build
python "$D/a02follow_cone_build.py" --output work/NEW-a02follow-cone-build
python "$D/a02follow_run.py" batch --output work/NEW-a02follow-formal --build work/NEW-a02follow-build
python "$D/a02follow_cone.py" --output work/NEW-a02follow-cone --build work/NEW-a02follow-cone-build
python "$D/a02follow_rounding.py" --output work/NEW-a02follow-rounding --build work/NEW-a02follow-cone-build
python "$D/a02follow_contact_exchange.py" batch --output work/NEW-a02follow-contact-exchange --build work/NEW-a02follow-build
python "$D/a02follow_report.py" --run work/NEW-a02follow-formal --output work/NEW-a02follow-reported
python "$D/a02follow_contact_report.py" --run work/NEW-a02follow-contact-exchange --main-run work/NEW-a02follow-formal --main-report work/NEW-a02follow-reported --output work/NEW-a02follow-contact-reported
python "$D/a02follow_cone_report.py" --cone work/NEW-a02follow-cone --rounding work/NEW-a02follow-rounding --output work/NEW-a02follow-cone-summary.json
```

専用ソースコピーだけにP/R・最小分割・小型観測を適用する。生産ソースの置換はしない。
変換位置は一致件数を検査して止める。生成ソース・patch・実行版SHAと入力の対応を保持する。
既定の最小1・P/R nativeを省略した通常CLIへこの切替は存在しない。

交差対照54主経路とP×R交換24主経路。接触側の実分割24/12を比較するB開閉8経路を加え、主求解は86経路。
本体43経路は観測なしwholeも別求解し、合計129経路として区別して計数する。
coneは公開局所則・投影だけの指定入力診断で、全体施工の新規平衡には数えない。
最初の確認は1経路だけでよい：

```powershell
python "$D/a02follow_run.py" worker --output work/NEW-a02follow-one --build work/NEW-a02follow-build --n 4 --parts 2 --case stress_free-cycle --engine native --predictor native --repartition native
```

Windows未実測の予算はビルド10〜40分/2〜8GB、1経路は5分以内をまず目安とする。
1経路の本体プロセス既定300秒・外側worker960秒で停止記録を保存する。全行列は1〜3時間/8GB以内を予算とする。
未完了なら失敗/途中出力を保持し、完走・scope拒否・物理破壊へ読み替えない。時間予算は合格基準ではない。
24 Tet10のN=2/4/8以外へ拡大しない。新しい接触則、クリップ、cone許容幅、空間参照を導入しない。
残るN=8のP/R×実接触分割と、切替区間の入力二分による全体平衡は今回未実施。
同じ小型入力でクラウド実行可能な残件であり、資源不足としてPCへ必須移管しない。

## 全件試験・P6

```powershell
python "$D/a02follow_validation.py" rust --output work/NEW-a02follow-full/rust
python "$D/a02follow_validation.py" python --output work/NEW-a02follow-full/python
python "$D/a02follow_validation.py" collect --output work/NEW-a02follow-full
python interim/2026-12-06-p6-01-status-reconciliation/audit_status.py --root . > work/NEW-a02follow-p601.json
```

fmt/clippy/build、Rust全target/doc、toolsを含むPython全7群。発見/割当/実行IDを照合し、成功/失敗/未実行/skipを記録する。
単独スモークや再確認を固有件数へ二重加算しない。新しい重い試験をrust-integration-1へ追加しない。
台帳の固定コミット・SHA-256・Git blobを検証し、SOURCE_NOT_INVENTORIED=0を確認する。
新記録を登録しないまま全件成功とはしない。push後のWindows CIは未確定。
A-01、A-03、A-04、B-01〜B-03、C-01〜C-04、NC-01、P6は未達のまま。
A-02は未達、passed=false。
