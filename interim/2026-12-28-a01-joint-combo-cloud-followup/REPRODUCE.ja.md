# A-01クラウド追補の小型再現

GPD WIN Max 2、使えるメモリ約48GB。**PC実行はしていない**。
見込みはWindows未実測の予算であり、受入条件ではない。旧出力を上書きせず、新しいworkパスを使う。
新規768/6144要素、40GB超、1経路1日超の計画は作らない。
上限超えの一般化は**規模上限のため確認不可・未達**。元の許容値を維持する。

## 取得と全件試験の準備

Python3.13、Rust1.99.1、requirements-reference.lock.txt、PowerShell。
生証拠のバイトを保つ取得はcore.autocrlf=false。文書検査だけuniversal newlineで読む。
PR #13の同じブランチを使い、マージしない。

```powershell
git switch codex/plastic-spatial-convergence
git pull --ff-only
git config core.autocrlf false
# core全件に必要なアーカイブを除外しない。旧postfold除外だけを保持する。
git sparse-checkout set --no-cone '/*' '!/interim/2026-11-10-nc-postfold-branches/evidence.tar.gz.part*'
python -m pip install -r requirements-reference.lock.txt
$env:CARGO_BUILD_JOBS='2'
$env:RUST_TEST_THREADS='2'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
$D='interim/2026-12-28-a01-joint-combo-cloud-followup'
```

全件試験は既存の約27GBの証拠を読み、取得/展開/ビルド用に70GB以上のディスク余裕が必要になる。
これは既存保存証拠の保管量であり、新しい大規模FEM計算の計画ではない。
partial cloneでは作業ファイルの原Git blob SHAと封印SHAの照合に加え、P6台帳の固定コミット/出典blobも取得する。
欠落をskip・参照値の置換・manifestの書換えで回避しない。生ファイルだけ取得した場合はGitオブジェクトもあると思い込まない。
完全な通常cloneなら追加の固定コミット取得は不要。shallow/partialで履歴全体を取得し直す必要はなく、台帳にある正確なOIDを限定取得できる。

## 保存証拠の検査と再集計

```powershell
python "$D/a01follow_evidence.py" --verify $D --extract work/NEW-a01follow-sealed
python "$D/a01follow_report.py" --formal work/NEW-a01follow-sealed/formal --extra work/NEW-a01follow-sealed/stacked-extra --separation work/NEW-a01follow-sealed/separation --output work/NEW-a01follow-reported
python -m unittest discover -s tests -p test_a01_joint_combo_cloud_followup.py -v
```

想定1〜5分、1〜4GB。全part/entry SHA、P×R全履歴、仕事の残り、枝別rank、18候補と未採用の追加観測を確認する。
再集計は新しい平衡求解に数えない。別OSで最終ビットが違ったら生差と演算環境を記録し、新しい許容値で吸収しない。

## 同じ小型入力を新規求解

```powershell
cargo build --locked --release
python "$D/a01follow_build.py" --output work/NEW-a01follow-build
python "$D/a01follow_kernel_build.py" --output work/NEW-a01follow-kernel
python "$D/a01follow_batch.py" --output work/NEW-a01follow-formal --build work/NEW-a01follow-build
python "$D/a01follow_batch.py" --output work/NEW-a01follow-extra --build work/NEW-a01follow-build --stacked-offdiagonal
python "$D/a01follow_separation.py" --formal work/NEW-a01follow-formal --output work/NEW-a01follow-separation --contact-probe work/NEW-a01follow-kernel/a01follow_contact_probe.exe
```

12 Tet4（辺共有A/B、板B）と18 Tet4（平行2面）だけ。診断ソースコピーのenv切替はこの専用ビルドだけが読む。
通常CLIの既定や生産ソースを書き換える操作ではない。
全体8分割・各エンジンP×R4通り、部分対照3本、平行2面の4/8/16を含む。
ソースコピーが変わった場合は一致数検査で止める。差分を無条件に当てず、失敗記録を残す。
生成ソースと実行版SHAを別のrunのものに付け替えない。

ビルド想定10〜40分/2〜8GB、通常の診断1〜15分/1〜4GB。
平行2面の独立P＋本体Rは分割50,000に達するため900秒上限の未完了記録となり得る。
全体を1〜2時間以内/8GB以内の予算とし、上限に達したら途中出力を保持して未完了とする。
元の候補18点・下降条件を変更しない。追加cut=18以降は観測だけで、解に採用しない。
counterfactualの解析接線/初期力/世代/刻みは残差・方向の評価だけで、施工の新しい完走経路ではない。

既存12 Tet10の任意範囲確認が必要な場合は、[12-20の小型手順](../2026-12-20-a01-joint-combo-cause-diagnosis/REPRODUCE.ja.md)の
辺共有Bの1ケース（54土節点、324材料GP）から始め、見込み8〜32GB・0.5〜2時間/経路の範囲だけとする。
今回の3点をPCへ必須移管する必要はない。大きいメッシュへ拡大しない。

## 全件ソフトウェア検証とP6台帳

```powershell
python "$D/a01follow_validation.py" rust --output work/NEW-a01follow-full/rust
python "$D/a01follow_validation.py" python --output work/NEW-a01follow-full/python
python "$D/a01follow_validation.py" collect --output work/NEW-a01follow-full
python interim/2026-12-06-p6-01-status-reconciliation/audit_status.py --root . > work/NEW-a01follow-p601.json
```

fmt/clippy/build、全Rustターゲット/doc、toolsを含むPython全7群。
全発見ID/実行ID、成功/失敗/error/未実行/skipを記録し、単独再確認を二重加算しない。
旧証拠を含む全件は数時間になる場合がある。重い新規経路をCIへ追加しない。
SOURCE_NOT_INVENTORIED=0、integrity_ok=trueを確認する。push後のWindows CIは未確定。

A-01は未達、passed=false。
