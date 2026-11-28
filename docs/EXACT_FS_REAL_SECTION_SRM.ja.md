# 正確なFsの証拠検査と実断面SRMのユーザーPC手順

起点は未マージPR #13の `9f6c228d3b40d02b32b2877f0ed3d0feb9279768`。その変更・数値基準・許容値・参照値・失敗記録を保持した追加です。新しい検証CLIと準備器を追加し、生産の材料計算、探索制御、工程ゲートの数値経路を変更していません。

## 今回実測した範囲

**一様応力・単軸15 kPa・c=10 kPa・φ=ψ=0・直線6 Tet10（27節点）のケースに限り、宣言絶対幅0.002で正確なFsの区間を実測しました。** Fsの値は区間 **[1.3320312500000002, 1.3337890625000002]**、幅 **0.0017578124999999556** です。理論比4/3を含みます。中点の点値は報告しません。

完成した未低減施工からのepoch探索と、元入力から各率を再生するfresh探索を別々に実行し、同じ区間でした。下端は完成した材料経路の数値釣合いを、上端は実際に低減した材料と元外力に結び付く平均応力の必要条件による排除を、同じ実行ファイルで再監査しました。釣合いと排除に矛盾はありません。連続体安定性、設計安全率、斜面Fs、実断面Fsの認定ではありません。

[測定・全コマンド・SHA-256](../interim/2026-11-08-exact-fs/uniaxial-measurement.json)と[元Run・全試行・入力・判定のZIP](../interim/2026-11-08-exact-fs/uniaxial-evidence.zip)を保存しました。アーカイブには実行ファイルを含めません。保存結果の再監査は生成時と完全に同じ実行ファイルを必要とします。別のビルド・OSでは小ケースを新規出力先に再実行してください。

過去の一様応力区間 **[1.3250000000000002, 1.3531250000000001]** は変更していません。今回の測定は別記録です。Tet4の6/48要素の過去記録も保持し、今回の6 Tet10の結果へ置き換えません。

## 検査条件

`verify-exact-fs` と `verify-exact-epoch-fs` は求解を進めません。別途保持した元入力・SPEC・結果ファイルSHA-256に結び付けて結果を検査します。

1. 下端の材料経路が完成し、その率で外力・活動領域・拘束・降伏・力の数値監査が成立する。
2. 上端は平均応力または仮想仕事を入力から再評価して独立に排除できる。保存された排除フラグだけでは受け入れない。
3. 試行の順序・率・要求率の実評価・停止理由・区間の添字と実幅が一致し、釣合いと排除が矛盾しない。
4. 幅が別途指定したSPECの `absolute_factor_tolerance` 以内で、探索が `factor_tolerance` で終わっている。

欠落・改変は終了コード1。証拠は整合していても試行上限、未完の材料経路、不収束、物理的上端なし等で未確定なら、`exact_fs_for_this_case=false` の判定JSONを新規保存して終了コード2です。条件を満たす場合のみ0です。元の探索・ジョブの終了コード0は「計算終了」で、正確なFsの合格とは別です。検証器は中点フィールドを出しません。

epoch結果は既存保存ジョブと同じ材料checkpointの封印・元Run・目標パラメータ・完成状態・全試行の再監査を使います。fresh結果は元施工接頭部と各率のGround/Caseを照合し、保存された完成Runの終点を再監査します。求解例外は判定不能としてのみ保持し、例外の文章から上端を作りません。結果全体の保存時SHAを別に保持することで、監査対象外の過去フレームや計測時間の改変も拒否します。SHAは署名ではないため、元の入力・SPEC・保存時SHAと生成実行ファイルを一緒に保全してください。

工程ゲートは要求率の釣合いと適切な探索終了で継続する機能です。たとえば `no_physical_upper` でも要求率のゲートは継続できますが、正確なFsには上端が必要なので、この報告器では未完了です。工程ゲートのpassedを正確なFsへ読み替えません。

## 既存の斜面証拠（再求解・上端再計算なし）

[公開JSONからの要約](../interim/2026-11-08-exact-fs/archived-evidence-summary.json)に原ファイルのSHAも保存しました。

|対象|公開記録|意味|
|---|---|---|
|固定1200 Tet10、F=1.45|全32400 GP、独立外力仕事173.11199359041558、容量164.68252295518616、差8.429470635229421 kN m|独立な離散釣合い排除。正確なFs・設計安全率ではない|
|固定1200 Tet10、F=2|全32400 GP、独立外力仕事217.9687881634925、容量164.80726244312086、差53.16152572037163 kN m|離散排除。斜面探索区間を完成させる下端・宣言幅の証拠ではない|
|Griffiths比較条件1.3500000000000003|過去の比較に指定した低減係数|物理的破壊Fsではない|

[1.45の原JSON](../interim/2026-10-30-Root-SRM-and-MCC-recovery/quadratic-proposals/v5/F1p45/independent-PASS.json)、[2の原JSON](../interim/2026-10-30-dry-SRM-search-and-slope-work-bound/slope1200-work-audit/terminal.json)を参照。クラウドでは今回、斜面1200要素級SRM探索も仮想仕事上端の再計算も行っていません。

Griffiths中央面L2 **15.582／15.886／15.593／15.672％の不合格**、連続側の判定不能は保持します。柱のTet4最細重力場 **6.828％不合格**、連続推定 **1.034％合格**、Tet10の5％合格は、それぞれの既存ケースの結果です。既存MCC掘削の空間差、押出し応力・振動モニター・有限幅奥行き応力の判定不能も保持します。[Griffiths](GRIFFITHS_2D_COMPARISON.ja.md)、[柱](PLANE_STRAIN_2D.ja.md)、[奥行き診断](DEPTH_SPATIAL_DIAGNOSTICS.ja.md)の判定は変更しません。

## 実断面入力の検索結果と準備器

起点の追跡2659ファイル、4 MiB以下のJSONの入力schema、ZIP 101個のメンバー名、実断面/CAD関連文書・入力候補を調べました。[検索台帳](../interim/2026-11-08-exact-fs/input-inventory.json)を保存しています。既存SRM仕様は一軸検証、数値失敗検証、Griffiths1200の準備入力でした。実CADの288 Tet10記録は二層地盤の乾燥MCCであり、乾燥MC実断面SRM入力として代用できません。測量断面等の出典と対応する実断面SRM入力の組は特定できませんでした。ZIP名の台帳は、全アーカイブの全ファイル内容を検証した証明ではありません。

**実断面SRMは未実行・未完了です。** `tools/prepare_real_section_srm.py` は、ユーザーが指定する形状・材料・施工・SPECと断面出典文書を既存の保存型fresh探索ジョブへ結び付ける準備器です。入力検査とジョブ作成だけを行い、探索を開始しません。形状やFsを生成せず、出典の測量上の真正性を自動認定もしません。テストは明記した一軸検証入力を使い、実断面実行の証拠にはしません。

準備記録の `search_executed_by_preparer`、`exact_fs_for_this_case`、`real_section_srm_complete` はすべてfalseです。検証時は出典文書・ジョブmanifest・全入力・保存実行ファイルを照合し、完成ジョブと保存時の結果SHAを確認してから上記の報告器を実行します。未実行、結果欠落、改変、幅未達を合格にしません。

## ユーザーPC：実断面の準備・探索・検証

以下はWindows PowerShell例です。実際の実断面入力が揃ったPCで実行します。今回クラウドでこれらの実断面コマンドは実行していません。先にこのブランチを取得し、Rust 1.99.1とPython 3.13を用意します。

```powershell
git switch codex/plastic-spatial-convergence
git pull --ff-only
git merge-base --is-ancestor 9f6c228d3b40d02b32b2877f0ed3d0feb9279768 HEAD
cargo build --locked --release
python -m pip install -r requirements-reference.lock.txt
$srmSolver = (Resolve-Path target/release/soil-fem-3d.exe).Path
$srmInput = "D:/ユーザーが保有する実断面入力"
$srmBundle = "D:/計算記録/実断面SRM-新規01"
python tools/prepare_real_section_srm.py prepare --solver $srmSolver --mesh "$srmInput/mesh.json" --ground "$srmInput/ground.json" --case "$srmInput/case.json" --spec "$srmInput/specification.json" --source "$srmInput/断面出典.pdf" --section-id "出典で特定した断面ID" --output $srmBundle
if ($LASTEXITCODE -ne 0) { throw "準備未完了" }
```

出典・断面ID・寸法・座標・奥行き・材料・荷重・施工接頭部を利用者が確認してください。既存SPECに必要な材料ID、ψ/引張方針、要求率、最大率、宣言絶対幅、最大試行数、工程数を明示します。許容幅や数値釣合い条件を結果に合わせて緩めません。元の探索は乾燥MC・対応済み工程・明示したゼロ変位拘束に限定され、水・構造・JOINT・Davis探索・MCC/IDHからの変換などを拒否します。

```powershell
$srmJob = "$srmBundle/job"
$srmJobSolver = "$srmJob/solver/soil-fem-3d.exe"
& $srmJobSolver verify-srm-search-job $srmJob
& $srmJobSolver run-srm-search-job $srmJob
& $srmJobSolver srm-search-job-status $srmJob
python tools/prepare_real_section_srm.py verify $srmBundle "$srmBundle/exact-fs-assessment.json"
if ($LASTEXITCODE -ne 0) { Write-Output "正確なFsは未完了。停止理由と証拠を確認する" }
```

出力先はすべて新規です。保存用solverはLinuxでもファイル名が `.exe` になります。Linuxでは同じコマンドを実行パス表記に置き換えて使えます。準備後は同じ保存solverを使い、探索の途中で再ビルド版へ置き換えません。検証器の0かつ `exact_fs_for_this_case=true` が得られた場合だけ、**指定した実断面・工程・離散モデル・宣言幅**に限る区間を報告します。Fsの点値、連続体崩壊、設計安全率の自動認定はしません。

長時間計算の停止は別ターミナルから `stop-srm-search-job JOB`。終了コード2は受理境界での休止です。継続は既存 `prepare-srm-search-recovery OLD_JOB NEW_JOB` で新しいジョブを作り、保存solverで `run-srm-search-job NEW_JOB` を行います。元ジョブと元結果を保存します。回復ジョブを新しい実断面準備記録へ黙って差し替えません。回復結果は新ジョブの保存時SHAと元の実断面出典を保持し、以下の直接検証を使います。

```powershell
$srmJob = "D:/計算記録/実断面SRM-回復01"
$srmJobSolver = "$srmJob/solver/soil-fem-3d.exe"
& $srmJobSolver verify-srm-search-job $srmJob
$srmStatus = Get-Content "$srmJob/status.json" -Raw | ConvertFrom-Json
& $srmJobSolver verify-exact-fs "$srmJob/inputs/mesh.json" "$srmJob/inputs/ground.json" "$srmJob/inputs/case.json" "$srmJob/inputs/specification.json" "$srmJob/report.json" $srmStatus.report_sha256 "$srmJob/exact-fs-assessment.json"
```

保存時SHAを別媒体へ保存します。改変後のファイルからSHAを再計算して「期待値」に差し替えないでください。

## 小ケースの独立した再測定

以下は6 Tet10の一様応力ケースだけであり、実断面の代用ではありません。`runs`を作成し、同じビルドで一連のコマンドを実行します。元の記録は上書きしません。

```powershell
New-Item -ItemType Directory -Force runs | Out-Null
$srmExample = "examples/exact-fs-uniaxial"
& $srmSolver solve-stages "$srmExample/mesh.json" "$srmExample/ground.json" "$srmExample/case.json" runs/uniaxial-origin-new.json
& $srmSolver search-mc-epoch-srm "$srmExample/mesh.json" runs/uniaxial-origin-new.json "$srmExample/epoch-spec.json" runs/uniaxial-epoch-new.json
$srmReportHash = (Get-FileHash runs/uniaxial-epoch-new.json -Algorithm SHA256).Hash.ToLowerInvariant()
& $srmSolver verify-exact-epoch-fs "$srmExample/mesh.json" runs/uniaxial-origin-new.json "$srmExample/epoch-spec.json" runs/uniaxial-epoch-new.json $srmReportHash runs/uniaxial-assessment-new.json
```

この直接探索では生成直後のSHAを別に保存します。fresh経路は `search-mc-srm MESH GROUND CASE fresh-spec.json NEW_REPORT`、検証は `verify-exact-fs MESH GROUND CASE fresh-spec.json REPORT EXPECTED_SHA256 NEW_ASSESSMENT` です。

## 試験と残件

新規試験はRust6件・Python4件。完成経路と平均応力排除、元入力/SPEC結合、試行欠落・順序・端点改変、偽造排除、不収束の読み替え、材料試行/探索試行上限、幅改変、上端なしゲートとの区別、準備のみ・出典欠落・ジョブ改変・実CLIを検査します。[検証記録](../interim/2026-11-08-exact-fs/README.ja.md)に全件数とWindows CIを記録します。

未完了は斜面の正確なFs、実断面SRM、一般3D空間収束、新斜面h系列の5％、既存MCC掘削の空間差、一般の水・構造・JOINT SRM、Native UI、別PC/GPU、完成配布、P6全体です。Griffiths1200/4800/19200 Tet10、乾燥MC4/8/16・8/16/32、既存MCC274〜14909要素をクラウドで再求解していません。
