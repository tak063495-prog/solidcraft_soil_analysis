# C-04 16/64×P/Rの小型再現

PR #13、codex/plastic-spatial-convergence、起点71655b7b09e1895e364ba476b9089075ba7130c0以降。マージしない。
GPD WIN Max 2・利用可能メモリ約48GB。Windows PCでは未実行。
12 Tet4・16節点・12材料GP・6接触GP・13工程、水位32、力学16/64だけ。
768/6144要素、メモリ40GB超、1経路1日超の新しい計画はない。
必要なら **規模上限のため確認不可・未達** と記録し、許容値を変えない。

## 取得・準備

初回全履歴転送を避け、未存在のNEW-c04followへ単一ブランチをfiltered取得する。
既存CIのpostfold除外以外の必要な旧アーカイブを欠落させない。
A-02〜B-02のgit-recoveryと、今回のvalidation/git-recoveryも参照する。
一時ディスクにGitキャッシュ/展開/ビルドの余裕を持たせる。
キャッシュ整理をするなら、旧実ファイルをGit blob/SHAで照合した後の再取得可能な単一blobだけとする。

```powershell
git clone --filter=blob:none --single-branch --branch codex/plastic-spatial-convergence --depth=1 --no-checkout https://github.com/tak063495-prog/solidcraft_soil_analysis.git NEW-c04follow
Set-Location NEW-c04follow
git config --local core.autocrlf false
git sparse-checkout set --no-cone '/*' '!/interim/2026-11-10-nc-postfold-branches/evidence.tar.gz.part*'
git switch codex/plastic-spatial-convergence
$Py='.venv/Scripts/python.exe'
if (-not (Test-Path $Py)) { py -3.13 -m venv .venv }
& $Py -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw '固定依存の準備未完了' }
$D='interim/2027-01-03-c04-transient-water-cloud-followup'
$env:CARGO_BUILD_JOBS='2'
$env:RUST_TEST_THREADS='2'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
```

Rust1.99.1/rustfmt/clippy、Python3.13、PowerShell。HOMEやNC固有設定は変更しない。
読取専用HOMEならXDG_CACHE_HOME/XDG_CONFIG_HOME/XDG_DATA_HOMEだけをworkの書込可能な場所へ向ける。
partial cloneでP6-03固定BASEのtreeと30文書パス/29固有blobも必要。
台帳の固定コミットをdepth=1/blob:none、必要blobをOID指定で取得し、SHAを照合する。
今回の取得器はB-02の方式を継承し、ログとソースを封印に保持する。
同じ起点の復旧を新cloneで使う場合：

```powershell
New-Item -ItemType Directory -Force work/c04follow-preparation | Out-Null
& $Py "$D/c04follow_evidence.py" --verify $D --extract work/NEW-c04follow-reopened
if ($LASTEXITCODE -ne 0) { throw '封印SHA不一致' }
& $Py work/NEW-c04follow-reopened/c04follow-preparation/c04follow_fetch_fixed.py
if ($LASTEXITCODE -ne 0) { throw 'P6固定出典取得未完了' }
git -c fetch.negotiationAlgorithm=noop fetch --depth=1 --filter=blob:none --no-tags --no-write-fetch-head origin e3d830d8d1bdffd727c4ab7271755f52e556006d
if ($LASTEXITCODE -ne 0) { throw 'P6-03固定BASE取得未完了' }
```

## 新規全体求解

出力は未存在のwork先にする。過去の場・失敗・拒否・参照値を上書きしない。

```powershell
& $Py "$D/c04follow_build.py" --output work/NEW-c04follow-build
if ($LASTEXITCODE -ne 0) { throw '隔離ビルド未完了' }
& $Py "$D/c04follow_run.py" --build work/NEW-c04follow-build --output work/NEW-c04follow-formal
if ($LASTEXITCODE -ne 0) { throw '数値未完了。部分場とfailure.jsonを保持' }
& $Py "$D/c04follow_report.py" --run work/NEW-c04follow-formal --output work/NEW-c04follow-observations
if ($LASTEXITCODE -ne 0) { throw '観測未完了。全体求解の完了とは分ける' }
```

計16履歴・7,552受理状態。各経路に元の力/線形基準を使う。
保守的PC予算はビルド10〜90分・2〜8GB、正式系列10〜90分・1〜4GB（Windows未実測）。
1子プロセスの上限3,600秒を延ばして大規模化しない。タイムアウト/部分場は数値未完了として保持する。
通常CLI・生産RustにはC04FOLLOW切替はない。隔離コピーの既定N/Nは保持し、P/Rだけ明示交換する。
接線/差分幅/許容値は交換しない。独立側へ本体の受理解・分割・接線を渡さない。

fields.jsonl、raw/accepted.jsonl、raw/constitutive.jsonl、raw/run.json、独立result.jsonを読む。
本体受理状態のc04follow_decisionは実closed枝、独立rowのc04follow_decisionsは既存last_flow分類呼出しの観測。
後者は元のpenultimate再構成に基づく分類であり、保存累積塑性すべりの差とは別である。
両辺がnull/unchangedなら判定を再評価していない旨をそのまま読む。欠損を0に置換しない。
各不一致レコードのactual_classification_operands.shortcut_conditionは本体の配列厳密一致の両辺、
inherited_decisionは保持元の受理座標と実closed判定の両辺。現在0増分と保持元の非0流れを混同しない。

measurements.jsonは絶対L2・生参照ノルム・分母・相対L2、全GP更新分割分布、P/R条件付き差と平均主効果を保存する。
*-signed-effects.jsonには符号付き全成分を保存する。ノルムを足して100%の原因割合にしない。
m*-native??-independent??.jsonは全32比較の不一致全件、工程→受理刻み→固定GP順である。
戻り後余裕scaledの1e-14帯と、保存塑性すべり増分の厳密0を別に確認する。
帯は観測だけで、カテゴリを一致へ置換したり新しい受入幅にしない。
A-04の力学32表は引用のみ。今回3水準の漸近次数/外挿を算出しない。

## 封印再読込とソフトウェア全件

上の新規求解とは別に、封印の全member SHAを確認し、保存場から観測ファイルを再現できる。
再集計は追加の全体求解に数えない。

```powershell
& $Py "$D/c04follow_evidence.py" --verify $D --extract work/NEW-c04follow-verified
if ($LASTEXITCODE -ne 0) { throw 'part/ZIP/memberのSHA不一致' }
& $Py "$D/c04follow_report.py" --run work/NEW-c04follow-verified/c04follow-formal --output work/NEW-c04follow-reported
if ($LASTEXITCODE -ne 0) { throw '観測の再現未完了' }
& $Py "$D/c04follow_validation.py" rust --output work/NEW-c04follow-full/rust
& $Py "$D/c04follow_validation.py" python --output work/NEW-c04follow-full/python
& $Py "$D/c04follow_validation.py" collect --output work/NEW-c04follow-full
& $Py interim/2026-12-06-p6-01-status-reconciliation/audit_status.py
if ($LASTEXITCODE -ne 0) { throw '固定出典/台帳整合不成立' }
```

fmt/clippy/build、Rust全ターゲット＋doc、tests/toolsの全発見Python7群を対象にする。
成功/失敗/error/未実行/skipをvalidation.jsonへ分離する。重い診断はCIに追加しない。
全件試験のPC予算30〜180分・2〜8GB、既存約29GB証拠の読取は小型求解のRAMとは別（Windows未実測）。
文書末尾はuniversal newlineと実CRLFで検査、ノルムは固定順fsum。Windows差を許容値で吸収しない。
PCへ回すのはこの上限内入力のWindows再現だけ。未定義仕様をPC計算で補完しない。
接合面水頭種/所有写像、界面透水・貯留と開閉、物理時刻・世代、連成方向/フィードバック、全場散逸/カテゴリ許容値はユーザー判断待ち。

push後のWindows CIは未確定。他のA/B/C、NC-01、P6も全体未達のまま。

保存済みLinux結果はRust467/Python929成功、失敗/error・未実行・skipは0。
正式16履歴と、全件/登録後確認/封印再集計の件数を混ぜない。
P6最終登録は350出典・SOURCE_NOT_INVENTORIED=0。Windows CIは当該push後未確定。

C-04は未達、passed=false。
