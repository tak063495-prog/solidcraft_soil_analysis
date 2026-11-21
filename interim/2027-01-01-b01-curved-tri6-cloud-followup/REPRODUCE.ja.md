# B-01クラウド追補の小型再現

PR #13、`codex/plastic-spatial-convergence`、起点c6b1176以降。マージしない。
GPD WIN Max 2、使えるメモリ約48GBを前提。Windows PCでは未実行。
最大8 Tri6・25節点・150 DOFだけであり、768/6144要素・40GB超・1経路1日超の新計画はない。
上限超の確認が必要な点は**規模上限のため確認不可・未達**。元の基準/参照値を変えない。

## 取得と準備

既存作業と旧証拠を上書きしない。出力先は毎回未存在のworkパスを使う。
初回取得なら `git clone --filter=blob:none --no-checkout` で開始し、チェックアウト前にsparse規則を設定する。
初回全履歴転送は行わず、P6台帳の固定コミットと出典Git blobは必要なOIDだけfiltered取得する。
A-02〜A-04追補のgit-recoveryと、本フォルダvalidation/git-recoveryを参照。
全件試験に必要な既存C-01群も取得し、新しいevidence.zip.part*を除外しない。
既存CIと同じpostfold除外だけを保持する。欠落を空ファイル/skip/保存値変更で代用しない。

```powershell
git config --local core.autocrlf false
git sparse-checkout set --no-cone '/*' '!/interim/2026-11-10-nc-postfold-branches/evidence.tar.gz.part*'
git fetch --filter=blob:none --no-tags origin codex/plastic-spatial-convergence
git switch codex/plastic-spatial-convergence
git merge-base --is-ancestor c6b1176b1147b8f242ef979f6c4a6a70e6829579 HEAD
if ($LASTEXITCODE -ne 0) { throw '指定起点の子孫でない' }
$Py = '.venv/Scripts/python.exe'
if (-not (Test-Path $Py)) { py -3.13 -m venv .venv }
& $Py -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw '固定依存の準備失敗' }
$env:CARGO_BUILD_JOBS='2'
$env:RUST_TEST_THREADS='2'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
$D='interim/2027-01-01-b01-curved-tri6-cloud-followup'
cargo build --locked --release --all-targets
if ($LASTEXITCODE -ne 0) { throw '本体ビルド未完了' }
```

Python3.13、Rust1.99.1/rustfmt/clippy、PowerShellを使う。
既存全件は約29GBの生証拠を読む。ビルド/証拠再展開/取得キャッシュにもディスクの余裕を確保する。
これは新しい大規模FEM計算ではない。今回のクラウドは実ファイルとGit blobの二重保持を考慮して空き容量を監視した。
LinuxでHOMEが読取専用ならHOMEを変えず、XDG_CACHE_HOME/XDG_CONFIG_HOME/XDG_DATA_HOMEをwork配下の書込可能な場所へ指定する。
NC固有設定は変更しない。PowerShellがなければ公式SHA確認済みportable版をworkへ置く。

partial cloneでは現先端の実ファイルだけでは過去コミットの文書blobが不足する場合がある。
今回の初回全件のP6-03はこの理由で1 errorとなった。固定出典を書き換えず、次を全件前に行える。
取得後はbranch-sources.jsonのSHA256/Git blobと照合する。復旧はvalidation/git-recovery/p603-fixed-blobs-recovery.json参照。

```powershell
& $Py -c "import json,subprocess; m=json.load(open('interim/2026-12-07-p6-03-cloud-reconciliation/branch-sources.json',encoding='utf-8')); ids=sorted({v['git_blob_sha1'] for v in m['files']}); subprocess.run(['git','-c','fetch.negotiationAlgorithm=noop','-c','gc.auto=0','-c','maintenance.auto=false','fetch','--no-tags','--no-write-fetch-head','--filter=blob:none','--stdin','origin'],input=('\n'.join(ids)+'\n').encode('ascii'),check=True)"
if ($LASTEXITCODE -ne 0) { throw '固定文書blobの取得未完了。全履歴取得やskipで代用しない' }
```


## 同じ小型診断を新規実行

```powershell
& $Py "$D/b01follow_build.py" --output work/NEW-b01follow-build
if ($LASTEXITCODE -ne 0) { throw '診断呼出器のビルド未完了' }
& $Py "$D/b01follow_run.py" --probe work/NEW-b01follow-build/b01follow_probe.exe --output work/NEW-b01follow-formal
if ($LASTEXITCODE -ne 0) { throw '数値未完了。途中入力・生出力・終了コードを保持' }
& $Py "$D/b01follow_report.py" --run work/NEW-b01follow-formal --output work/NEW-b01follow-measurements.json --effects work/NEW-b01follow-formal/cross-effects.json
if ($LASTEXITCODE -ne 0) { throw '集計未完了' }
& $Py "$D/b01follow_precision.py" --formal work/NEW-b01follow-formal --measurements work/NEW-b01follow-measurements.json --output work/NEW-b01follow-precision
if ($LASTEXITCODE -ne 0) { throw '同じ保存行列の高精度観測未完了' }
# 精度観測をreport出力へ追加したmeasurementsを指定して、保存差の法線/接線分解も再現できる。
# b01follow_components.py --formal ... --measurements ... --output work/NEW-components.json
```

元の生産Kernelを直接呼ぶだけ。曲面ガードを外す手順はなく、新しいdirector/曲率/物理はない。
ビルドを除く保守的PC予算は1〜10分・メモリ1GB以内（Windows未実測、受入条件ではない）。
ビルドは10〜90分・2〜8GBを見込む。プローブ1呼出し上限300秒。上限を延ばして1経路1日超へ進めない。
数値が終わらなければ未完了として残し、拒否/物理破壊/合格へ変換しない。

plan.jsonの32パッチ×6場=192直接解、許可256/拒否16要求、同入力別PID一致を見る。
現行呼出器は回転前の形状からexact-patch付記を分類する。封印の初回出力では全体z=0の付記だったため、
回転平面28場のメタデータだけが異なる。初回ソース/生フラグは封印に保持し、数値は訂正していない。
assemblies/に全K、全固有値/零空間、実6剛体場、全変位/回転・力/偶力・反力を保存する。
recovery-responses.jsonは全要素6点のひずみ・N/M/Q・drilling・4エネルギー。
並進力kNと回転モーメントkN mを分離し、mixed_dof_*は代数観測と読む。
一般回転は列ベクトルにRx17→Ry(-31)→Rz29の固定順。全支持空間と荷重も共変に移す。
再現用行列は求解前plan.jsonのround-trip定数に固定し、Windowsのlibm差で別の入力行列を作らない。
中間節点だけの有限移動は拒否。歪みは頂点変更後の直線中点配置に限る。

Decimalは同じ丸め行列のRHS/直接解を50/80桁にするだけ。高精度値で通常値を置換して合格にしない。
最悪平面膜場の通常回転差約1.21e−12が約1.05e−15へ縮むことと、基底/K差・条件数を別々に読む。
S6の投影外割合約7.55e−8とK r正規化残差約2.57e−17は別量。
固有空間の高精度診断はクラウドで可能な残件であり、この手順に新しい求解として追加していない。
外挿・Richardson補正・新しい許容値は使わない。

## 封印からの再集計と全件検証

```powershell
& $Py "$D/b01follow_evidence.py" --verify $D --extract work/NEW-b01follow-reopened
if ($LASTEXITCODE -ne 0) { throw 'part/ZIP/entryのSHA不一致' }
& $Py "$D/b01follow_report.py" --run work/NEW-b01follow-reopened/b01follow-formal --output work/NEW-b01follow-reported.json --effects work/NEW-b01follow-reopened-effects.json
if ($LASTEXITCODE -ne 0) { throw '保存生値の集計未完了' }
& $Py "$D/b01follow_validation.py" rust --output work/NEW-b01follow-validation/rust
& $Py "$D/b01follow_validation.py" python --output work/NEW-b01follow-validation/python
& $Py "$D/b01follow_validation.py" collect --output work/NEW-b01follow-validation
& $Py interim/2026-12-06-p6-01-status-reconciliation/audit_status.py
if ($LASTEXITCODE -ne 0) { throw '固定出典または台帳不一致' }
```

measurements.jsonには別時点のprecision_accountingも含むので、再現時はreopened/b01follow-precision/summary.jsonを
report出力のprecision_accountingへ追加してから全バイトを比較する。この操作は新規物理求解ではない。
同じ環境の固定順で再現し、異なるOSの最終ビット差は生値と環境を別記する。新許容値で吸収しない。
新試験は軽い契約/固定順/実CRLF/封印だけ。rust-integration-1への追加は0。
全件は既存7 Python群（tests/tools）とRust全ターゲット＋docを実行する。
PC予算30〜180分・2〜8GB（Windows未実測）。成功/失敗/未実行/skipをvalidation.jsonへ正確に残す。
SOURCE_NOT_INVENTORIED=0の最終監査と固定コミット/SHA確認を行う。

Windows CIはpush後未確定。ソフトウェア試験成功をB-01の数値受入へ変えない。
曲面候補・独立曲面全体解は7項目の仕様/実装決定待ち。A群・他B/C・NC-01・P6も未達。

B-01は未達、passed=false。
