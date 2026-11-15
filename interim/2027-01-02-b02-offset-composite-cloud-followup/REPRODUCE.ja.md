# B-02クラウド追補の小型再現

PR #13、`codex/plastic-spatial-convergence`、起点23c7b0ea48906723620b47eefb50e7481074d875以降。マージしない。
GPD WIN Max 2、使えるメモリ約48GB。Windows PCでは未実行。
最大96 Tet10・225土節点・8 Tri6＋梁1本。768/6144・40GB超・1経路1日超の新計画はない。
上限超が必要な一般確認は**規模上限のため確認不可・未達**。基準・参照値は変更しない。

## 取得と準備

初回は`git clone --filter=blob:none --no-checkout`、対象ブランチだけfiltered/depth=1で取得する。
必要なP6固定コミット/出典blobはOID単位でfiltered取得し、初回全履歴転送を避ける。
小型入力の計算と全件試験の既存約29GBの証拠読取は別である。
A-02〜B-01追補のgit-recoveryと今回validation/git-recoveryを参照する。
作業ディスクにはビルド・展開・Gitキャッシュの余裕を持たせる。キャッシュ削減は元Git blob/実バイトの再照合後だけ。
既存CIのpostfold除外以外の必要アーカイブ、新旧evidence.zip.part*を欠落させない。
Git worktree設定がある場合は`git config --show-origin --get-regexp 'sparse|partial|promisor'`も確認する。
今回、config.worktreeのsparse=falseがlocal設定を上書きしていたため、一括checkoutを中断して限定取得へ変更した。
失敗した取得を合格扱いせず、最終的に全必要実ファイルのGit blob一致を確認する。

```powershell
git config --local core.autocrlf false
git sparse-checkout set --no-cone '/*' '!/interim/2026-11-10-nc-postfold-branches/evidence.tar.gz.part*'
git fetch --filter=blob:none --no-tags --depth=1 origin codex/plastic-spatial-convergence
git switch codex/plastic-spatial-convergence
$Py='.venv/Scripts/python.exe'
if (-not (Test-Path $Py)) { py -3.13 -m venv .venv }
& $Py -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw '固定依存の準備未完了' }
$env:CARGO_BUILD_JOBS='2'
$env:RUST_TEST_THREADS='2'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
$D='interim/2027-01-02-b02-offset-composite-cloud-followup'
```

Python3.13、Rust1.99.1/rustfmt/clippy、PowerShellを用いる。HOME/NC固有設定は変えない。
読取専用HOMEではXDG_CACHE_HOME/XDG_CONFIG_HOME/XDG_DATA_HOMEだけをwork内の書込可能な場所へ指定する。
partial cloneのP6-03過去文書blob不足は、B-01追補REPRODUCEの固定blob filtered取得で先に防ぐ。
P6-01の固定コミットtreeも必要分だけ取得し、台帳のSHA256/Git blobと照合する。skipで代用しない。

## 同じ感度対照を新規実行

出力は毎回未存在のworkディレクトリを使い、過去の生値や失敗を上書きしない。

```powershell
& $Py interim/2026-12-01-tri6-offset-composite-b02/build_observer.py work/NEW-b02follow-build
if ($LASTEXITCODE -ne 0) { throw '公開API観測器のビルド未完了' }
& $Py "$D/b02follow_run.py" --observer target/release/b02-observer.exe --output work/NEW-b02follow-formal
if ($LASTEXITCODE -ne 0) { throw '診断未完了。入力・生出力・終了コードを保持' }
& $Py "$D/b02follow_report.py" --run work/NEW-b02follow-formal --output work/NEW-b02follow-measurements.json
if ($LASTEXITCODE -ne 0) { throw '観測集計未完了' }
```

18本体＋54入力のみの独立履歴。板厚t=.05/.1/.2 m、e=.025/.05/.1 mでe/t=.5を保持する。
元のe固定では薄板側が許可域を出るため、腕長/板自重も変わる系列と明記する。板厚だけの因果とはしない。
要素/Tet10/材料/領域・面・板の順序を同時に既定/逆順/固定seed置換する。
要素内接続と物理設置順は変えない。order-map.jsonの逆写像とcanonical入力SHAを確認する。
GL4/5/6は独立Tri6だけ。土・梁・本体の積分/ソルバー/停止条件を変えない。
新しい低減積分要素・MITC処方・接続・物理を導入する手順はない。

保守的PC予算は全系列10〜120分・メモリ1〜4GB、ビルド10〜90分・2〜8GB（Windows未実測）。
1子プロセス上限7200秒を延ばして大規模化しない。未完了は未完了として保存する。
**今回の薄板96要素は再組立1e−10を超えた未達の再現対象である。** 場ゲートの成功と分け、残差を消す設定へ変更しない。
各GLの全生場、observations.jsonの共通点・仕事・残差、assembly-gl*.npzの通常/固定順行列を読む。
符号付き全主効果/二因子/三因子差はsigned-factorial-effects.json。
反力/残差の並進力kNと構造回転モーメントkN mは分離し、混在量の相対値は代数観測として読む。
瞬間F·dと累積端点仕事は別量。固定順集計のWindows差を新許容値で吸収しない。

## 封印と全件検証

```powershell
& $Py "$D/b02follow_evidence.py" --verify $D --extract work/NEW-b02follow-reopened
if ($LASTEXITCODE -ne 0) { throw 'part/ZIP/entry SHA不一致' }
& $Py "$D/b02follow_report.py" --run work/NEW-b02follow-reopened/b02follow-formal --output work/NEW-b02follow-reported.json
if ($LASTEXITCODE -ne 0) { throw '保存場の観測再現未完了' }
& $Py "$D/b02follow_validation.py" rust --output work/NEW-b02follow-full/rust
& $Py "$D/b02follow_validation.py" python --output work/NEW-b02follow-full/python
& $Py "$D/b02follow_validation.py" collect --output work/NEW-b02follow-full
& $Py interim/2026-12-06-p6-01-status-reconciliation/audit_status.py
if ($LASTEXITCODE -ne 0) { throw '固定出典/台帳不一致' }
```

Rust全ターゲット＋doc、tests/toolsのPython全発見7群を実行し、成功/失敗/未実行/skipをvalidation.jsonへ記録する。
全件のPC予算30〜180分・2〜8GB（Windows未実測、受入条件ではない）。
新試験は契約/実CRLF/SHAだけ、重い再求解をCIへ追加しない。
同じ環境の再集計を比較し、OS間の最終ビットは環境と生差を別記する。
B-02の包括基準、拒否先複合、仕様判断待ちをソフトウェア成功で埋めない。
768残差は5e3d191の引用のみ。2空間水準の漸近次数・外挿・Richardson補正は計算しない。

push後のWindows CIは未確定。A-01〜A-04、B-01/B-03、C-01〜C-04、NC-01、P6も未達のまま。

B-02は未達、passed=false。
