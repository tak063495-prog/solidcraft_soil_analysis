# A-04クラウド追補の再現

A-04は未達、passed=false。PR #13へ積み、マージしない。
元の受入基準・拒否・既定を変更せず、土だけの施工と別の乾燥/静水圧JOINTだけを使う。
Windows手順はGPD WIN Max 2・利用可能メモリ約48GBの上限内。Windowsは今回未実測。
新しい768/6144要素、メモリ40GB超、1経路1日超の計算を追加しない。
上限を超える一般保証は「規模上限のため確認不可・未達」とする。

## 取得と環境

既存リポジトリのPR #13ブランチで行う。起点はd9ca6634a7c394d0188a502c11c9038a4a0b98e5。
初回取得は `--filter=blob:none --no-checkout` とし、先に通常CIのsparse規則を適用する。
通常CIでも除外している `interim/2026-11-10-nc-postfold-branches/evidence.tar.gz.part*` を一括先読みしない。
新しい本フォルダーのevidence.zip.part*は除外しない。全件試験では既存のC-01各群を含む必要な旧partも揃える。
不足したアーカイブを空ファイル・生成値・skipで代用しない。

A-02/A-03追補のgit-recoveryと、今回validation/git-recoveryの実行ログを参照する。
必要な固定Git blob IDだけをfiltered取得し、Git blob/SHA/サイズを照合する。
今回の一括checkoutは未index化packを中断し、完全prefixだけをstrict indexで回収した。
原header・未完suffix・SHAを保存し、旧の作業ファイルやコミットを書き換えていない。
その後の必要blobは個別取得した。初回の全履歴転送は行っていない。

Python 3.13とrequirementsの既存固定版、Rust 1.99.1/rustfmt/clippyを使う。
PowerShellがないLinuxでは公式SHA照合済みportable 7.5.4をwork/へ置ける。
NC固有の環境設定を変更せず、PowerShell不在をskipで通さない。 LinuxのHOMEが読取専用の場合はHOMEを変更せず、
`XDG_CACHE_HOME` / `XDG_CONFIG_HOME` / `XDG_DATA_HOME` をそれぞれwork/配下の書込可能な
cache/config/dataへ指定してからPowerShell試験を実行する。これらは一般のランタイム格納先で、NC固有設定ではない。
証拠再展開用の一時空き容量も必要。今回の初回全件は空き不足で失敗し、旧作業ファイルのGit blob一致を確認した上で
再取得可能な単一blob promisorキャッシュだけを整理して約19.3GBを確保した。旧のファイル/参照値/コミットは保持した。

```powershell
$Py = '.venv/Scripts/python.exe'
$D = 'interim/2026-12-31-a04-richards-wet-joint-cloud-followup'
$env:RUST_TEST_THREADS='2'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
cargo build --locked --release --all-targets
if ($LASTEXITCODE -ne 0) { throw '本体ビルド未完了。ログ保持' }
& $Py "$D/a04follow_build.py" --output work/NEW-a04follow-build
if ($LASTEXITCODE -ne 0) { throw '隔離診断ビルド未完了。ログ保持' }
```

本体ソースはコピーして読むだけ。交換とQ逆順は専用ビルド内にしか存在しない。
通常CLIへ新しい環境変数・物理・要素・既定を追加しない。
各経路の実行版とソースはbindingでSHA固定し、指紋を付け替えない。
ビルドの保守的予算は10〜90分・メモリ2〜8GB（Windows未実測）。

## 3点の新規全体求解

```powershell
$env:A04FOLLOW_PROCESS_SECONDS='3600'
try {
  & $Py "$D/a04follow_run.py" --build work/NEW-a04follow-build --output runs/NEW-a04follow
  if ($LASTEXITCODE -ne 0) { throw '診断未完了。途中状態・PID・失敗ログを保持' }
  & $Py "$D/a04follow_surface.py" --formal runs/NEW-a04follow --output runs/NEW-a04follow-full-surface
  if ($LASTEXITCODE -ne 0) { throw '専用静水圧ray則の順序対照が未完了' }
  & $Py "$D/a04follow_input_audit.py" --run runs/NEW-a04follow --output runs/NEW-a04follow-input-audit.json
  if ($LASTEXITCODE -ne 0) { throw '物理入力の同一性が確認できない' }
  & $Py "$D/a04follow_report.py" --run runs/NEW-a04follow --surface runs/NEW-a04follow-full-surface --output runs/NEW-a04follow-measurements.json
  if ($LASTEXITCODE -ne 0) { throw '集計未完了。求解完了へ読み替えない' }
} finally { Remove-Item Env:A04FOLLOW_PROCESS_SECONDS -ErrorAction SilentlyContinue }
```

土・順序・交換を別プロセス群にする場合は `--part soil` / `--part order` / `--part exchange` を使い、
それぞれ別の新規outputへ保存する。同じフォルダーへの上書き再開をしない。
今回の実行は土・順序・交換を別にソース固定し、正式集計用フォルダーへバイト同一でまとめた。
各元プロセスのログとsource-snapshot/bindingを保持している。

- 土は鉛直2/4/8層＝12/24/48 Tet4、dt=.25/.125/.0625秒、独立Jacobi3/5/7の27経路。
  本体既定積分3で全9格子も解く。共通16時刻を補間なしで比較する。
  3D等方細分化の率ではなくh_z限定。積分の次数を変えたGPを最近接対応にしない。
- 土VG、既存乾燥先端B、完全水没静水圧CONTINUOUS/BのE/F/Q全8組合せを両実装で解く。
  Eでは明示crack-tip owner番号も同じ物理要素へ付け替える。接続・material・regionの組を崩さない。
  Qでは点/重み/形状関数の組を逆順にする。追加4経路は専用ray則の内外点も逆順。
- C-04の12 Tet4・水位32/力学32で各実装にP/Rの全4組合せを載せる8経路。
  全体平衡を入力から解き、相手の解・分割・受理列を渡さない。
  特殊な中立/冷間予測子が必要な経路は交換未完了として止め、未使用分岐まで移植したとはしない。

全行列のPC予算は30〜240分・メモリ1〜8GB、各プロセス上限1時間（Windows未実測）。
予算内で完了しなければ数値未完了・未達を保持する。経路予算を1日超へ延ばす新手順はない。
土の独立積分キャッシュは同一水頭引数だけに使い、元の式・許容値を変えない。
元の本体力/線形基準、独立残差、厳密Coulomb/no-tension、元の端点散逸を維持する。

|出力|見る値|
|---|---|
|plan.json、source-binding.json、source-snapshot/、isolated-build-binding.json|求解前の行列・実行版とソースの固定SHA|
|mesh/ground/case.json、*-process.json、*.log|生入力・PID・終了コード・時間切れ。拒否と数値未完了を分ける|
|土native.json、*-fields.jsonl、*-layout.json|P0セルH/面H、各水理GP、χp内力・水重量、全水量台帳と全受理場|
|orderのnative/run.json、accepted/constitutive/progress.jsonl|全材料/接触履歴、反力、Newton停止列、世代。GPを付け替えない|
|exchangeのfields.jsonl、result.json、constitutive.jsonl|全416状態、散逸/すべり、厳密カテゴリ、実分割数、独立の再分割イベント|
|measurements、order-input-audit、order-point-set-audit|条件付き主効果、混合有限差、各因子4対、同じ物理集合の証明|

ノルムと節点仕事の観測集計は固定順math.fsum。反力の順序差を新許容値で吸収しない。
今回、乾燥open/.25の1置換が元の力残差で1反復早く受理され、反力差が元ゲートを超えた。
PCで同じ点を見る場合も失敗を残し、残差基準の変更で受入を作らない。
水理を持たないJOINTの流量・貯留・水収支はnullであり、土の台帳を流用しない。

## 封印の読取と全件検証

```powershell
& $Py "$D/a04follow_evidence.py" --verify $D --extract runs/NEW-a04follow-reopened
if ($LASTEXITCODE -ne 0) { throw 'part/結合ZIP/entryのSHA不一致' }
& $Py "$D/a04follow_report.py" --run runs/NEW-a04follow-reopened/formal --surface runs/NEW-a04follow-reopened/retained/a04follow-full-surface --output runs/NEW-a04follow-reopened-measurements.json
if ($LASTEXITCODE -ne 0) { throw '保存場の再集計失敗' }
& $Py "$D/a04follow_validation.py" rust --output runs/NEW-a04follow-validation/rust
& $Py "$D/a04follow_validation.py" python --output runs/NEW-a04follow-validation/python
& $Py "$D/a04follow_validation.py" collect --output runs/NEW-a04follow-validation
& $Py interim/2026-12-06-p6-01-status-reconciliation/audit_status.py
if ($LASTEXITCODE -ne 0) { throw 'P6固定出典/台帳の不一致' }
```

全件はRust全ターゲット＋doc、Pythonは既存7群（testsとtools）を全て実行する。
新試験は軽い契約・CRLF・固定順・封印SHAだけで、重い数値再計算をCIへ追加しない。
全件のPC予算は30〜180分・メモリ2〜8GB（Windows未実測）。実測の成功/失敗/未実行/skipをvalidation.jsonへ記録する。
単独スモーク・登録後再確認は全件の固有IDへ重複加算しない。
文書はuniversal newlineで読み、実CRLFファイルでも末尾一致を検査する。
封印読取・再集計は新規求解ではない。外挿・Richardson補正・新定義だけで合格を作らない。

この手順は同じ小型入力のOS間確認で、S1〜S7/S12の仕様決定や双方向連成実装ではない。
当該push後のWindows CIは未確定。A-01〜A-03、B、C、NC-01、P6も全体未達を保持する。

A-04は未達、passed=false。
