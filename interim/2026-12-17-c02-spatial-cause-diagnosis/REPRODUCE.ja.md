# C-02 原因切り分けの再現・PC手順

2026-12-17は依頼されたフォルダ名。実施日は2026-11-14（日本時間）。PR #13 の起点は `026a4302fc5ed9b65eead8a493b1e11850b8b9ee`。マージしない。以下の `NEW-...` は未存在の出力先へ置換し、非0終了・途中停止・欠落を保存する。既存の保存値・封印・checkpointのSHAを付け替えない。PC欄の大規模計算はクラウドでは実行していない。

## 環境・取得

PowerShell、Python 3.13、Rust 1.99.1。ディスクは全件試験用に80〜120 GB、PCの細空間・全受理保存には500 GB以上を計画する。時間・資源のPC値は計画値で、Windows実測ではない。

```powershell
git config core.autocrlf false
git fetch --filter=blob:none origin codex/plastic-spatial-convergence
git switch codex/plastic-spatial-convergence
git pull --ff-only origin codex/plastic-spatial-convergence
git merge-base --is-ancestor 026a4302fc5ed9b65eead8a493b1e11850b8b9ee HEAD
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
py -3.13 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements-reference.lock.txt
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:CARGO_BUILD_JOBS='2'
$env:RUST_TEST_THREADS='2'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
$D='interim/2026-12-17-c02-spatial-cause-diagnosis'
$Old='interim/2026-11-19-wet-joint-spatial-convergence'
```

Linux/WSLではvenvと実行版名から`.exe`を外す。クラウド環境のactivationは `/workspace/soil-environment/activate.sh`。途中ファイルは `work/` に置く。

Linuxの全件試験にも既存のPowerShell試験が含まれる。今回と前回C-01で使ったportable 7.5.4を `work/powershell-7.5.4/` へ配置し、PATHに追加した。公式packageは [PowerShell 7.5.4 linux-x64](https://github.com/PowerShell/PowerShell/releases/download/v7.5.4/powershell-7.5.4-linux-x64.tar.gz)、SHA-256は `1fd7983fe56ca9e6233f126925edb24bf6b6b33e356b69996d925c4db94e2fef`。`XDG_CONFIG_HOME`・`XDG_DATA_HOME`・`XDG_CACHE_HOME`・`DOTNET_CLI_HOME`・`DOTNET_BUNDLE_EXTRACT_BASE_DIR` はwork配下へ向け、HOMEは変更しない。NC固有の環境変数・設定は変更しない。PowerShell未配置による最初の通常群のskipは正式成功にせず、ログを保持して再実行した。これらは既存の小さいPowerShell単体試験で、大規模PC計算の実行ではない。

今回の再実行で、既存の `test_main_thread_runtime_sigterm_reaps_real_child` がworker-pidの存在だけを待ち、内容書込み前の空文字を読む競合が一度発生した。失敗した全群ログを `validation/retained-core-pid-race/` に保持し、該当単独試験と全core0を逐次で再実行した。fixture・アサーション・実装を変えずに通ったが、競合を修正したとは扱わない。再現時にもログと全件結果を分けて記録する。

通常の全件試験は既存CIの7 Python群・3 Rust群の分割を維持する。C-01専用群だけがZIPを絞り込む既存CI設定は変更していない。通常core群には新しい `evidence.zip.part*` が必要。全件試験前には少なくとも、各試験が読む旧アーカイブと `interim/2026-12-05-nc-independent-arc-residual/reference400-fixture.zip` の存在・SHAを確認する。ZIP除外を使ったcheckoutの後にfixtureを取り忘れない。未使用の `interim/2026-11-10-nc-postfold-branches/evidence.tar.gz.part*` だけは従来どおり除外できる。これは試験のskipではない。

この作業では巨大な履歴一括fetchとbulk lazy fetchを中断し、固定commitのtreeから必要ファイルを選び、raw配信の生バイトをGit blob SHA-1と照合して取得した。旧アーカイブは改変しない。取得結果と固定台帳の17 commitの取得を `validation/` に記録した。通信・取得未完の失敗も残した。

## 軽い検査・封印証拠

```powershell
python -m unittest discover -s tests -p test_c02_spatial_cause_diagnosis.py -v
python "$D/c02diag_evidence.py"
python "$D/c02diag_driver.py" --mode manufactured --output work/NEW-c02-manufactured
```

通常CIはこの軽い数学的試験・封印認証を行い、全体平衡8経路や20有限再生を再起動しない。`c02diag_evidence.py` は全part、ZIP全体、全エントリー、診断ソースを検証する。ドキュメント末尾は`read_text`でCRLFを正規化し、封印ZIPとソースSHAは生バイトで照合する。

封印ZIPの各エントリー名とSHAは `archive-index.json`。ZIPを展開する場合は、認証を先に行い、partを番号順に連結して新規フォルダへ展開する。`global/` は新規全体平衡の全場・結果・配置・コマンド、`measure/` は同一保存場の点差、`labels/` は開区間と望遠和、`replay/` は有限局所再生、`variant-measure/` は離散化感度、`provenance/` は旧封印の選択来歴、`source/` は新規診断ソースである。

## 同一保存場からの再計測

```powershell
$Saved='work/NEW-c02-saved-fields'
python "$D/c02diag_extract.py" --output $Saved
if ($LASTEXITCODE -ne 0) { throw 'old seal extraction failed' }
python "$D/c02diag_driver.py" --mode measure --source $Saved --output work/NEW-c02-measure
python "$D/c02diag_label_study.py" --source $Saved --output work/NEW-c02-labels
python "$D/c02_label_mask_robustness_check.py" --labels work/NEW-c02-labels --output work/NEW-c02-label-mask-check
python "$D/c02diag_replay.py" --source $Saved --output work/NEW-c02-replay
```

旧86part・ZIP全体と選択465エントリーを検証する。大きい旧JSONLを再生成せず、接触履歴に必要なNPZだけを新規展開する。選択は約1.61 GB。再構成は全672物理点、全7工程、共通実在時刻のみ。ノルムと内積は固定C順・`math.fsum`。重い処理の計画時間は5〜20分、RAM 4〜8 GB、追加ディスク3〜8 GB。

見る値は `measurements.json` の点/面積L2・p・方向余弦・補間の法則違反、`label-budget.json` の符号付き投影と交差項、`local-budgets.json` の有限16→32差と3成分恒等式。開区間の違いはカテゴリ端点の違いと別。局所再生を全体平衡・厳密解・誤差上界として扱わない。旧測定との表示一致と、旧Linuxバイトの完全再現は別に記録し、Windows差を新許容値で吸収しない。

## 接触離散化の新規全体平衡

```powershell
$Variants='work/NEW-c02-variants'
foreach ($variant in @('legendre4-tri6','jacobi4-corner3')) {
  $rule=if ($variant -eq 'legendre4-tri6') { 'legendre4' } else { 'jacobi4' }
  $interp=if ($variant -eq 'legendre4-tri6') { 'tri6' } else { 'corner3' }
  foreach ($h in 0,1) {
    foreach ($g in 32,64) {
      python "$D/c02diag_driver.py" --mode solve-variant --source $Saved --output "$Variants/$variant/h$h-g$g" --h $h --steps $g --rule $rule --interpolation $interp
      if ($LASTEXITCODE -ne 0) { throw 'input-only global equilibrium failed' }
    }
  }
}
python "$D/c02diag_variant_measure.py" --source $Saved --variants $Variants --output work/NEW-c02-variant-measure
```

12/96 Tet10、32/128接触GP、元の独立残差1e-10、接触最大5e-5 mを維持する。8経路のクラウド実測は個別に封印し、全て完了。Windows計画10〜30分、RAM 4〜8 GB、ディスク1〜3 GB。見る値は各 `result.json` の全7N受理列・最大残差、`layout.json` の点・重み・補間、`discretization-summary.json` の空間差の変化と時間との相互作用。角節点補間は改善案ではなく感度対照。mortarを本体へ導入しない。

## PCのみ：独立768 dense系・最細の時間差

既存ドライバで**入力から**768を独立に解く。旧保存履歴を初期値へ転写しない。まず64/128、資源と終了コードを見て256、必要なら512へ。計画RAM 32〜64 GB、空き100〜200 GB、各経路1〜12時間以上（未実測、機種依存）。

```powershell
cargo build --locked --release --example wet_joint_spatial_probe
$FineTime='work/NEW-PC-c02-768-time'
python "$Old/study_c02.py" --mode init --max-level 2 --user-pc --solver target/release/examples/wet_joint_spatial_probe.exe --output $FineTime
foreach ($g in 64,128,256) {
  foreach ($h in 0,1,2) {
    python "$Old/study_c02.py" --mode native --h $h --steps $g --user-pc --timeout 86400 --output $FineTime
    if ($LASTEXITCODE -ne 0) { throw 'native incomplete' }
    python "$Old/study_c02.py" --mode independent --h $h --steps $g --user-pc --output $FineTime
    if ($LASTEXITCODE -ne 0) { throw 'independent dense incomplete' }
  }
}
python "$Old/study_c02.py" --mode restart --h 2 --steps 256 --user-pc --timeout 86400 --output $FineTime
python "$Old/measure_c02.py" $FineTime work/NEW-PC-c02-original-measure --max-level 2
foreach ($g in 64,128,256) {
  python "$D/c02diag_pc_measure.py" --source $FineTime --output "work/NEW-PC-c02-fixed-order-g$g" --max-level 2 --steps $g --user-pc
}
```

旧 `measure_c02.py` は元の受入定義をそのまま評価する保存済みのドライバ。新 `c02diag_pc_measure.py` は同じ全場・共通座標について固定順のノルムを別に記録し、結果を`passed=false`で出す監査である。旧JSONを上書きしない。PC監査CLIは今回クラウドで実行していない。

見る値：本体/独立の同一768・同一Nの全節点変位、全材料GP、全接触GP、全7工程の差、N64→128→256の全量の差と方向、接触GP散逸と塑性すべり、残差・受理列・実checkpoint再開一致。元の力1e-8・真の線形1e-10・独立1e-10を維持し、Dの局所分割差を隠さない。独立同一メッシュ一致だけでは連続解を認定しない。N256を時間収束済みと決めず、差が再増加したら512の別経路と失敗を残す。

## PCのみ：6144 Tet10 / 2048接触GP

さらに細い4水準を別出力へ入力生成する。計画RAM 64〜128 GB、空き500 GB以上、1〜数日。独立denseの6144は既存の15000自由度guardを越えるため、このコマンドでは起動しない。guardを外す仕様変更を隠して行わず、独立最細は前節の768で確認する。

```powershell
$Space='work/NEW-PC-c02-6144'
python "$Old/study_c02.py" --mode init --max-level 3 --user-pc --solver target/release/examples/wet_joint_spatial_probe.exe --output $Space
foreach ($g in 64,256) {
  foreach ($h in 0,1,2,3) {
    python "$Old/study_c02.py" --mode native --h $h --steps $g --user-pc --timeout 172800 --output $Space
    if ($LASTEXITCODE -ne 0) { throw 'spatial path incomplete; retain failure and checkpoint' }
  }
  python "$D/c02diag_pc_measure.py" --source $Space --output "work/NEW-PC-c02-h0123-g$g" --max-level 3 --steps $g --user-pc
}
python "$Old/study_c02.py" --mode restart --h 3 --steps 256 --user-pc --timeout 172800 --output $Space
python "$Old/measure_c02.py" $Space work/NEW-PC-c02-original-h0123 --max-level 3
```

見る値：12/96/768と96/768/6144の固定Nの全場L2・p・差方向、N64/256で空間差がどれだけ変わるか、開閉区間と点別寄与、全工程欠落なし、元の残差と別プロセス再開一致。反転・非単調・有限参照の不確かさが残れば未達。外挿やRichardson補正を使わず、新しい許容値・新しいL2定義だけで合格を作らない。

## 全件検証とP6-01登録

```powershell
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release
cargo test --locked --release --all-targets
cargo test --locked --release --doc
foreach ($group in @('python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence','python-c01-switch','python-c01-path','python-c01-cause')) {
  python tools/run_ci_python_evidence_shards.py --group $group --output "work/NEW-c02-tests/$group"
  if ($LASTEXITCODE -ne 0) { throw 'full Python group failed' }
}
python interim/2026-12-06-p6-01-status-reconciliation/audit_status.py
```

試験を始める前に必要な旧part・fixtureと固定commitのtreeがあることを確認する。台帳監査はlazy fetchを使わないので、`inventory.toml`の全固定commitと元BASEを取得してから行う。全発見IDと実行IDの和集合、重複・未実行・skip・expected failureを数え、準備中の再試行を正式成功へ二重計上しない。今回の新規JSONは科学内容をcommit固定後、既存 `register_sources.py --commit FULL_SHA --id UNIQUE_ID interim/.../acceptance.json` 等で出力したTOMLを追記し、`audit_status.py`の補助台帳SHAだけを更新する。旧登録を再生成しない。

C-02は未達、passed=false。C-01、C-03、C-04、NC-01、P6、A/B各項も未達。push後の新先端Windows CIは未確定。
