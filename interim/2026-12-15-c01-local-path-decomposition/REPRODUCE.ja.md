# C-01 局所経路検証・寄与分解の再現

PR #13、`codex/plastic-spatial-convergence`、起点 `73e29c327745fbc045b84b23b405408a67e91fe8`。マージしない。
以下はリポジトリルートの PowerShell。`NEW` は未存在の名前に置き換える。PC実行は今回行っていない。
フォルダ名の2026-12-15は指定された作業識別子である。

## 準備

生SHAと固定ソースを検査する再現は、次の新規cloneを推奨する（既存の作業変更を上書きしない）。
`core.autocrlf=false`はこのcloneだけに設定し、JSON・ソース・証拠のバイトを保持する。
既にCRLFへ変換されたソースを持つcheckoutでは、変更を保持したまま別cloneを用意する。文書末尾の比較自体はCRLFにも対応する。

```powershell
git clone --single-branch --branch codex/plastic-spatial-convergence --config core.autocrlf=false https://github.com/tak063495-prog/solidcraft_soil_analysis.git NEW-C01
Set-Location NEW-C01
```

Git履歴と作業treeだけで計50 GB規模になり得る。以下の出力・一時領域の予算はこれと別で、更新履歴により増減する。
リポジトリルートから続行する。

```powershell
git fetch origin codex/plastic-spatial-convergence
git switch codex/plastic-spatial-convergence
git pull --ff-only origin codex/plastic-spatial-convergence
git merge-base --is-ancestor 73e29c327745fbc045b84b23b405408a67e91fe8 HEAD
# sparse checkout中なら全必要証拠を確保する。過去のZIP/part/失敗記録も必要。
git sparse-checkout disable
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
cargo build --locked --release --example wet_joint_local_path_probe
```

Linux/WSLは `.venv/bin/python`、実行版に`.exe`を付けない。文書比較だけは`read_text`で改行を正規化する。
証拠・入力・JSON・ZIP/part・ソース指紋のSHAは生バイトで検証し、改行変換しない。
Git履歴と全アーカイブの取得自体に数十分かかる場合がある。新しい診断コードのimport名は`c01_`で始め、旧NC診断の`diagnose`等を使わない。

## 4つの解析可能な局所経路（全体求解なし）

```powershell
python interim/2026-12-15-c01-local-path-decomposition/c01_local_study.py --solver target/release/examples/wet_joint_local_path_probe.exe --output work/NEW-c01-local
```

入力は4経路、各線形区間の外側分割1/2/4/8/16/32、内部固定1/2/4/8、端点下限、prescribed trial列による凍結拡張規則再生の計168経路。
`request.json`、本体の`response.json`、独立解析解との全増分比較`local-results.json`を出す。
想定数秒〜1分、RAM1 GB以下、ディスク数十MB。Windows実測ではない。
見る値は塑性すべり・traction・増分／累積散逸の誤差、符号、観測次数、非単調性、端点下限と実分割の差、非負散逸、仕事収支。
凍結規則の局所再生はNewton試算を入力として与える診断であり、全体平衡を新しく解く操作ではない。

## 保存済み全GP履歴からの寄与分解（全体求解なし）

```powershell
python interim/2026-12-15-c01-local-path-decomposition/c01_decomposition.py --prepare --inputs work/NEW-c01-saved-inputs --output work/NEW-c01-decomposition --refinement 32
python interim/2026-12-15-c01-local-path-decomposition/c01_refinement_audit.py --report work/NEW-c01-decomposition/decomposition.json --points work/NEW-c01-decomposition/point-increments.jsonl --output work/NEW-c01-decomposition/finite-reference.json
```

12-11/12-12の全part・ZIP・選択エントリーのSHAを確認して保存済み受理場を展開し、12-13の全partとZIPを検証して実分割数を取り出す。
同じ保存jumpを独立実装で実分割、端点下限、実分割×16/×32で再積分する。本体・独立の16/32/64/128/256、2制御の全20経路と、既に12-13で保存済みの切替制御・独立N512を加えた全21経路・全7工程・全32GPを使う。N512の新規求解は行わない。
`decomposition.json`は元の全履歴L2、符号付き寄与・交差項・差方向の分解、工程・対応増分の集計、再生参照の細分差。
`point-increments.jsonl`は工程・進行度・GP順に対応した実分割／下限・旧新法線圧・水圧・元散逸増分・局所再生増分。
`finite-reference.json`は再生×16/×32の分解の変化と欠損が大きい対応増分。`provenance.json`は旧索引・既存partition報告のSHA。保存場の再計算を新規全体平衡解と呼ばない。

想定10〜60分、RAM8〜16 GB、ディスク6〜12 GB（取得済みGit履歴を除く）。Windows実測ではない。
見る値は実分割再生と保存場の一致、細分16→32の差、元の散逸非単調性・すべり差方向余弦、局所・変位経路・開閉区間の符号付き投影と交差項。
投影比は負や1超にもなる。二乗ノルムを足すだけの寄与率に置き換えない。
局所再生の参照細分差を調べるPC候補は`--refinement 64`で別出力へ再計算し、32/64の寄与が安定するかを見る（20〜120分、同程度RAM・追加数GB）。
外挿・Richardson補正は使わない。

## PCのみ：独立512/1024と本体・独立の元の全場系列

以下は前回の未完・未収束を引き継ぐ手順。上限・許容値を変えない。クラウドでもPCでも今回新たには実行していない。

```powershell
python interim/2026-12-13-c01-cause-diagnosis/c01_cause_diagnostics.py --extract switch --output work/NEW-PC-c01-switch-saved
python interim/2026-12-13-c01-cause-diagnosis/c01_cause_diagnostics.py --extract path --output work/NEW-PC-c01-path-saved
foreach ($n in 512,1024) {
  python interim/2026-12-13-c01-cause-diagnosis/c01_cause_experiment.py --control switch --steps $n --contact-order 4 --pc-series --output work/NEW-PC-c01-reference/switch-n$n-q4
}
python interim/2026-12-13-c01-cause-diagnosis/c01_cause_compare.py --saved-switch work/NEW-PC-c01-switch-saved/paths --saved-path work/NEW-PC-c01-path-saved/paths --experiments work/NEW-PC-c01-reference --reference-levels 256 512 1024 --output work/NEW-PC-c01-reference-comparison
python interim/2026-12-13-c01-cause-diagnosis/c01_cause_experiment.py --control path --steps 512 --contact-order 4 --output work/NEW-PC-c01-reference/path-n512-q4
cargo build --locked --release --example wet_joint_path_probe
python interim/2026-12-12-c01-switch-event-localization/study.py --mode init --solver target/release/examples/wet_joint_path_probe.exe --output work/NEW-PC-c01-full
python interim/2026-12-12-c01-switch-event-localization/study.py --mode native --pc-series --steps 256 512 1024 --timeout 14400 --output work/NEW-PC-c01-full
python interim/2026-12-12-c01-switch-event-localization/study.py --mode independent --pc-series --steps 256 512 1024 --output work/NEW-PC-c01-full
python interim/2026-12-12-c01-switch-event-localization/measure.py work/NEW-PC-c01-full work/NEW-PC-c01-full-measure --steps 256 512 1024
```

独立切替512/1024は想定5〜60分以上、RAM16 GB、5〜15 GB。誤差制御512は既存2000試行/工程上限に達した失敗の再現でもある。
`failure.json`・途中NPZ・manifest・全ログを保存し、完了系列に混ぜない。
全場系列は想定1〜6時間以上、RAM16 GB、30〜60 GB。見る値は全18場の元の全履歴L2・絶対差方向・次数、元の力／線形残差、全工程完了。
保存済み512を収束済み参照や厳密解として使わない。1024も既存上限内で完了した場合だけ比較する。
PCで全体平衡を比較して初めて、局所再生では分離できない平衡との連成を検討できる。

## 全件検査・証拠・台帳

```powershell
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release
cargo test --locked --release
foreach ($group in 'python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence','python-c01-switch','python-c01-path','python-c01-cause') {
  python tools/run_ci_python_evidence_shards.py --group $group --output work/NEW-c01-validation/$group
  if ($LASTEXITCODE -ne 0) { throw "failed: $group" }
}
python interim/2026-12-15-c01-local-path-decomposition/c01_local_evidence.py
```

全発見IDと実行IDの重複なし和集合を確認し、成功／失敗／未実行／skipを`validation.json`に記録する。
新しい試験は軽量の通常Python群に配分され、過去のC-01の3専用シャードと75分予算を保持する。通常群で新ZIPが必要であり、全必要アーカイブを取得する。
旧NCコード／設定・旧シャード実装`run_ci_python_tests.py`は変更しない。
固定コミットへ保存した状態JSONを`interim/2026-12-10-ci-p601-inventory-fix/register_sources.py --commit FULL_SHA --id UNIQUE_ID interim/2026-12-15-c01-local-path-decomposition/FILE.json`で登録し、既存台帳を再生成せず追記する。
commit／SHA-256／Git blob／全boolean-null pointer検証を維持する。SOURCE_NOT_INVENTORIEDを例外・skipで逃がさない。
push後の新先端Windows CIは未確定。

C-01は未達、passed=false。
