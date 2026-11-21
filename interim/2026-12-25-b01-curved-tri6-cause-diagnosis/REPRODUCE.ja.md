# B-01の再現・ユーザーPCの小規模手順

GPD WIN Max 2、利用可能メモリ約48GBを想定する。以下はPCで実行していない。
各出力先は新規にし、旧証拠を上書きしない。許容値・支持・拘束は元のまま。
新たな768/6144要素や40GB超・1経路1日超のPC計算は計画しない。
上限を超える確認は**規模上限のため確認不可・未達**。許容値は緩めない。

## 1. 環境と封印の読込（求解なし）

リポジトリ直下のPowerShellで実行する。Python3.13、requirements-reference.lock.txt、
rust-toolchain.tomlのRust1.99.1を用いる。Git checkoutは旧封印のためautocrlf=falseを推奨する。
文書末尾を読む新しい試験は実CRLFでも一致する。既存の数値/封印検査は変更しない。

```powershell
$b01 = 'interim/2026-12-25-b01-curved-tri6-cause-diagnosis'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
python "$b01/b01cause_evidence.py" --verify $b01 --extract work/b01cause-pc-unsealed
python "$b01/b01cause_report.py" --formal work/b01cause-pc-unsealed/b01cause-formal --output work/b01cause-pc-remeasured.json
python -m unittest discover -s tests -p test_b01_curved_tri6_cause_diagnosis.py -v
```

想定1–2分、Python作業メモリ1GB未満、出力work/b01cause-pc-*。
見る値はpart/全entry SHA、保存行列からの物理剛体残差、9全体解の変形/反力、参照系列の次数null。
再集計は新規求解ではない。別OSで最終浮動小数ビットが異なれば差を別記し、許容値で吸収しない。

## 2. 同じ小さい内部診断を新規実行

```powershell
cargo build --locked --release
python "$b01/b01cause_build.py" --output work/b01cause-pc-probe
python "$b01/b01cause_driver.py" --probe work/b01cause-pc-probe/b01cause_probe.exe --output work/b01cause-pc-formal
python "$b01/b01cause_report.py" --formal work/b01cause-pc-formal --output work/b01cause-pc-observations.json
```

buildは未測定の目安10–30分/2–6GB、診断は1–2分/1GB未満。
最大8 Tri6、25構造節点、150DOF、全体直接解9個。新しい曲面殻の求解ではない。
元srcはバイト保持し、生成moduleはworkだけ。生成差分がガード1ブロック以外を変えるなら停止し、
新しい物理として続行しない。生産Kernelの同時ガード結果を必ず確認する。
平面/受理入力での剛性完全一致、曲面/直線対照の一致、実剛体回転の人工エネルギーを見る。

## 3. 旧21拒否＋受理対照4だけを標準CLIで再確認

```powershell
python "$b01/b01cause_cli_controls.py" --solver target/release/soil-fem-3d.exe --output work/b01cause-pc-cli-controls
```

想定1–5分/2GB未満。旧B-01の全part/entryをSHA確認して同一入力を取り出す。
見る値は拒否21、既存受理4、完成Run未生成の拒否、既存診断文字、受理Runの全比較。
異なる実行版間では認証済み実行版SHAと時間だけを別扱いする。
これを曲面ベンチマーク、改変拒否、曲面の実checkpoint再開へ数えない。
同一実行版の再開で時間だけを除く旧検査は既存全件試験で保持する。

## 4. 任意の平面多面体h系列（真の曲面の受入ではない）

旧B-01生成器のまま、屋根n=2/4/8＝8/32/128枚の平面直線Tri6に限る。
土は既存12 Tet10で屋根には未付着。点/拘束/drilling条件を旧入力から変更しない。
各出力の途中Run/不完了も残し、時間切れを幾何拒否・破壊・合格にしない。

```powershell
$oldB01 = 'interim/2026-12-03-tri6-curved-b01'
foreach ($n in @(2,4,8)) {
  $caseDir = "work/b01cause-pc-faceted-roof-n$n"
  python "$oldB01/b01_cases.py" $caseDir --kind roof --n $n --steps 4
  if ($LASTEXITCODE -ne 0) { throw '入力準備失敗。出力を保持する' }
  target/release/soil-fem-3d.exe solve-stages "$caseDir/mesh.json" "$caseDir/ground.json" "$caseDir/stages.json" "$caseDir/run.json"
  if ($LASTEXITCODE -ne 0) { throw '求解未完了。失敗を保持する' }
  python -c "import json; r=json.load(open('$caseDir/run.json')); f=next(f for f in r['frames'] if f['state'].get('stage_index')==2 and f['state'].get('stage_progress')==1); print('complete',r['complete'],'n',$n,'uz',f['state']['beams']['node_displacement_m'][$n][2],'force_residual',f['relative_force_residual'])"
}
```

想定各1–20分/1–8GB（PC未測定の保守的目安）。n=8でも約1,800構造DOF以下。
元の大規模土メッシュへは進めない。元n=2/4の保存値、符号付き自由辺中央uz、残差、
3水準差の向き/参考次数を見る。点変位の一致だけで全場収束を認定しない。
想定範囲を超えてメモリが増える場合は旧記録を残して中止し、同じ予算を拡大しない。

## 5. 真の曲面Tri6・標準殻のh系列

**実装・仕様決定待ち。現行の生産CLIには実行可能な曲面要素コマンドがない。**
したがって架空のsolverオプション、ガード解除した本番CLI、密行列上限解除は書かない。
同じ起点で上の手順3を実行すれば現行拒否だけを確認できる。

将来の候補は各対称領域でn=2/4/8＝8/32/128 Tri6の3水準を上限とし、
最初に単要素/8要素の物理剛体6モード、余剰モード、平面退化、適合曲面製造解と荷重仕事を確認する。
屋根は元R25/L50/t.25、E4.32e8、nu0、鉛直面荷重−90、元の自由辺中央uzを固定。
円筒は元R300/L600/t3、E3e6、nu.3、全荷重1/対称領域.25と元荷重点を固定。
開口半球は穴18度・R10/t.04、E6.825e7、nu.3、元荷重/支持/観測点を固定し、
穴なしやLE3へ変更しない。KL/MRと点荷重直下変位の意味は実装前に判断する必要がある。

目標のPC予算は各経路1時間以内/8GB以下、最長でも1日/40GBを超えない範囲だが、
**実装がないため時間/メモリは実測も確定見積りもできない**。実装後に8要素で測定し、
既存記録を残して次の水準の見込みが上限内の場合だけ手順を具体化する。
現段階の出力先は将来のwork/b01cause-pc-curved-{benchmark}-n{n}という計画名だけで、求解済みファイルはない。
必要なコマンドは実装・仕様決定後に固定し、今は「実装・仕様決定待ち」とする。
元の未定義許容値を便宜的に設定せず、元の判定を残す。
上限内の系列で確認できない連続極限/全場保証は、必要規模が上限を超えるなら
**規模上限のため確認不可・未達**。文献の値や外挿で合格を補わない。

## 6. ソフトウェア全件と台帳監査

既存分割を用い、全群を省略しない。下はPC未実行のコマンドで、Linuxの実績はvalidation.jsonを参照。
新しい診断試験だけを選んだ成功と、全件結果を混同しない。

```powershell
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release --all-targets
foreach ($g in @('rust-core','rust-integration-0','rust-integration-1')) {
  python tools/run_ci_rust_tests.py --group $g --output "work/b01cause-pc-validation/$g"
}
cargo test --locked --release --examples
foreach ($g in @('python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence','python-c01-switch','python-c01-path','python-c01-cause')) {
  python tools/run_ci_python_evidence_shards.py --group $g --output "work/b01cause-pc-validation/$g"
}
python interim/2026-12-06-p6-01-status-reconciliation/audit_status.py --help
```

想定全体1–4時間/2–12GB。1経路1日超の新しい数値検証ではなく、既存のソフトウェア全件である。
旧証拠約29GBのディスクと、テストの一時展開分の空きを確保する。
SOURCE_NOT_INVENTORIED、失敗/error、未実行、skipを別記し、許容値やNC設定を変えて通さない。
push後のWindows CIは今回のLinux結果とは別に未確定。
