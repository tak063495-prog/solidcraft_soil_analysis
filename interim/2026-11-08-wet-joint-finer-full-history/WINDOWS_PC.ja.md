# WindowsユーザーPCでの再現と残件

今回のクラウド計算はLinuxである。Windowsでの数値再現と、新しい先端のWindows CI成功は、この手順の記載だけでは確認済みにならない。追加クレジット購入は不要。PR #13はマージしない。

## ブランチと環境

未保存の作業がないことを確認し、リポジトリのルートで実行する。既存の作業を強制リセットしない。

```powershell
git fetch origin refs/heads/codex/plastic-spatial-convergence:refs/remotes/origin/codex/plastic-spatial-convergence
if ($LASTEXITCODE -ne 0) { throw 'fetch failed' }
git switch codex/plastic-spatial-convergence
if ($LASTEXITCODE -ne 0) { throw 'switch failed' }
git merge --ff-only origin/codex/plastic-spatial-convergence
if ($LASTEXITCODE -ne 0) { throw 'non fast-forward' }
git merge-base --is-ancestor c4de9fdafa4d7147f4f2bd8c8c877ca969b269b6 HEAD
if ($LASTEXITCODE -ne 0) { throw 'wrong history' }
git rev-parse HEAD
# Python 3.13と既存lockを使う。既存venvが適切なら作成は省略する。
py -3.13 -m venv .venv
.venv/Scripts/python.exe -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw 'dependencies failed' }
$Study = 'interim/2026-11-08-wet-joint-finer-full-history'
```

## 保存証拠の再計測（求解なし）

証拠は19 ZIP・約415 MiBである。全19ファイルをリポジトリから取得しておく。各 `NEW-*` は存在しない出力先を指定する。`restore.py` は全ZIPのSHA256を照合した後に展開し、全エントリーのSHA256も照合する。旧2048/32・1024/32のZIPはリポジトリ内の元の場所が必要である。

```powershell
.venv/Scripts/python.exe "$Study/restore.py" --output runs/NEW-wet-finer-evidence
if ($LASTEXITCODE -ne 0) { throw 'archive verification failed' }
.venv/Scripts/python.exe "$Study/study.py" --remeasure runs/NEW-wet-finer-evidence --output runs/NEW-wet-finer-remeasurement.json
if ($LASTEXITCODE -ne 0) { throw 'remeasurement failed' }
.venv/Scripts/python.exe "$Study/validate.py" --output runs/NEW-wet-finer-tests
if ($LASTEXITCODE -ne 0) { throw 'software checks failed' }
```

保存形式は、各場のlittle-endian float64のビット列をuint64と見て、時間方向の整数予測値 `2*previous-previous_previous (mod 2^64)` とのXORを取り、バイト順を並べ替えたNPZである。先頭は生ビット列、2番目は直前とのXOR。浮動小数の丸め・切捨て・閾値による省略はない。`study.decode` で全18場を復元できる。工程ごと・最大128受理点ごとに区切り、各チャンクの全受理メタデータとSHAを `accepted/manifest.json` に持つ。

再計測は324材料GP・32接触GPの全場、全受理点の進行度と残差、7工程端とのビット一致、既存の接触・エネルギー恒等式を検査する。2048/64と4096/64は共通進行点の全場差も計測する。細かい側だけの中間点を補間して比較したとはしない。

## 同じ小ケースの新規求解

```powershell
.venv/Scripts/python.exe "$Study/study.py" --output runs/NEW-wet-finer-solve 2>&1 | Tee-Object runs/NEW-wet-finer-solve.log
if ($LASTEXITCODE -ne 0) { throw 'diagnostic failed' }
```

既存3入力のみから、2048/64、4096/64の順に求解する。各経路の初期平衡を含む全7工程について材料・接触場を全受理点で保存する。過去の力学結果は新規2経路の求解完了後に比較のために読む。入力・力学解の生成へ過去の結果を使わない。計算時間やメモリー使用はPCによって異なる。終了コード0は診断完了を示し、`passed=false`、`history_converged=false`、`production_adopted=false` を変更する理由にはならない。

## 数値受入に残る作業

既存の全履歴誤差基準による証明はなく、観測差から新しい許容値を作らない。旧2048/32には工程途中の全場がないため、その接触内比較は7工程端の18場と受理列に限られる。今回の新規2経路では全場保存を実現したが、独立参照の全履歴受入や生産採用を意味しない。

生産側の解析水圧体積項・材料GP表示・水理遷移・掘削保持荷重・保存再開との整合は別途必要である。従来保存場を解析体積項で再組立した約9.55e-5の残差は前回の拒否証拠のまま保持する。今回の新規参照解でその記録を置換しない。

低プレロードNCは、既存 `make_free_contact.py --mcc --pc 100 --preload-kpa 5 --force 1 .5 0` の保存点0.9514000415802002・266増分から、別の検証済み改善根拠が得られた場合だけ検討する。今回は根拠を新たに得ていないので同条件の再試行を行っていない。Linuxの封印をWindows版へ付け替えず、力1e-8・真の線形1e-10等を保持する。100 kPaプレロード合格を転用しない。全上部撤去・全解放後凍結・独立塑性全場・両実装の刻み感度・完了までの再開は未達である。

Griffiths大規模系列、乾燥MC大規模系列、MCC274〜14909、斜面1200要素級SRM、仮想仕事上端再計算、82944降雨、105456性能、2807要素級実CAD、長時間Richards、実断面SRM本計算、Native UI、別PC/GPU、完成配布は今回クラウドで実行しない。必要性と予算をユーザーPC側で別途判断する。P3b/P4各行、構造・水・JOINT全体、P6全体は未完了のままである。
