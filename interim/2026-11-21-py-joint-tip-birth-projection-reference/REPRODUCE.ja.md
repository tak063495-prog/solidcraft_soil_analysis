# Windows PCでの再現

PowerShell 7、Python 3.13、リポジトリルートから実行する。既存参照値・許容値・テスト・CI設定は変更しない。パスは実PCの配置に読み替える。下記の時間はPC性能による未実測の目安。A-02の拡張範囲は全手順後も未達、`accepted=null`。

## 環境と小型参照（通常1分未満、依存導入1〜5分）

```powershell
git switch codex/py-joint-tip-birth-projection-reference
$d = 'interim/2026-11-21-py-joint-tip-birth-projection-reference'
py -3.13 -m venv work/joint-tip-venv
$py = (Resolve-Path 'work/joint-tip-venv/Scripts/python.exe').Path
& $py -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw 'dependency installation failed' }
$env:OPENBLAS_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$env:PYTHONUTF8 = '1'
& $py "$d/reference.py" --output work/joint-tip-pc/reference.json
if ($LASTEXITCODE -ne 0) { throw 'reference calculation failed' }
& $py "$d/checks.py" 2>&1 | Tee-Object work/joint-tip-pc/checks.log
if ($LASTEXITCODE -ne 0) { throw 'new checks failed' }
```

出力：`reference.json` のprojection、birth、numerical。合否：新規34試験が全件成功、skip/expected failure等0。これは参照実装の検証のみ。各部分活動参照の全節点力差、片側モーメント差、仮想仕事差、同入力/順序一致、引張/円外/再現不能の拒否を読む。既存1e-10 kNを投影再現に適用し、数値結果を見て新しい許容値を作らない。科学的受入はnull。

## PC用の長めの分離系列（目安1〜30秒、クラウド未実行）

```powershell
Measure-Command {
    & $py "$d/reference.py" --pc-series --output work/joint-tip-pc/series-16-128.json
    if ($LASTEXITCODE -ne 0) { throw 'PC series failed' }
}
```

GP系列は256/1024/4096/16384点、要素系列は16/32/64/128帯×固定GL8。出力は誤差、観測次数、弱形式の力残差と仕事欠損。粗/細の同じ定義を保つ。丸め域では次数が不安定になり得るのでnullや悪化を保存し、外挿補正しない。最低次数・本体受入誤差の許容値は未定義、`accepted=null`。独立3D地盤の要素系列・特異先端点値・ボディBIRTHの本体求解の代用にはしない。

## 人工CSV/JSONの往復と将来本体export比較（目安1〜5秒）

```powershell
& $py "$d/compare.py" --write-synthetic --output work/joint-tip-pc/synthetic.json
if ($LASTEXITCODE -ne 0) { throw 'fixture generation failed' }
& $py "$d/compare.py" --input work/joint-tip-pc/synthetic.json --output work/joint-tip-pc/synthetic-comparison.json
if ($LASTEXITCODE -ne 2) { throw 'expected undefined acceptance, exit 2' }
& $py "$d/compare.py" --write-synthetic --output work/joint-tip-pc/synthetic.csv
& $py "$d/compare.py" --input work/joint-tip-pc/synthetic.csv --output work/joint-tip-pc/csv-comparison.json
if ($LASTEXITCODE -ne 2) { throw 'expected undefined acceptance, exit 2' }
```

CSVには同名stemの `.meta.json` が必要。人工データの `source_kind` は`synthetic_reference`のまま保持する。これを本体比較済みの証拠にしない。

本体の将来exportがREADMEの固定乾燥ケースと列・GP・節点・フレーム・入力指紋まで一致する交換契約を実装した後だけ、次を実行する。**現行本体にこのexportコマンドはなく、架空のCLIオプションは記載しない。** 変換器が必要なら別承認の実装とし、原本の力を補間・クリップ・平滑化しない。`native_export`はラベルだけなので、原入力・実行版SHA・原本CSV/JSONのSHA・実行コマンドもPC上に保存する。

```powershell
& $py "$d/compare.py" --input work/joint-tip-pc/native-export.json --output work/joint-tip-pc/native-comparison.json
$status = $LASTEXITCODE
if ($status -eq 1) { throw 'export rejected: inspect errors' }
if ($status -ne 2) { throw 'unexpected completion status' }
```

正常診断は2、拒否は1。`accepted=null`、A-02未達。範囲外6入力（部分再誕生、先端跨ぎ、引張、円外、再現不能、水圧組合せ）の拒否記録が必要。欠落・逆行世代、旧世代改変、非ゼロの新世代すべり/仕事、同入力順序依存、全節点力欠損も拒否する。試験内の故意改変を無効化して通さない。

## 既存Python全件（PC目安5〜30分、保存証拠の展開空間を確保）

既存試験はunittest。4群のrunnerは既存設定をそのまま使い、収集件数条件を変更しない。保存証拠の再計測は重い新規求解とは異なる。全追跡ファイルを取得し、部分チェックアウトなら必要な既存archiveも元Git blobで復元する。

```powershell
cargo build --locked --release
if ($LASTEXITCODE -ne 0) { throw 'release CLI build failed' }
foreach ($g in @('python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence')) {
    & $py tools/run_ci_python_tests.py --group $g --output "work/joint-tip-pc/existing/$g"
    if ($LASTEXITCODE -ne 0) { throw "existing Python group failed: $g" }
}
```

合否は全発見IDの重複・欠落0、全群の実行数が割当数と一致し、fail/error/skip/expected failure/unexpected successがすべて0。時間切れ・中断・キャンセルは未完了で成功扱いにしない。新規34件は別起動で、既存件数に紛れ込ませない。Linux実測は `results/validation.json`。Windows本手順の未実行と、新先端Windows CIの結果は区別する。

今回のLinux最終実績は既存534件（265/261/4/4）＋新規34件＝568件、上記異常各0。群の実測は約170/386/809/126秒。PCでこれより速く終わることを保証しない。初回の未配置ZIPによる通常群の6エラーは保持し、元Git blobを復元した後に265件を全件再実行した。ファイルを見つけられないときはskipや新規参照値の生成で置換せず、指定起点の元データの配置とSHAを確認する。

大規模3D、実CAD、一般バルク塑性との結合、湿潤JOINT、水理、長時間BIRTH経路はこの参照の対象外。許可された入力・仕様・比較契約が揃うまで実行や合格判定を補完しない。
