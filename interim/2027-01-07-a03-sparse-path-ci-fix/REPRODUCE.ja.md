# 疎選択のパスプロトコルの再現

起点6ee03451c49b86bc9f238ba0bed2f21b15f4299c、PR #13。マージしない。
新しいFEM求解・新しい大規模メッシュ・閾値変更はない。

必要なアーカイブと固定Git出典を全て配置したcheckoutで実行する。
取得はA-03 01-06のgit-recoveryと同じfiltered/固定OID方式を使い、初回全履歴転送を避ける。
旧封印を変更せず、新しいフォルダもP6台帳へ不変コミット・SHA付きで登録する。
Windowsではcheckout前にcore.autocrlf=falseとし、文書比較はuniversal newlineを使う。

```powershell
$Py='.venv/Scripts/python.exe'
& $Py -m unittest discover -s tests -p test_a03_unequal_thickness_anisotropic_k.py -v
```

上はWindows/Linux共通の12件の契約検査。次は本セッションで使用したLinux用の原因再現スクリプトで、
Windows形式の論理パスをLinuxのGitへ直接渡す。保存JSONのLinux再現ラベルも含めた実行ソースを維持している。

```bash
.venv/bin/python interim/2027-01-07-a03-sparse-path-ci-fix/a03path_reproduce.py --output work/NEW-a03path-reproduction
```

このLinux対照はWindows OSでの実行を代用しない。stdin/stdoutの生bytesと結果を新規保存先へ残す。
既存の保存先は上書きしない。NUL bytesはシェルのパイプへ流さずPythonからGitへ直接渡す。

全件は既存の検証器を再使用する。PowerShellを必要とする既存試験にはPowerShellを準備し、
読取専用HOMEではXDG_CACHE_HOME/CONFIG_HOME/DATA_HOMEをwork内にする。HOMEやNC設定は変更しない。

```powershell
$D='interim/2027-01-06-a03-unequal-thickness-anisotropic-k'
& $Py "$D/a03aniso_validation.py" rust --output work/NEW-a03path-full/rust
& $Py "$D/a03aniso_validation.py" python --output work/NEW-a03path-full/python
& $Py "$D/a03aniso_validation.py" collect --output work/NEW-a03path-full
```

collectの起点フィールドは旧検証器のラベルなので、今回の固定起点・実行ソースSHAは今回validation.jsonと全件証拠を参照する。
成功・失敗/error・未実行・skipと発見/実行IDを照合し、skipや許容値で差を吸収しない。
新しいWindows CIの結果は未確定。A-03を含む従来の全体未達を維持する。

A-03は未達、passed=false。
