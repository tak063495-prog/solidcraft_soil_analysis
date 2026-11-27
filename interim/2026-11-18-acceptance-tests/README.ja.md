# TODOの受入試験を先に固定する

起点は PR #13 の `codex/plastic-spatial-convergence` の
`5f564621c7396c2a605edbfe5451cd800cf318d3`。追加先は
`codex/acceptance-tests-todo`、PRのbaseは `codex/plastic-spatial-convergence`。
既存ソース、既存試験、CI設定、`docs/P6_STATUS.ja.md`、既存interimフォルダーは変更しない。

対象の一覧は [todo.ja.md](../2026-11-16-issue-triage/todo.ja.md)。試験中の(a)は
[JOINT/BIRTHの基準](../2026-11-13-joint-ext-combos/README.ja.md)、(b)は
[Tri6/オフセットの基準](../2026-11-14-plate6-curved-offset/README.ja.md)を指す。
`tests/acceptance_todo/` の試験名のIDが下表と対応する。

## 対応表と判定の意味

|項目|試験名（同じ接頭辞のパラメーター展開を含む）|固定した基準・参照|xfailの理由／残る範囲|
|---|---|---|---|
|A-01|`test_A01_hydrostatic_all_node_force_and_work_reference` / `test_A01_wet_tip_birth_nodal_water_work_and_real_restart` / `test_A01_multiple_tip_birth_input_only_full_fields_and_restart`|静水圧 `p=γw(H-z)`、Tri6全節点荷重・仮想仕事、全節点再組立、接触仕事、A/B両BIRTH、並進した独立2組の162節点/全GP場、別CLI保存再開|(a)の先端BIRTHは単一・乾燥のみ。現在の範囲拒否をstrict xfailにする。相互作用/湿潤を含む複数JOINT・構造複合の独立3D全場は未整備|
|A-02|`test_A02_undefined_reactivation_remains_refused`|先端を跨ぐ土BIRTH、所有要素が全活動でないJOINT BIRTHの個別拒否|(a)が現行拒否を要求。成功時の投影・活動面仕様は未定義なので、成功を期待するxfailは創作しない。引張・摩擦円外・改変拒否の既存試験は維持|
|A-03|`test_A03_darcy_head_flux_and_birth_subdomain_reference` / `test_A03_supported_hydrostatic_birth_audit_positive_control` / `test_A03_darcy_birth_hydrostatic_manufactured_limit`|1D Darcyの線形水頭・一定流束、土BIRTH前後の静水圧極限、CONTINUOUS/SOIL/PRESCRIBED、全節点荷重、実再開|CONTINUOUS/PRESCRIBEDは定常流後BIRTHの明示水理初期化ガード、SOILは土/不透水構造と新設土を指定してBIRTH→Darcyを試す（誕生前は水理土領域が空）。非ゼロ勾配のBIRTH初期化形式は現行型にない。陽性対照は既に許可された静水圧BIRTHのみ|
|A-04|`test_A04_saturated_richards_transient_storage_balance_reference` / `test_A04_soil_only_saturated_limit_positive_control` / `test_A04_wet_joint_saturated_pressure_mapping_and_balance`|飽和Richardsの正弦過渡参照、水収支恒等式、一定全水頭の2水準、物理時刻到達、水圧写像・面荷重、元の流量/水頭ゲート|ONE_WAYの湿潤JOINTを拒否する現行ガード。土だけのONE_WAYを湿潤JOINTの成功へ転用しない。双方向・開閉・世代・非飽和過渡の結合収束は未整備|
|B-01|`test_B01_curved_tri6_metric_and_rigid_motion_reference` / `test_B01_planar_rigid_patch_positive_control` / `test_B01_curved_rigid_six_modes`|(b)の二次幾何・計量・幾何法線、剛体6モード。二次曲線辺・面外中間節点・円筒/球面の24パッチ。平面6モードを陽性対照|平面・直線以外を現行要素が拒否。曲面のガードが消えただけでは生値・回転・エネルギー検査を通過できない。余剰零モード、非剛体適合製造解、古典殻h系列は未整備|
|B-02|`test_B02_offset_strip_exact_reference` / `test_B02_offset_water_plate_fields_and_real_restart` / `test_B02_dry_full_field_reference_positive_control` / `test_B02_wet_zero_offset_legacy_limit_positive_control`|(b)の `κ=eN/D`、`εs=N/(Et)+eκ`、片持ち/単純支持・正負/ゼロe、正負e×初期水/工程水、入力だけのGL125全節点変位/全GP応力・設置参照・中立面回復、実再開|新方式と水の複合を現行ガードが拒否。JOINT・塑性MC/MCC/IDH・アンカー・土BIRTH/DEATH/MATSETを含む独立全場はこの試験群だけでは受入しない|
|B-03|`test_B03_outside_half_thickness_remains_refused` / `test_B03_inconsistent_shared_signed_normal_offset_remains_refused`|有限 `abs(e)<=t/2`、共有節点の符号付き法線オフセット。正負範囲外・段差・法線反転の個別拒否|現行仕様は拒否を維持する。(b)に段差/折れの成功仕様はなく、平均法線で受理する成功条件は作らない。Native/CAD/表示はPCで別検証|

参照式だけの試験、既存の許可範囲を使う陽性対照、現行拒否を維持する試験がpassしても、
各項目の完了を意味しない。曲面節点は円筒/球面上の点を通る二次補間であり、厳密な円弧・球面そのものではない。
古典殻ベンチマークの実測合格/不合格をこの幾何拒否から作らない。

新規71件の最終小型実行：

|項目|試験数|strict xfail|pass|
|---|---:|---:|---:|
|A-01|6|4|2|
|A-02|2|0|2|
|A-03|7|3|4|
|A-04|6|2|4|
|B-01|32|24|8|
|B-02|14|4|10|
|B-03|4|0|4|
|合計|71|37|34|

FAIL/ERROR/通常skip/XPASSはいずれも0。passには参照式の自己検査、陽性対照、拒否維持を含む。
パラメーター展開は各1件として数え、再実行を加算しない。

小型参照の具体値：DarcyはL=1 m、K=1e-6 m/s、下端全水頭2 m、上端2/2.5/3 mで
上向き流束0/−5e-7/−1e-6 m/s。板ストリップはE=2e8 kPa、t=.1 m、N=10 kN/m、L=2 m、ν=0。
e=+.05 mでκ=3e-5/m、εs=2e-6、片持ち先端(u,w,θy)=(4e-6 m,−6e-5 m,6e-5 rad)、
単純支持中央w=1.5e-5 m。負のeでは曲率/たわみ/回転の符号だけが反転し、e=0ではκ=w=θy=0、先端u=1e-6 m。

## ゲートとxfailの扱い

(a)/(b)の既存値をそのまま使う：再組立1e-10、非線形平衡1e-8、真の線形残差1e-10、
独立弾性相対1e-7、構造並進/回転絶対1e-9、板生値絶対1e-8/相対1e-7、
付着拘束1e-10 m、接触仕事絶対1e-10。試験が扱う量に該当するゲートのみ適用する。
B-02の新規 `elastic_offset.py` は乾燥正負eと既存方式の湿潤e=0極限で陽性対照を通す。
非ゼロe＋水の本体ガードは維持し、旧参照の許可範囲も変更しない。
新しい独立3D全場解を実装していない他の項目に、独立弾性相対1e-7の「合格」を記録しない。
Richardsの水頭1e-8 m・相対流量1e-8・絶対流量1e-16 m3/sは
`src/richards.rs::Settings::default` の既存定義。静水圧平衡極限では参照流量/蓄積が0なので絶対流量ゲートを用いる。
参照の過渡正弦解は飽和線形極限に限り、非飽和van Genuchten流れの参照へ転用しない。

未実装項目はすべて `pytest.mark.xfail(strict=True, raises=CapabilityUnavailable)`。
実際のCLIが文書に対応する既知の範囲拒否を返し、完成Runを出していない場合だけこの例外へ変換する。
参照値、履歴、ソースを改変してガードを迂回することはない。
未知の入力エラー、数値不収束、シグナル停止、120秒タイムアウト、欠落バイナリ、参照検査失敗は通常の失敗。
実装後に全後続アサーションを満たせばstrict XPASSで失敗し、実装をレビューしてから当該xfailを削除する。
`skip`、動的な `pytest.xfail()`、無条件の失敗/xfail、`--maxfail`による部分実行は使わない。

保存再開は実バイナリによる `checkpoint-stages` と別プロセスの `resume-stages` を使う。
計測時間 `elapsed_seconds` だけを除いた全Runをcanonical JSONで比較し、符号付き0も区別する。
保存JSONの再監査だけを実再開と呼ばない。

## 収集と検証

新規ディレクトリーには意図的に `__init__.py` を置かず、pytest関数で実装する。
既存CIは `tools/run_ci_python_tests.py` のunittest収集器でtests/toolsを収集し、4群へ分割する。
その収集器も `python -m unittest discover -s tests -p 'test_*.py' -v` もこのディレクトリーへ降りない。
既存のpytest設定、CIの件数ゲート、既存TestCaseの挙動を変更しない。
新規試験は既存設定のまま `python -m pytest tests/acceptance_todo -ra` で全件実行する。
ローカルconftestは新規試験用バイナリのパスと既存toolsのimportだけを用意する。
4群に含まれる既存toolsの2件も実行対象である。

追加前後の既存512 IDと全4群の割当は完全一致した（270/234/4/4）。
新規pytest収集・実行は71件、34 pass / 37 strict xfail / 通常skip・FAIL・ERROR・XPASS 0、終了コード0。
既存512件も4群すべて終了コード0、失敗・エラー・skip・expected failure・unexpected success 0。
既存tests 510件＋tools 2件を含み、NC mean/spatialの各4件もすべて実行した。
合計583件（546 pass / 37 xfail）。再実行や途中実行は重複加算しない。
全ログ・SHA・71件の項目別内訳は [validation.json](validation.json)、既存512 IDと群割当は [existing-inventory.json](existing-inventory.json)。
既存core-0/core-1/mean/spatialの実測時間はそれぞれ353.303/245.393/810.942/136.838秒、新規71件は6.86秒。
元の封印証拠の取得不足があった早期実行は最終成功件数に加えていない。追加ファイル以外の全既存Gitパスのmode/blob IDを起点と照合する。
NC-01、C-01〜C-04、P6の全行、過去の失敗・拒否、元の許容値と参照値は保持する。
新規試験成功やCI成功をこれらの完了へ転用しない。Windows CIの時間切れ・キャンセルは成功ではない。

## 試験化し切れていない基準

- A-01：相互作用/湿潤を含む複数JOINT・構造複合の全場独立解と全世代・改変checkpointの組合せ（独立2組の乾燥先端モデルは試験化済み）。
- A-02：先端を跨ぐBIRTHと部分活動面再誕生の成功仕様、特異点を平滑化しない初期投影の独立解。
- A-03：非ゼロDarcy勾配の新領域の水理初期化形式、部分水没・MCC/旧履歴の許可仕様。
- A-04：双方向の入力/連成方程式、開閉/世代の水頭写像・水量移送、非飽和過渡の独立解と刻み収束基準。
- B-01：曲面剛性の余剰零モード検査API、非剛体適合製造解、3古典殻の固定入力/出力量/参照値/許容差/3水準h系列。
- B-02：水以外の各複合と独立3D全場/改変拒否、塑性履歴/工程の設置基準。
- B-03：段差/折れ・範囲外e・コンパクト/Native/実CAD/表示を成功として受理する幾何/付着/表示仕様。

これらは基準不明と未実装の監査手段を区別して残す。全体の未達を埋めるための新許容値は定義しない。
大規模メッシュ・長系列・GUI/GPU/実CADはクラウドで実行しない。
[Windows手順](REPRODUCE.ja.md)に小型再現とPC上での追加証拠作成手順を記す。
