# 範囲拡張(a)：先端JOINT＋BIRTH／湿潤土BIRTH／Richards湿潤JOINT

PR #13 の `codex/plastic-spatial-convergence`、指定先端 `8d62e7b2636be67e3b6095383ca9f98f0424bad5` を起点とする。同じブランチへpushし、PRはマージしない。最終の限定受入は `acceptance.json`、全件チェックは `validation.json`、独立測定は `measurements.json` を正とする。

## 受理する範囲と残す拒否

| 組合せ | 今回の範囲 | 範囲外 |
|---|---|---|
| 1 先端宣言付きJOINT＋BIRTH | 乾燥・単一の適合Tri6面集合、全所有要素が活動するA/B、接合専用 `joint_death` と再誕生 | 湿潤、複数JOINT、構造部材との複合、先端を跨ぐボディBIRTH、部分活動面の再誕生 |
| 2 湿潤土BIRTH | 完全水没の新規MC土1ボディ、`cad_new_nodes`、明示した同一静水位、両側が土の単一CONTINUOUS JOINT・B | SOIL／CONTINUOUS／PRESCRIBEDの定常流を伴うBIRTH、部分水没、新設MCC、古い材料履歴、未定義初期水圧、活動中の土、他の配置・構造材・先端との複合 |
| 3 Richards湿潤JOINT | **受理しない** | ONE_WAY・双方向とも拒否。接合面の水頭写像、時刻・開閉・世代と過渡水量の整合、質量収支と水圧荷重の独立再組立、時間刻み細分化の結合検証が未実装・未検証 |

2で使うCONTINUOUSは静水圧の接合面方針であり、定常Darcy BIRTHの合格を意味しない。Richardsの単体機能・既存ONE_WAY機能は変更していない。本追補の新しいBIRTH入力は `stages::Stage::Birth.hydraulic_initialization` であり、既に別経路で検証されている `one_way_construction::Action::Birth`（明示水頭を持つRichards施工）は従来どおりである。Richards湿潤JOINTと土だけの既存Richards施工を区別する。結合の水量収支や刻み収束を測定していない経路を合格にしない。

## 実装と明示入力

先端節点の `−N_i + N_i = 0`、縮退要素、部分接合と未分離側の接着を維持する。先端は節点対 `[n,n]` の乗数未知数ではなく、一つの物理節点である。誕生前の `Transform::with_node_pairs` の推移的同値類と整合するよう、Aの伝達力回収から自己対を除く。共有先端は接触力0の行として残し、再現不能な共有先端の自由残差を絶対1e-10で拒否する。既存の面積重み付きGP投影・2000GP／6000節点成分の上限・全節点絶対1e-10 kN検査は維持する。先端付近の非一様初期力を平滑化・クリップしない。引張、摩擦円外、再現不能な投影を拒否する。これは既存の縮退Tet10近似の検証であり、連続体の特異応力点値やXFEMの高精度検証ではない。

Bは計算した誕生時ジャンプをそのまま保存する。先端近傍の1e-12 mの非ゼロジャンプを0へ置き換えない試験も含む。新世代のすべり・累積すべり・塑性散逸・数値散逸・誕生後仕事は厳密に0から開始し、旧世代を再積分しない。全所有要素が活動しない接合の再誕生と、先端を跨ぐボディBIRTHは拒否する。

湿潤土は例えば次を指定する。`states` は選択ボディと厳密に一致する新規MC材料状態であり、有効応力・塑性ひずみ・間隙比・履歴を省略できない。

```json
{
  "kind": "birth", "name": "fill", "bodies": ["UpperClay"],
  "placement": "cad_new_nodes",
  "hydraulic_initialization": {"kind": "hydrostatic", "level_m": 1.5},
  "states": {"UpperClay": {
    "effective_stress_compression_positive_kpa": [5,5,5,0,0,0],
    "plastic_strain_compression_positive": [0,0,0,0,0,0],
    "void_ratio": 0.7, "equivalent_plastic_strain": 0,
    "plastic_work_kpa": 0, "last_return": "elastic"
  }}
}
```

初期水圧は `p0(x)=γw max(level_m−z,0)`。既存の静水位と一致し、選択ボディ全体が水没することを検査する。材料履歴は新規入力、有効応力はこの入力の生値を誕生準備状態へ保存する。既存土の変形を新設土の初期ひずみにしない。古い塑性ひずみ・累積塑性ひずみ・塑性仕事・非弾性returnは拒否する。水圧のクリップや未定義値の自動補完で受理しない。

体積水圧は内部力、露出面・接合面水圧は物理右辺へ各1回だけ入る。活動領域が変わるBIRTHでは、既存 `joint_water::regenerate_birth` で面荷重を再積分し、露出／被覆面も活動要素から再組立する。`retain` の0化を新しい誕生面の生成に使わない。静水圧以外の初期化形式は型と実行時ガードで拒否する。

## 独立参照と保存再開

新規 `tools/joint_ext_reference.py` は入力だけからGL125体積剛性・Tri6接合・全体拘束のSVD零空間を構成する。先端A/Bの全節点変位・全体積GP応力・全接触GP力と、配置時変位を独立計算した湿潤新設土の後載荷を比較する。保存された材料応力を用いる全節点再組立は、独立弾性解と分けて記録する。既存 `verify_wet_structures.py::reference` と `joint_birth_reference.py` の許可範囲は拡張していない。

各受理経路は、誕生前・誕生途中・すべり後・撤去後・再誕生途中／終点から別CLIで再開する。JSONの符号付き0も保持し、`elapsed_seconds` だけを除くRun全体を完全一致させる。意図的棄却と同じ受理刻み列の棄却なし計算は、全状態・フレーム・遷移監査・受理Attemptを完全一致させる。全受理Runについても、棄却Attempt、意図的棄却指定・初回刻み指定とその入力指紋、計測時間だけを区別して完全一致を検査する。物理入力、初期平衡、掘削監査を除外しない。再開同士の比較は計測時間以外を除外しない。

チェックポイントは外側ハッシュをRustで再計算した後も、初期力・基準・世代・凍結世代・水圧・遷移監査の改変を拒否する。新しい湿潤BIRTH入力でも保存地点までの再計算による認証を有効にした。CLI検証は固定コピーした同一バイナリで行い、ビルド中の実行ファイル置換を避ける。

再組立1e-10、非線形釣合い1e-8、線形残差1e-10、独立弾性相対1e-7、拘束1e-10 m、接触仕事収支絶対1e-10を維持する。構造並進・回転1e-9、板生値の絶対1e-8／相対1e-7、水頭差1e-8 mの既存基準も変更していない。

## 証拠と未達

生入力・全Run・棄却なし経路・範囲拒否ログを `evidence.zip`、再開と再ハッシュ改変を `restart-*.zip`、旧版との比較を `legacy-evidence.zip` に保存する。全エントリーのSHAは `evidence-manifest.json` に記録する。不合格・途中の試行は `retained-failures/` に残し、成功した最終実行と混ぜない。全件試験・専用CLIの最終件数と終了コードは `validation.json`、受入フラグと測定値は `acceptance.json` と `measurements.json` を参照する。

案1・案2・案3の過去記録と全ての失敗証拠を保持する。低プレロードNC全上部撤去、湿潤JOINT塑性すべり全履歴の刻み収束、P6全体は未達である。却下済みNC調整は行っていない。未収束・時間切れ・キャンセル・未検証条件を合格にしない。大規模・実CAD・GUI/GPUは [ユーザーPC手順](REPRODUCE.ja.md) に分ける。Windows CIのpush後成否は未確定として報告する。

## 最終の小型測定

- 6ケース・36フレームを全節点再組立。最大相対残差 `3.4763740517771664e-12`、物理右辺の最大成分差 `2.1316282072803006e-14 kN`。
- 入力からの独立弾性解は3ケース・14フレーム。最大相対差は変位 `4.3578372245105025e-12`、全体積GP応力 `4.7202987112678604e-12`、接触GP力 `5.079876774267678e-12`。
- 湿潤BIRTHの誕生準備状態から、新領域の体積水圧・露出面水圧・接合面水圧と仮保持荷重を別組立。初回・再盛土の全節点保持荷重の最大成分差 `1.7763568394002505e-14 kN`。準備スナップショットに残る直前工程の荷重／面水圧を、新領域の右辺の代用にしない。
- 接触仕事収支の最大絶対差はGPで `5.117434254131581e-17 kPa m`、全体で `3.2263530105670698e-18 kN m`。
- 別CLI再開41地点、意図的棄却18回と棄却なしの同じ受理刻み列、外側ハッシュを再計算した18改変の拒否。局所エネルギー恒等式も保つよう改変した凍結世代は、受理経路の再計算不一致で拒否された。
- 指定先端の元ソース70ファイルをGitの内容とバイト照合して隔離ビルド。既存の乾燥JOINT、湿潤JOINT、先端JOINT、Mindlin Tri6板の4完成Runは、計測時間と実行版SHAだけを除き全体が完全一致した。

`metrics.json` と `independent.json` はこれらの実測値を記録する。`source-snapshot.json/.zip` の817ソースは最終全件チェック前後で照合する。新しい結合のRichards水量収支・時間刻み収束は測定しておらず、0誤差や合格に置き換えない。

追加の `audit_checkpoints.py` で、完成した親Runに属する41保存状態の水圧・保持／解放荷重・接触力・拘束・仕事収支を別組立した。途中状態の最大相対残差は `6.238730734961501e-12`、全節点右辺の最大成分差は `1.4210854715202004e-14 kN`。結果は `checkpoint-reassembly.json`。途中保存だけを完成解析の合格にしていない。

## Linux最終全件チェックと判定

修正版で `cargo fmt --all --check`、`cargo clippy --locked --all-targets -- -D warnings`、release全ターゲットビルド、`cargo test --locked --release` が成功した。Rustは全444件（68集計、失敗・無視・除外0）。Pythonはtests全490件とtools2件、計492件が成功し、失敗・エラー・スキップ0。最後の完全実行は `logs/complete-python-tests.log` と `logs/complete-python-tools-tests.log`。途中の不合格・未配置・部分実行はこの件数へ加算していない。

`acceptance.json` の `passed=true` は表の限定受入範囲に対する判定で、`requested_all_combinations_passed=false` を保持する。Richards湿潤JOINTは `passed=false`、P6全体も未達である。試験後も817ソースのSHAが一致する。Windows CIのpush後成否は未確定。
