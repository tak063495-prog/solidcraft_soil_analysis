# NC-01・案1：PC fold 再計算結果の取り込みと確定条件

**結論は「未確定・局所foldを強く支持」。PC結果によって支持は強まったが、refinement／独立実装の全場照合が不足する。NC-01は未達、P6その他の未達も維持する。**

PR #13、`codex/plastic-spatial-convergence` の指定先端 `e6e4e735495951661bcc7b1d72168c14866597b6` からの追記。PRはマージしない。フォルダー日付2026-11-30は依頼の保存先ラベルであり、転記対象のPC実施日は **2026-11-12 07:57–08:01 JST**。b069ec8の[案1の結論](../2026-11-17-nc-fold-stability/README.ja.md)を確定済みへ変更しない。数値コード・物理・ソルバー・既定・許容値・旧参照・旧hash/snapshotを変更していない。

## 出典と証拠の限界

[pc-results-transcribed.json](pc-results-transcribed.json) は**ユーザーがPC出力をチャットへ書き写した値**を保存したもの。クラウドの新規求解ではない。Windows、worktree `D:\Work\scsa_b069ec8`、HEAD `702a19a7ce189f1173caed727ee2a0059298a607`、Python 3.13.14。[PC再実行手順](../2026-11-24-nc-pc-rerun-fix/REPRODUCE.ja.md)の節1–2に基づく。元ZIPはクラウドから見えず、原JSON・manifest・全場・バイナリ/runtime指紋を直接検証していない。概数はJSONの`approx`として保持し、未転記値を補間しない。

両系列の固定入力検査はrc=0、`raw_equal=false`、`fixed_inputs_match=true`。既知の保存スキーマ差を除いた値の一致という**PC報告**である。元の入力・履歴を書き換えたという意味ではない。

共通コマンドは `tools/nc_fold_path.py --map native --fresh-native runs/NEW-nc-fold-PC/series/<128|256>/history --binary <observer>/material_diagnostic.exe --step-size <ds> --count <N> --budget-seconds 1800 --output <FixRun>/<id>`、続いて `tools/analyze_nc_fold_path.py --path ... --output ...-audit`。`<observer>`／`<FixRun>`は未転記の実パスを表す。4経路の計算・解析はいずれもrc=0。

## PCの数値と支持が強まった点

**128／256はメッシュ分割数ではなく、同一メッシュでの撤去初期刻み1/128・1/256の分母**である。`prepare_nc_fold_series.py` は `stages[2].initial_increment` だけを変える。固定54節点・12 Tet10・324材料GP・32接触GP、撤去後活動162材料GPという既存案1の対象を維持する。転記にある「メッシュ」の呼称は出典注記に残すが、この4本を空間収束の証拠には数えない。

|PC経路|ds m／点数|時間 s|最大サンプルλ|転記のbracket点|最終λ|最小σ(K)概数|
|---|---:|---:|---:|---|---:|---:|
|native128-coarse|2e-5／128|14.1|0.9514000963742456|7–8|0.9513899099079405|3.4242e-4|
|native128-fine|1e-5／256|28.4|0.9514000963737252|16–17|0.9513899099022921|3.4253e-4|
|native256-coarse|2e-5／128|14.4|0.951482232752456|14–15|0.9514733053795776|3.9585e-4|
|native256-fine|1e-5／256|28.3|0.9514822327804049|30–31|0.9514733053740204|2.8559e-4|

bracket点番号は転記表記のままであり、原JSONのstep／0始まりindexとの対応は未検証。bracket両端のλ、転記された接線・横断性・特異値はJSONに保持した。

全4本で、全点accepted、向きを継承した`dλ/ds`の正→負が1回、`verdict.conclusion=local_fold_bracket`。拡張系の最小特異値は約0.99999999（単位平面行の尺度）、rank(K)=54／rank(拡張系)=55、bracketの相対力残差は約5e-16、撤去済み全GP・接触は凍結を保持したとの報告。横断性の範囲は128粗0.372026–0.374364、128細0.372003–0.374351、256粗0.371489–0.373817、256細0.371466–0.373804。

転記小数の差を十進演算すると、局所ds半減のpeak差は128で **5.204e-13**、256で **2.79489e-11**。比例初期刻み系列128↔256の差は粗 **8.21363782104e-5**、細 **8.21364066797e-5**。同じ前史内の局所ds感度が小さく、両初期刻み系列で同じ局所折返しの形が得られた点は支持を強める。一方、局所ds差より系列間差が大きい。これを無視して全前史収束・独立一致・空間収束を宣言しない。値は最大**サンプル**であり、厳密な特異fold点そのものの受理ではない。

## 既存判定基準の点検

詳細な出典・実装上の区別・転記証拠の限界は[判定表JSON](criteria-assessment.json)を参照。b069ec8の案1README／REPRODUCE、702a19aの再実行README／REPRODUCEと、現在の`nc_fold_path.py`／`analyze_nc_fold_path.py`／`fold_stability.py`／`compare_nc_fold_reports.py`／既存`continue_path.py`を確認した。702a19aから指定起点e6e4e735まで、対象の数値・判定コードの差はない。既存基準は変更していない。

|対象|既存の要件／実装|今回のPC結果の判定|
|---|---|---|
|固定入力・元の共通λ|既知スキーマ既定値だけ補完した完全一致。全成分で`Ffinal+(1−λ)H`|固定入力一致の報告。analyze rc=0。ただし原入力・全荷重・前史の直接照合なし|
|局所残差|judgeは力≤1e-8、拡張系／力ブロック各≤1e-10。経路受理は力<1e-8、弧長拘束<1e-12 m|`local_fold_bracket`／acceptedは既存検査通過を報告。bracket力≈5e-16。全点の線形残差最大値・実行列は未転記|
|局所fold bracket|3点以上、受理連続枝、向き付き接線の正→負、横断性>64eps、近傍rank(拡張系)=n+1／rank(K)≥n−1|4本とも報告上充足。剛性のrank落ちを直接受理した意味ではない|
|荷重・接線・凍結の監査|全成分比例荷重、GP接線の組立一致、撤去済みGP／全接触の履歴不変、checkpoint/hash監査|凍結true、analyze成功の報告。原manifest／全場／GP接線のクラウド再監査なし|
|左右弱モード・特異値分離・低応力分布|実非対称Kの左右特異方程式、次の特異値との分離、全GPのp′・体積弾性率・位置・相関を保存|最小σとrank以外は未転記。PC原low-stress／モードを未監査。旧b069ec8の分布を今回のPC場へ転用しない|
|局所ds細分化|半減でpeak・モード・全GP／接触・荷重・滑らかさを比較|peak差は小さく支持が強まる。全場比較は未取得。**一致のpeak許容値は基準未定義**|
|比例初期刻み1/128↔1/256|元の入力・前史・比例経路を保って別比較|peak差≈8.21e-5。**一致基準未定義**のため一致合格にも不合格にも認定しない|
|メッシュ一致|この2系列は同一メッシュ。FreshNativeに幾何変更を通すオプションはない|真のメッシュ比較は未実施。メッシュ間peak一致基準も未定義|
|独立実装一致|入力由来前史400/800、局所ds半減、材料内4/8、本体との全場・モード・荷重比較|今回のPC転記はnativeだけ。旧本体対独立peak差≈1.78546e-4を保持。新しい独立一致は未達|
|新規診断restart|同じソース／binary／runtime／ds／入力から発行したcheckpoint＋sidecar、直接／分割の全保存場・実系比較|今回4本の再開照合は未報告。旧解析例・短経路の一致で代用しない|
|`static_limit_confirmed`|`follow()`で固定false。refinementからtrueへ昇格する実装なし|falseを維持。全撤去／本番再開を局所foldの必須条件と混同しない|
|`passed`／NC-01|診断report・analyze・比較器はいずれも固定false。生産受入器ではない|`passed=false`／`nc01_complete=false`を維持。転記の`fold_passed=false`も保持（この2ツールの標準fieldではなく集計元未確認）|
|全撤去・全解放後凍結・本番再開|NC-01には元条件の全撤去と独立全場・刻み感度・実checkpointで完了までの受入証拠が必要|全撤去false、全解放凍結false、完了再開false。**NC-01未達**|

`compare_nc_fold_reports.py` は数値差を表にするだけで、結論は未確定、peakの新しい一致許容値を作らない。したがって8.21e-5を既存の力／線形残差の許容値と比較することもしない。別の量の基準を流用しない。

全4本の `report.passed`、`history_converged`、`static_limit_confirmed`、`production_adopted`、`full_removal_reached`、`full_release_freeze_verified`、`restart_to_completion_verified` はfalse、`native_committed_increments=0`、弧長合計0.00256 m、停止は **bounded step count; inconclusive**。点数上限を物理的な到達不能・時間切れ・安定性の証明にしない。局所foldが確定しても別の非接続静的解の不存在や動的／熱力学的安定性を主張しない。

## 確定に不足するものと追加PC計算

不足するのは、まずPCの原JSON・SHA・全場・実系の持ち帰りと直接監査、次に同じ前史での追加ds細分化の全場照合、比例初期刻み系列の前史感度の評価、入力由来の独立400/800・材料内4/8と本体の比較、新しい実診断checkpointでの直接／分割一致である。数値の一致許容値が未定義の項目は、計算を増やすだけで自動合格にはならない。既存の判断範囲内で差と収束の状況をレビューする必要があり、今回新しい基準は設けない。

[REPRODUCE.ja.md](REPRODUCE.ja.md)に、固定worktree／HEAD、未存在の出力先、Windowsコマンド、想定時間、出力ファイル、各監査の既存合否条件を記載した。局所dsを2e-5→1e-5→5e-6 m、同じ弧長0.00256 mで比較し、同じdsの直接512点対256＋256点を照合する。独立前史は各5〜120分、native局所追加は各数十秒〜約1分が目安。PCからpushしない。真の細かいメッシュは現行FreshNativeの固定入力範囲外なので、この手順へ新設しない。初期刻み1/512等もFreshNativeの既存許可範囲外であり、新しいオプションは追加しない。

## Linuxでの再計算と全件検証

Linuxでは、SHA検証済みの旧1/128前史だけを既存`--evidence` adapterで読み、新しい材料観測器で粗32点／細64点を求解した。前史は再計算していない。両方とも弧長0.00064 mでbracket近傍を含む短区間であり、PCの128点／256点・弧長0.00256 mの全区間再現ではない。PCの1/256前史も再計算していない。PC側は別ファイルの`--fresh-native`履歴で、元ZIPがないため、以下のpeak差を純粋なOS差や全場一致の認定とはしない。PC転記JSONは上書きしていない。

|新規Linux経路|点数／ds m|経路実測 s|最大サンプルλ|対応するPC転記peakとの差（絶対値）|判定|
|---|---:|---:|---:|---:|---|
|linux-native128-coarse|32／2e-05|2.627|0.9514000963742423|3.3E-15|local_fold_bracket|
|linux-native128-fine|64／1e-05|5.184|0.9514000963737215|3.7E-15|local_fold_bracket|

2経路の計算・解析はそれぞれrc=0。力残差最大7.81076e-16、拡張系残差最大2.42539e-15、力ブロック残差最大1.39333e-15。いずれも元基準内で、全受理点の撤去済みGP／接触の凍結を確認した。両経路の`passed`／`static_limit_confirmed`／`full_removal_reached`はfalse、停止は点数上限による未確定である。これは既存の固定前史で局所bracketを再現した支持であり、新しい独立実装一致ではない。

[実測・PC差・コマンド・runtime・source binding](linux-observations.json)、[新しい全場・材料要求/応答・実行列・GP接線・checkpoint・sidecarと観測器の原出力ZIP](linux-fold-path-evidence.zip)、[ZIP全メンバーのmanifest](linux-evidence-manifest.json)を保存した。アーカイブ取得・全旧メンバーのSHA検証と観測器ビルドの時間は、表の局所経路時間に加算していない。[旧アーカイブの検証](observations/old-fold-extraction.json)、[取得ファイルのGit blob/SHA照合](observations/artifact-acquisition.json)、[ソースと同一メッシュの確認](source-review.json)も別記録とした。新しいLinux checkpointをPCのcheckpointとして扱わない。

**Python全件成功は未達。** 最終4群は573 IDを欠落・重複なく割り当て、skip／expected failure／unexpected successは0。[validation.json](validation.json)と`validation/`に全ログ・全ID・終了コード・現行ソース577件の指紋を保存した。

|Linux検証|最終実数・結果|
|---|---|
|fmt / Clippy locked all-targets warnings-as-errors / release build|終了コード 0 / 0 / 0|
|Rust locked release all-targets|454件成功、0失敗、ignored 0、filtered 0、87集計群|
|Rust doc|0件、終了コード0|
|Python通常群|300＋265＝565件成功、失敗・error・skip 0|
|Python旧NC mean|本体0件、準備エラー1、割当4件は未実行|
|Python旧NC spatial|本体0件、準備エラー1、割当4件は未実行|
|Python全発見集合|573件、実行本体565件、成功565件、本体未実行8件、準備エラー計2、skip 0|

旧meanは`ValueError: saved mean-stress/complement partition differs from own history`、旧spatialは`ValueError: remeasured JSON differs`で元の厳密検査が拒否した。Linuxにfloat差許可オプションを付けていない。元の参照値・許容値・封印source・hash/snapshot・拒否試験は維持した。各4件の本体を実行済みや成功へ数えず、過去の8件未実行記録も保持する。今回の失敗をNCの物理やfoldの新しい失敗判定へ結び付けない。

通常core-0の初回は単一ZIP証拠2件の取得漏れにより、296件実行・error 2・本体4件未実行だった。[初回ログとID](validation/retained-first-core-0/)を保持した。元manifest SHAとHEADのGit blob hashが一致するZIPだけを復元し、core-0全300件を再実行して成功した。初回の成功分と再実行を最終件数へ重複加算していない。数値コードや検査は変更していない。

push後のWindows CIの成否は未確定。Linux成功や過去のCIを新先端Windows成功へ転用しない。時間切れ・取消を成功にしない。
