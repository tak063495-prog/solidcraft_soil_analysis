# Issue #35 C-01 global-step の物理監査

PR #13 の起点 `8cf97d2bb4c3e5b25667445f6c0033379d565ec1` で、実際の `study.audit` に Tet4 体積要素 1 個、接触積分点 1 個、全 7 工程に各 1 状態を与えた。native と independent の両方で `audit.py:38` の最初の体積要素の応力照査に達し、`NameError: name 'fem' is not defined` を確認した。 [Issue #35](https://github.com/tak063495-prog/solidcraft_soil_analysis/issues/35)・[再現手順](REPRODUCE.ja.md)

## 修正と回帰試験

global-step の `audit.py` で、2026-12-12 の switch-event 版と同じく `from study import …, fem` を使う。動的ロード先で `fem.elastic_tensor` を解決できる。変更は import 1 行。FEM の計算、監査の許容値、参照値と保存済みの判定は従来どおり。

同じ合成履歴で両経路とも 7 状態を読み、再組立力残差と弾性応力の整合誤差は 0。応力の 1 成分を 1 kPa ずらすと既存の応力照査が拒否した。試験では合成の釣合い演算子を渡し、実際の保存・読込、動的ロード、体積要素の応力照査を通す。

## 保存記録と P6-01

旧監査ソースの blob は `c270508abe7a4e9f8286cd3c2320a4f632315333`、修正後は `f3052b8b4af3ef7f3eb36579eff948fcfcc25afe`。旧コミットの実バイトを保存済みの 2 本のソース指紋と照合し、現行の実バイトは今回の Git blob ID と SHA-256 で個別に照合する。 [修正前後の指紋](acceptance.json)

P6-01 の追補台帳に `c01_global_step_issue35_20270111_acceptance` を追加した。修正コミット `04bc2e6d128289febd82fdd2e5a87f5103a83d11` の `acceptance.json` を 40 桁コミット、Git blob ID、SHA-256、型付きの状態値で固定した。旧追補台帳のバイト列を接頭部として保持した。 [照合手順](REPRODUCE.ja.md)

## Linux の実測

| 対象 | passed | failed・error | skipped |
|---|---:|---:|---:|
| Rust 全ターゲット | 467 | 0 | 0 |
| Rust doc-test | 0 | 0 | 0 |
| Python CI 7 群（tests・tools） | 986 | 0 | 2 |
| Linux 専用 pytest | 7 | 0 | 0 |
| Python 合計 | 993 | 0 | 2 |

今回の回帰試験は pytest で 5 passed、0 failed、0 skipped。6 サブテストも成功した。上表の Python CI 7 群に含まれるため、合計では重複計上していない。Python の skip 2 件は Linux に PowerShell がないため。全件の未実行は 0。P6-01 の 373 出典で整合性を確認した。 [全件の内訳](REPRODUCE.ja.md)

A-01〜A-04、B-01〜B-03、C-01〜C-04、NC-01、P6 は未達、passed=false。
