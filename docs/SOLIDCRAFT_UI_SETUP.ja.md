# SolidCraft画面内連携の作業環境

2026-10-24。計画書P5用に、固定版SolidCraftから独立したGit作業ツリーを用意した。画面内連携の実装・動作完了を示す記録ではない。

- 固定元：`external/solidcraft`、`edf57cec293dcea301201f1d5173b33e16351466`（sources.lock.jsonと一致）。
- 作業先：`external/solidcraft-soil-fem`。
- 作業ブランチ：`codex/soil-fem-ui`。固定版から開始したローカル専用ブランチで、公開PR/issueへの変更はしていない。
- 現行の取得・解析は引き続きportable 0.20.0と外部スクリプトを使う。新しいUI版を既存検証の実行版へ置き換えていない。

SolidCraftはRust workspaceで、`crates/ui-egui`が画面、`crates/engine`がコマンド、`apps/solidcraft`がネイティブ起動を担う。画面追加は共通のUIコンテナとテーマを用い、操作はengineの登録コマンドを通す。プログラム呼出しからダイアログを開かない。WASMにネイティブプロセス依存を持ち込まない。

元の参照用checkoutはsparse設定で`xtask`が省略されていた。新しい作業ツリーだけに`xtask/docs/scripts/tests/examples/packaging`を追加し、`cargo metadata --locked --no-deps`で17個のworkspace packageを確認した。元の参照checkoutとportable実行版には変更していない。UI版のビルド・画面試験はこれから行う。

次の実装対象は、材料・初期状態・工程とCADボディの対応、版を固定した解析ジョブの起動/監視/中断、CADまたは解析入力の変更時の結果失効、地盤メッシュ上の3Dコンター/変形/断面。積分点の生値・体積平均・表示用平滑化を区別する。既存のJSON入出力・ハッシュ証跡を基に、UIは薄いRust層として構築する。メッシャーとソルバーは外部ジョブ境界を保つ。

この作業ツリーでは、upstreamの`AGENTS.md`にあるpanic禁止、コマンド登録、レイヤ分離、UIコンテナ、フォーマット移行と画面証跡の規則を適用する。gitignoredの`plan/`は取得元に存在しないため、ユーザーが指定した本プロジェクトの`docs/PLAN.ja.md`と`docs/P6_STATUS.ja.md`を実装範囲と受入ゲートの根拠とする。

将来の受入では、このブランチの差分・固定元・ビルド版・実画面操作と結果表示を保存する。ジョブUIの設計やソース準備だけでP5完了とは判定しない。

固定版の`cargo check --locked -p solidcraft-ui-egui`がWindows x64で合格した（4分48秒、2ビルドworker）。ログは`runs/solidcraft-ui-baseline-check-20261024.log`。コード変更なしの基準ビルド確認であり、解析パネル・ジョブ操作・描画の実装や画面受入はまだ行っていない。

## 解析属性・STEP出力の実装中

同作業ツリーに`analysis.soil.*`コマンド、属性編集窓、論理ボディ名と作成フィーチャー/出力ファミリーの保存、改名時の対応更新、CAD/入力変更の検出と新規フォルダーへのSTEP出力を追加した。ネイティブ型チェック、engine6試験、WASM型チェック、Rust書式・依存層・資産/アイコン・UIコンテナ・一覧順序の検査は合格。[初回の版と成功ログ](../runs/solidcraft-analysis-quality-progress-20261024-v2.json)はその時点の記録として保持する。[実GUI出力と実SolidCraft→FEMの限定検証](SOLIDCRAFT_NATIVE_VERIFICATION.ja.md)に続き、UI3試験・clippyと小さい許容値の表示修正を確認した。外部ジョブワーカーも限定検証したが、nativeジョブ操作/結果UIは未接続でP5完了とは判定しない。

外部取込アダプターの5試験と、合成STEPによる多層MCC局所掘削のコントローラー事前試験は合格した。後者を実SolidCraft出力の証拠にしない。[属性・出力形式・検証の区別](SOLIDCRAFT_ANALYSIS_MODEL.ja.md)を参照。解析ジョブの起動/監視/中断と、画面内の結果コンターは引き続き未実装である。
