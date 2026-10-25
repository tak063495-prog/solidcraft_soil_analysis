# 3D Mohr–Coulombの操作と範囲

MCCに加え、等方弾性・完全塑性の3D Mohr–Coulomb（MC）を材料モデルとして選択できる。Tet4/Tet10、共通の全体Newton・非対称CSR/GMRES、一般境界条件、載除荷・局所掘削、履歴切戻し・保存再開に接続した。MCCとMCは同じメッシュ内で混在可能。

## 実行

```powershell
cd D:\Work\soil-fem-3d
# SolidCraftの2層＋PitをMCで解析
.\scripts\Invoke-LocalExcavation.ps1 -GroundFile examples\ground-mc-local.json -StageFile examples\stages-mc-local.json -MeshSizeM 4
# 材料点の定体積せん断と除荷
target\release\soil-fem-3d.exe mc-path examples\mc-shear-path.json runs\NEW_MC_PATH.json
# 材料点、柱理論、CAD施工、MCC互換の統合検証
.\scripts\Test-MohrCoulomb.ps1
```

各実行で新しいフォルダを作る。初期地盤・施工スクリプトは使用するexeを`solver/soil-fem-3d.exe`へコピーし、SHA256を確認して、そのジョブ内の版で解析・検証する。計算中に開発版をビルドしても実行ファイルが置き換わらない。`-SolverPath`で別の検証済みexeを指定して比較することもできる。保存再開には保存ジョブ内の同じ版を使用する。

## 材料と初期状態

`ground.materials`の各材料に`mcc`または`mohr_coulomb`のどちらか一方を指定する。同時指定・未指定・未知フィールドを拒否する。既存MCC入力のキー・状態JSONは維持した。

```json
{
  "id": 1, "name": "MC soil", "unit_weight_kn_m3": 20.0,
  "mohr_coulomb": {
    "young_kpa": 50000.0, "poisson": 0.25,
    "cohesion_kpa": 10.0, "friction_deg": 30.0, "dilation_deg": 5.0,
    "tension_cutoff_kpa": 0.0, "flow": "non_associated"
  }
}
```

E>0、−1<ν<0.5、c≥0、0≤ψ≤φ<89°を要求する。負のダイレイタンシーや圧縮キャップ、応力依存剛性、せん断/体積硬化はこのモデルに含めない。これらは計画P2bの別モデル。

MCの地層にはK0x/K0yと`consolidation: {"definition":"not_applicable"}`、`void_ratio: {"mode":"explicit","value":0.7}`を指定する。MCへ架空のpc/OCRを設定しない。MCの間隙比は全体積ひずみによる比体積の変化を追跡するが、E・νや強度を圧縮硬化させる量ではない。K0初期応力が降伏面・引張カットオフを超える場合は入力を拒否する。

水平直方体層のK0初期化に加え、2026-10-24に一般形状の種応力読込・非線形地圧平衡を追加した。MC/MCCの初期化中の塑性履歴を施工へ引き継ぐ。[一般初期平衡](GEOSTATIC.ja.md)を参照。

## 構成式と帰着

応力は有効応力・圧縮正、配列はxx,yy,zz,xy,yz,zx。せん断ひずみは工学せん断ひずみ。固有値を圧縮の大きい順σ1≥σ2≥σ3に並べる。

MCの6平面を`f_ij=(1−sinφ)σi−(1+sinφ)σj−2c cosφ≤0`（i≠j）とし、塑性ポテンシャルの方向は同じ式のφをψへ置き換えた勾配を使う。弾性予測応力の主軸内で、非負の塑性乗数と全平面の降伏条件・相補条件を満たすactive setを探索する。3次元の補正ベクトルは最大3個の生成方向で表現できるため、1〜3個の組合せを列挙する。よく使う面/稜線を先に判定するが、残る組合せを省略しない。

面、圧縮稜線、伸張稜線、MC頂点、引張面/稜線/頂点、MCと引張の混合帰着を区別する。弾性予測主軸から応力と塑性ひずみテンソルを全体座標へ戻す。φ=0はTrescaとして扱う。有限のサブステップに分け、全増分写像の6×6数値接線を、固定分割数で中心差分する。非関連流れの非対称性を保存する。切替点の接線は一般化された近似であり、一意な古典微分を仮定しない。

`tension_cutoff_kpa`が数値なら、各主応力に`−σi−t≤0`という関連Rankine面を追加する。t=0は引張強度0。φ>0ではt≤c cotφを要求し、上限を超える入力を黙って丸めない。nullは純MC。特にψ=0の純MCで頂点より引張側へ進む経路は、非関連流れの補正が存在しない場合がある。その場合は材料エラーとして返し、圧力のクリッピングや別の流れ則で隠さない。

`flow`は次の3方式。

|値|挙動|
|---|---|
|`non_associated`|指定したφとψを使う。ψ=0も可能|
|`associated`|ψ=φを要求する|
|`davis`|β=cosφ cosψ/(1−sinφ sinψ)、c*=βc、tanφ*=β tanφ、ψ*=φ*へ変換して関連流れを使う|

DAVISは降伏強度と流れ則を変える比較方式。元の非関連モデルと同一物性・同一安全率として扱わない。保存した入力には元のc/φ/ψと`flow`を残す。

## 履歴、出力、再開

MC状態は応力、塑性ひずみ、間隙比、相当塑性ひずみ、塑性仕事、最終帰着種別を保存する。塑性仕事はサブステップ終点の応力と塑性ひずみ増分で評価し、有限刻みに依存する近似。材料点・全体試行とも確定状態を直接変更しない。棄却時の塑性ひずみ変化を`rejected_trial_plastic_changed`へ記録し、実際に塑性化した試行を棄却したことも検査する。

MCを含む施工結果は`soil-fem-stage-result/3`、初期結果は`soil-fem-initial-state/2`。初期点は`mohr_coulomb_state`を持ち、MCCは従来の`mcc_state`。施工の状態配列は各材料に一致する履歴型を持つ。保存再開で型が不一致なら拒否する。MCCのpcをMCへ流用しない。

VTUには全3方向変位と体積平均状態、生の積分点状態を分けて出力する。`material_model`はMCC=0、MC=1。混在モデルで他の材料に存在しない履歴値はNaNで、0という物性値を意味しない。Z上向き正。沈下は負、隆起は正で、監視平均は節点算術平均。

## 検証の解釈

材料点は別実装の反復相補条件法、全体解は横ひずみを全節点で拘束した2材料柱の閉形式、3D局所掘削は独立した力/モーメント・荷重解放・降伏条件・履歴凍結で検査する。既存MCC結果との完全一致も確認する。[実測と未解決事項](MC_VERIFICATION.ja.md)を参照。

完全塑性MCは、塑性機構や局所化で変位が一意に定まらない場合がある。力の収束だけで変位精度や安全率を認定しない。柱の1D理論比較は横ひずみを明示拘束した基準問題であり、一般の3D塑性変位の独立解照合とは区別する。数値非収束から物理的な崩壊荷重やFs上端を作る機能は実装していない。2D版の実行結果との全面照合、一般地形の初期平衡、非圧縮/局所化の評価、SRMは後続要件。

定式化の参考：[Sysala・Cermak, arXiv:1508.07435v3](https://arxiv.org/abs/1508.07435v3)、[ItascaのMCと引張カットオフ](https://docs.itascacg.com/itasca950/common/models/mohr/doc/modelmohr.html)。2D版DAVISの比較対象は固定コミット`440d246d9e32f5892c79539d4a7b797ad29442ee`の`FEMEngine.bas`。上流コードは変更していない。
