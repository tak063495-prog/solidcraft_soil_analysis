# P6-01 追補　C-03 MCC 仕事の検証範囲

|範囲 ID|判定|根拠・適用範囲|
|---|---|---|
|work_evaluated|達成|本体・独立の6経路、324 GP・224受理座標。全材料内分割の終端状態を保存状態と完全一致照合。|
|decomposition_evaluated|達成|独立に計算した体積・せん断の和との差を、材料内分割・受理増分・全履歴で計測。|
|one_axis_comparisons|達成|同じ時刻の材料2対、接触2対、実装間3対の差を計測。|
|quantitative_work_acceptance|未達|参照なし・未達。仕事量・分解残差・実装差の承認済み定量許容幅なし。|
|force_all_points|未達|元基準1e−10を本体で各経路117/224点超過。独立の3経路は超過0。|
|general_c03|未達|一般材料・別荷重・開口・全場の刻み収束および元 C-03 全体。|

[全受理点の力残差監査](results/force-audit.json)・[数値と機械可読の状態](results/validation.json)・[P6-01 固定証拠台帳](../2026-12-10-ci-p601-inventory-fix/inventory.toml)
