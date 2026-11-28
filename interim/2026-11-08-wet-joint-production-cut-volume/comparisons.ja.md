# 生産4経路と既存4096/64参照：全18場

全7工程端における最大値。相対L2は百分率ではない。絶対差最大の工程と相対差最大の工程は異なる場合がある。全履歴の合格基準を新設しない。

## original

|場|最大絶対差|その工程|最大相対L2|その工程|
|---|---:|---|---:|---|
|変位 m|2.41615189758e-06|recontact|0.00396473181455|recontact|
|有効応力 kPa|0.778942508298|recontact|0.0020268256312|recontact|
|ひずみ|1.25017235481e-05|recontact|0.0077530700091|recontact|
|間隙比|9.55378844569e-06|recontact|3.4801389142e-06|recontact|
|材料GP水圧 kPa|3.33066907388e-15|initial_equilibrated|1.65597269382e-16|initial_equilibrated|
|塑性ひずみ|0|initial_equilibrated|0|initial_equilibrated|
|jump m|4.30454632177e-06|recontact|0.00285747160193|shear|
|traction kPa|1.0762899448|recontact|0.00690162069534|shear|
|塑性すべりベクトル m|0.00116451742923|reload|4.52618847412|reverse|
|累積すべり m|0.00289583458599|return|0.713949116286|recontact|
|弾性エネルギー kPa·m|2.90533558795e-05|recontact|0.00270520975366|shear|
|塑性散逸 kPa·m|0.03055061936|reload|0.728744855708|recontact|
|数値散逸 kPa·m|0.0228928843754|return|75974.8112009|shear|
|端点仕事 kPa·m|0.0167910594367|recontact|0.462613169961|open|
|接合面物理水圧 kPa|0|initial_equilibrated|0|initial_equilibrated|
|節点水圧荷重 kN|4.4408920985e-16|initial_equilibrated|5.23021967397e-16|initial_equilibrated|
|節点内力−外力 kN|0.0179112679656|recontact|0.00219551077179|recontact|
|gap m|2.12923284414e-07|recontact|0.00131906886816|shear|

## g256-m5e-6

|場|最大絶対差|その工程|最大相対L2|その工程|
|---|---:|---|---:|---|
|変位 m|8.41885839081e-09|shear|4.60218573745e-06|shear|
|有効応力 kPa|0.00215392227614|shear|5.4400118478e-06|shear|
|ひずみ|2.68371234633e-08|shear|1.86127443816e-05|shear|
|間隙比|2.4361965778e-08|shear|8.38603876829e-09|shear|
|材料GP水圧 kPa|3.33066907388e-15|initial_equilibrated|1.65597269382e-16|initial_equilibrated|
|塑性ひずみ|0|initial_equilibrated|0|initial_equilibrated|
|jump m|1.45474697654e-08|shear|1.45879318708e-05|shear|
|traction kPa|0.00277596293654|shear|2.86234852204e-05|shear|
|塑性すべりベクトル m|2.73548762664e-05|open|0.132053177317|reverse|
|累積すべり m|4.66903447573e-05|return|0.0116687109334|recontact|
|弾性エネルギー kPa·m|3.27591761712e-08|shear|9.81863436216e-06|shear|
|塑性散逸 kPa·m|0.000174239199911|recontact|0.00410782240093|open|
|数値散逸 kPa·m|0.000305584789945|return|1022.99972786|initial_equilibrated|
|端点仕事 kPa·m|0.000395131370267|return|0.00390242041615|initial_equilibrated|
|接合面物理水圧 kPa|0|initial_equilibrated|0|initial_equilibrated|
|節点水圧荷重 kN|4.4408920985e-16|initial_equilibrated|5.23021967397e-16|initial_equilibrated|
|節点内力−外力 kN|6.08782691423e-05|shear|3.70412824519e-06|shear|
|gap m|4.68135489798e-10|shear|4.80322586702e-06|shear|

## g512-m5e-6

|場|最大絶対差|その工程|最大相対L2|その工程|
|---|---:|---|---:|---|
|変位 m|1.57318741182e-08|shear|7.89502648136e-06|shear|
|有効応力 kPa|0.00375534164098|shear|9.37704716895e-06|shear|
|ひずみ|4.71605841764e-08|shear|3.23910383226e-05|shear|
|間隙比|4.30473292568e-08|shear|1.44773015046e-08|shear|
|材料GP水圧 kPa|3.33066907388e-15|initial_equilibrated|1.65597269382e-16|initial_equilibrated|
|塑性ひずみ|0|initial_equilibrated|0|initial_equilibrated|
|jump m|2.7506314727e-08|shear|2.53406178955e-05|shear|
|traction kPa|0.00503704347111|shear|5.11464261616e-05|shear|
|塑性すべりベクトル m|1.21751486046e-05|open|0.0508328859799|reverse|
|累積すべり m|2.23481063122e-05|return|0.00492154620414|recontact|
|弾性エネルギー kPa·m|5.81526859689e-08|shear|1.75213036467e-05|shear|
|塑性散逸 kPa·m|5.24739917009e-05|open|0.00194223205316|open|
|数値散逸 kPa·m|0.000172071572159|return|510.999863932|initial_equilibrated|
|端点仕事 kPa·m|0.000202305382125|return|0.0019493028667|initial_equilibrated|
|接合面物理水圧 kPa|0|initial_equilibrated|0|initial_equilibrated|
|節点水圧荷重 kN|4.4408920985e-16|initial_equilibrated|5.23021967397e-16|initial_equilibrated|
|節点内力−外力 kN|0.000101818542277|shear|5.63980545396e-06|shear|
|gap m|8.31017188426e-10|shear|8.57166412616e-06|shear|

## g512-m2.5e-6

|場|最大絶対差|その工程|最大相対L2|その工程|
|---|---:|---|---:|---|
|変位 m|4.15418244652e-09|shear|2.26567910122e-06|shear|
|有効応力 kPa|0.00106354043939|shear|2.67695125321e-06|shear|
|ひずみ|1.31999759205e-08|shear|9.16735013521e-06|shear|
|間隙比|1.20127416903e-08|shear|4.12499385473e-09|shear|
|材料GP水圧 kPa|3.33066907388e-15|initial_equilibrated|1.65597269382e-16|initial_equilibrated|
|塑性ひずみ|0|initial_equilibrated|0|initial_equilibrated|
|jump m|7.17833306782e-09|shear|7.17203197371e-06|shear|
|traction kPa|0.00137213958447|shear|1.41398908303e-05|shear|
|塑性すべりベクトル m|1.21751506718e-05|open|0.0508316608326|reverse|
|累積すべり m|2.2320816965e-05|return|0.00492082866981|recontact|
|弾性エネルギー kPa·m|1.61695889731e-08|shear|4.85017612969e-06|shear|
|塑性散逸 kPa·m|5.01490318296e-05|return|0.00119666073291|open|
|数値散逸 kPa·m|0.000152395390892|return|510.999863932|initial_equilibrated|
|端点仕事 kPa·m|0.000201446376901|return|0.0019493028667|initial_equilibrated|
|接合面物理水圧 kPa|0|initial_equilibrated|0|initial_equilibrated|
|節点水圧荷重 kN|4.4408920985e-16|initial_equilibrated|5.23021967397e-16|initial_equilibrated|
|節点内力−外力 kN|2.98707810344e-05|shear|1.79824109029e-06|shear|
|gap m|2.310662461e-10|shear|2.3721659544e-06|shear|
