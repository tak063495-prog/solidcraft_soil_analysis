"""Render the supplied diagnostic data; this plot asserts no convergence."""
import argparse,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
parser=argparse.ArgumentParser();parser.add_argument('data',type=Path);parser.add_argument('output',type=Path);args=parser.parse_args()
assert not args.output.exists();data=json.loads(args.data.read_text(encoding='utf-8'))
font=FontProperties(fname='C:/Windows/Fonts/meiryo.ttc')
plt.rcParams.update({'font.family':font.get_name(),'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'axes.labelcolor':'#243746','text.color':'#243746','xtick.color':'#243746','ytick.color':'#243746'})
fig,(left,right)=plt.subplots(1,2,figsize=(12,5),gridspec_kw={'width_ratios':[1,1.35]})
fig.suptitle('MCC 局所掘削：空間収束は未達',fontproperties=font,fontsize=18,y=.98)
left.bar([0,1],data['relative_l2_percent'],color=['#2166ac','#4d8abe'],width=.55,zorder=3)
left.set_xticks([0,1],data['mesh_pairs']);left.set_ylim(0,22)
left.set_ylabel('掘削増分変位の体積L2差（％）',fontproperties=font)
left.set_title('地盤と壁を同時に細分化',fontproperties=font,pad=14)
left.grid(axis='y',alpha=.2,zorder=0)
for i,value in enumerate(data['relative_l2_percent']):left.text(i,value+.7,f'{value:.2f}%',ha='center',fontweight='bold')
left.set_xlabel('四面体要素数',fontproperties=font,labelpad=10)
colors=['#2166ac','#8fa8bd','#d3dae2'];offset=[0.,0.]
labels=['底面直下 Lower11','上層 UpperSoil','その他の下層8ボディ']
for i,label in enumerate(labels):
 values=[data['volume_percent'][i],data['squared_difference_percent'][i]]
 right.barh([0,1],values,left=offset,height=.5,color=colors[i],label=label)
 if i==0:
  for j,v in enumerate(values):right.text(v/2,j,f'{v:.1f}%',color='white',ha='center',va='center',fontweight='bold',fontsize=10)
 offset=[x+y for x,y in zip(offset,values)]
right.set_yticks([0,1],['有効地盤の体積','二乗変位差の積分']);right.set_xlim(0,100)
right.set_xlabel('全体に占める割合（％）',fontproperties=font,labelpad=10)
right.set_title('945→1324要素：差は底面直下に集中',fontproperties=font,pad=14)
right.legend(loc='upper center',bbox_to_anchor=(.5,-.22),frameon=False,prop=FontProperties(fname='C:/Windows/Fonts/meiryo.ttc',size=9))
fig.text(.05,.025,'仮想材料・排水・小ひずみ。離散解どうしの差であり、連続体誤差の認定値ではありません。',fontproperties=font,fontsize=9)
fig.subplots_adjust(top=.78,bottom=.29,left=.075,right=.98,wspace=.56)
fig.savefig(args.output,dpi=175,facecolor='white');plt.close(fig)
