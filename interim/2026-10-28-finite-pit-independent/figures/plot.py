from pathlib import Path
import json,sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
p=Path(__file__).resolve().parent;data=json.loads((p/'figure-data.json').read_text(encoding='utf-8'));output=p/'finite-pit-independent.png';assert not output.exists()
font=FontProperties(fname='C:/Windows/Fonts/meiryo.ttc');legend=FontProperties(fname='C:/Windows/Fonts/meiryo.ttc',size=8)
plt.rcParams.update({'font.family':font.get_name(),'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'text.color':'#243746','axes.labelcolor':'#243746','xtick.color':'#243746','ytick.color':'#243746'})
fig,axes=plt.subplots(1,3,figsize=(14.5,5.7));x=[0,1,2]
fig.suptitle('有限幅・二層の局所掘削：同一メッシュの独立比較',fontproperties=font,fontsize=17,y=.97)
names=['総変位','応力','先行圧密圧力 pc','間隙比 e','塑性ひずみ'];colors=['#2166ac','#dc7b28','#5e8b3e','#8053a3','#bf4957']
for (key,values),name,color in zip(data['state_error_percent'].items(),names,colors):axes[0].plot(x,values,'o-',color=color,label=name,lw=1.8,ms=5)
axes[0].set_title('全工程の最大体積L2誤差',fontproperties=font,pad=14);axes[0].set_yscale('log');axes[0].set_ylim(1e-6,.3);axes[0].set_ylabel('独立200分割に対する相対差（％）',fontproperties=font)
axes[0].legend(prop=legend,loc='upper center',bbox_to_anchor=(.5,-.18),ncols=2,frameon=False)
v=data['increment_error_percent'];axes[1].plot(x,v,'o-',color='#2166ac',lw=2,ms=6,label='全工程の増分変位')
axes[1].plot([1,2],[v[0]/2,v[1]/2],'s--',color='#9b9b9b',ms=5,label='前の誤差の1/2')
axes[1].set_title('追加の工程別増分ゲート',fontproperties=font,pad=14);axes[1].set_yscale('log');axes[1].set_ylim(.0008,.03);axes[1].set_ylabel('最大体積L2差（％）',fontproperties=font)
for i,z in enumerate(v):axes[1].annotate(f'{z:.6f}%',(i,z),xytext=(0,10),textcoords='offset points',ha='center',fontsize=9)
axes[1].text(.02,.05,'16→32：未達を保存\n32→64：原閾値で合格',transform=axes[1].transAxes,fontproperties=legend)
axes[1].legend(prop=legend,loc='upper center',bbox_to_anchor=(.5,-.18),frameon=False)
for key,label,color in [('excavation_increment_mm','掘削増分の底面平均','#2166ac'),('cumulative_display_mm','累積表示の底面平均','#bf4957')]:axes[2].plot(x,data['floor_error_nm'][key],'o-',color=color,lw=1.8,ms=5,label=label)
axes[2].set_title('累積平均は工程誤差が相殺',fontproperties=font,pad=14);axes[2].set_yscale('log');axes[2].set_ylim(.015,100);axes[2].set_ylabel('独立解との差の絶対値（nm）',fontproperties=font)
axes[2].annotate('1/32で約208倍の相殺\n累積値の単調半減は未達',(1,data['floor_error_nm']['cumulative_display_mm'][1]),xytext=(.02,.3),textcoords='axes fraction',arrowprops={'arrowstyle':'-','color':'#bf4957'},fontproperties=legend,color='#9e3542')
axes[2].legend(prop=legend,loc='upper center',bbox_to_anchor=(.5,-.18),frameon=False)
for ax in axes:
 ax.set_xticks(x,['1/16','1/32','1/64']);ax.set_xlim(-.2,2.2);ax.set_xlabel('施工増分幅',fontproperties=font,labelpad=8);ax.grid(axis='y',alpha=.2)
fig.text(.055,.045,'274 Tet10・128板。独立RK4 100/200分割。仮想材料・排水・小ひずみ。同一離散モデルの照合で、空間収束の認定ではありません。',fontproperties=legend)
fig.subplots_adjust(left=.065,right=.98,top=.78,bottom=.29,wspace=.42)
fig.savefig(output,dpi=175,facecolor='white');plt.close(fig)
