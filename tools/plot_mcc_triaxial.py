"""Standalone scientific figure from sealed material-point verification data.

Requires Matplotlib and NumPy only for plotting; solver/runtime dependencies
are unchanged. Keep figures in a new directory separate from sealed evidence.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def plot(run,output):
    output.mkdir(parents=True,exist_ok=False);sealed=json.loads((run/'final-validation.json').read_text());assert sealed['passed']
    data=run/'paths';report=json.loads((data/'verification.json').read_text());m=1.2
    figure,axes=plt.subplots(2,2,figsize=(10,8),constrained_layout=True)
    for row in report['cases']:
        name=row['case'];oc=name.startswith('oc');drained=name.endswith('-drained');axis=axes[int(oc),int(drained)]
        native=json.loads((data/(name+'-fine-result.json')).read_text())['states'];reference=json.loads((data/(name+'-fine-reference.json')).read_text())['states']
        actual=[]
        for s in native:
            stress=np.asarray(s['effective_stress_compression_positive_kpa']);p=stress[:3].mean();dev=stress.copy();dev[:3]-=p
            q=np.sqrt(1.5*(dev[:3]@dev[:3]+2*(dev[3:]@dev[3:])));actual.append([p,q])
        actual=np.asarray(actual);reference=np.vstack(([100.,0.,row['ocr']*100.,native[0]['void_ratio'],0.],reference))
        pc0=row['ocr']*100;pc1=native[-1]['preconsolidation_pressure_kpa'];limit=max(pc0,pc1)*1.04
        x=np.linspace(0,limit,500);axis.plot(x,m*x,color='#555555',linestyle='-.',linewidth=1,label='Critical state q = M p')
        for pc,style,label in ((pc0,'--','Initial yield ellipse'),(pc1,':','End yield ellipse')):
            p=np.linspace(0,pc,500);axis.plot(p,m*np.sqrt(np.maximum(p*(pc-p),0)),color='#888888',linestyle=style,linewidth=1,label=label)
        axis.plot(actual[:,0],actual[:,1],color='#176bb0',linewidth=2,label='Rust, fine increments')
        axis.plot(reference[:,0],reference[:,1],color='#d46100',linestyle='--',linewidth=1.3,label='Independent midpoint reference')
        axis.scatter([row['theoretical_critical_p_kpa']],[row['theoretical_critical_q_kpa']],marker='x',s=48,color='black',zorder=6,label='Theoretical asymptotic limit')
        axis.scatter([100],[0],color='#176bb0',s=18,zorder=6)
        axis.set(xlim=(0,limit),ylim=(0,max(m*max(pc0,pc1)/2,row['theoretical_critical_q_kpa'])*1.15),
                 xlabel="Mean effective pressure p (kPa)",ylabel='Deviatoric stress q (kPa)',
                 title=('OC, pc0/p0 = 4' if oc else 'NC, pc0/p0 = 1')+(' | effective sigma3 = 100 kPa' if drained else ' | constant volume'))
        axis.grid(alpha=.18)
    handles,labels=axes[0,0].get_legend_handles_labels();figure.legend(handles,labels,loc='outside lower center',ncol=3,fontsize=8)
    figure.suptitle('Modified Cam-Clay: independent triaxial material-point verification',fontsize=13)
    figure.savefig(output/'mcc-triaxial.png',dpi=170);figure.savefig(output/'mcc-triaxial.svg');plt.close(figure)
    shutil.copyfile(Path(__file__),output/'plot-code.py')
    metadata={'source_final_validation_sha256':hashlib.sha256((run/'final-validation.json').read_bytes()).hexdigest(),
        'source_run':str(run.resolve()),'matplotlib_version':matplotlib.__version__,
        'meaning':'p-q material-point paths. End state is at axial strain 0.5; theoretical markers are asymptotic limits. This figure is not a global small-strain FEM deformation analysis.',
        'sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in output.iterdir() if p.is_file()}}
    with (output/'figure-provenance.json').open('x',encoding='utf-8') as f:json.dump(metadata,f,indent=2)
    print(output/'mcc-triaxial.png')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();plot(a.run.resolve(),a.output.resolve())
