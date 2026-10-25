"""Deterministic hydraulic examples on existing SolidCraft-compatible meshes."""
from pathlib import Path
import copy
import json


def cases():
    retention=dict(residual_water_content=.08,saturated_water_content=.45,alpha_per_m=1.2,n=1.6,saturated_specific_storage_per_m=1e-5)
    materials=[dict(id=i,saturated_conductivity_m_s=[[1e-6,0.,0.],[0.,1e-6,0.],[0.,0.,1e-6]],retention=retention) for i in (1,2)]
    def boundary(head,rain=0.):return [dict(face_set='bottom',condition=dict(kind='head',value=dict(head_m=head))),
        dict(face_set='top',condition=dict(kind='atmospheric_no_ponding',rainfall_vertical_m_s=rain))]
    base=dict(schema='soil-fem-richards/1',units='m-s',materials=materials,initial_head=dict(head_m=2.),analysis=dict(kind='steady',boundaries=boundary(2.)))
    examples={'richards-steady-wet':copy.deepcopy(base)}
    base['analysis']['boundaries']=boundary(0.)
    examples['richards-steady-dry']=copy.deepcopy(base)
    base['initial_head']['head_m']=0.
    base['analysis']=dict(kind='transient',periods=[dict(name=name,end_time_s=end,boundaries=boundary(0.,rain))
        for name,end,rain in [('rain',3600.,3e-6),('pause',7200.,0.),('rain_again',10800.,3e-6)]])
    examples['richards-rain-pause']=copy.deepcopy(base)
    return examples


if __name__=='__main__':
    folder=Path(__file__).resolve().parents[1]/'examples'
    for name,case in cases().items():
        (folder/(name+'.json')).write_text(json.dumps(case,indent=2)+'\n',encoding='utf-8')
