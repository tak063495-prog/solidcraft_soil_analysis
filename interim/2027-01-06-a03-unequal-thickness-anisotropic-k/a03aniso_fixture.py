"""A-priori unequal layers and explicit constant directional initial stresses."""
import hashlib, itertools, json, math
import numpy as np
from a03aniso_binding import study
HEIGHTS=(.4,.5,.6)
STRESSES={'k55':(.5,.5),'k58':(.5,.8),'k85':(.8,.5)}

def plan():
    specs=[]
    for h,s,n in itertools.product(HEIGHTS,STRESSES,(4,8,16)):
        specs.append(dict(name=f'h{int(h*10)}-{s}-n{n}-gl5',height=h,stress=s,n=n,volume=5,
                          minimum=1,reverse=False,element=False,face=False,engines=['native','independent']))
    for h,s in itertools.product(HEIGHTS,STRESSES):
        specs.append(dict(name=f'h{int(h*10)}-{s}-n8-gl7',height=h,stress=s,n=8,volume=7,
                          minimum=1,reverse=False,element=False,face=False,engines=['independent']))
    return dict(passed=False,a03_complete=False,diagnostic_only=True,planned_before_solves=True,
                formal_paths=63,formal_states=2304,common_coordinates='4 phases x i/4, i=1..4; no interpolation',specs=specs)

def fixture(n,height,stress):
    if n not in (4,8,16) or height not in HEIGHTS or stress not in STRESSES:
        raise ValueError('outside fixed diagnostic matrix; do not solve')
    m,g,c=study.fixture(n,2.)
    xyz=np.asarray(m['nodes'],float); z=xyz[:,2].copy()
    xyz[:,2]=np.where(z<=.5,2*height*z,height+2*(1-height)*(z-.5))
    m['nodes']=xyz.tolist()
    for r in m['regions']:
        cells=[v for v,rid in zip(m['tetrahedra'],m['region_ids']) if rid==r['id']]
        nodes=sorted({v for cell in cells for v in cell})
        r['bounds_m']=[xyz[nodes].min(axis=0).tolist(),xyz[nodes].max(axis=0).tolist()]
        r['volume_m3']=math.fsum(abs(float(np.linalg.det((xyz[cell][1:]-xyz[cell][0]).T)))/6 for cell in cells)
        r['source_sha256']=hashlib.sha256(json.dumps(dict(diagnostic='a03aniso input z transform',height=height,region=r['id']),sort_keys=True).encode()).hexdigest()
    kx,ky=STRESSES[stress]
    def sigma(v):return [kx*v,ky*v,v,0.,0.,0.]
    states=g['initialization']['geostatic']['seed']['states']
    states['LowerClay']['effective_stress_compression_positive_kpa']=sigma(90.)
    states['UpperClay']['effective_stress_compression_positive_kpa']=sigma(70.)
    c['stages'][0]['states']['UpperClay']['effective_stress_compression_positive_kpa']=sigma(5.)
    audit(m,g,c)
    return m,g,c

def audit(m,g,c):
    xyz=np.asarray(m['nodes']); assert xyz.shape==(54,3) and len(m['tetrahedra10'])==12
    base=study.fixture(4,2.)
    for key in ('tetrahedra','tetrahedra10','material_ids','region_ids','face_sets'): assert m[key]==base[0][key]
    edges=((0,1),(1,2),(2,0),(0,3),(1,3),(2,3))
    defect=max(float(np.max(np.abs(xyz[cell[4+i]]-(xyz[cell[a]]+xyz[cell[b]])/2))) for cell in m['tetrahedra10'] for i,(a,b) in enumerate(edges))
    assert defect<1e-15
    assert c['settings']['kinematics']=='side_rollers' and c.get('boundary') is None
    assert g['initialization']['water_level_m']==c['stages'][0]['hydraulic_initialization']['level_m']==2.
    assert g['initialization']['geostatic']['hydrostatic_face_sets']==['exterior'] and xyz[:,2].max()==1.
    assert c['joint_system']==base[2]['joint_system'] and c['initially_inactive_bodies']==['UpperClay']
    assert g['materials']==base[1]['materials']
    assert c['stages'][0]['states']['UpperClay']['last_return']=='elastic'
    volumes={r['body']:r['volume_m3'] for r in m['regions']}
    return dict(allowed=True,nodes=54,tet10=12,connectivity_unchanged=True,midpoint_maximum_defect_m=defect,
                bounds_m=[xyz.min(axis=0).tolist(),xyz.max(axis=0).tolist()],body_volumes_m3=volumes,
                total_volume_m3=math.fsum(volumes.values()),regions=m['regions'],explicit_water_level_m=2.,
                constant_stresses=g['initialization']['geostatic']['seed']['states'],fresh=c['stages'][0]['states'],
                fresh_mc=True,single_continuous_b=True,side_rollers=True,automatic_k0=False,depth_continuous_stress=False)
