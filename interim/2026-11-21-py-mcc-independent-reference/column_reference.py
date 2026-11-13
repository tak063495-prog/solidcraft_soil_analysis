"""Independent 1D constrained MCC column, documented mixed-MCC test inputs.
Sources S8--S10 in README. No old reference/production solver is imported.
Fixed geometry/self weight and small-strain displacement u=-integral(ev dz).
"""
import math
import numpy as np
from mcc_reference import Material, initial, rate, invariants, richardson, yield_value


def load_point(z, n=1024, mat=Material()):
    if not 0 <= z <= 1 or n < 1:
        raise ValueError('invalid height/increment count')
    # Identical material in both .5 m layers, gamma=20, P0=100, K0=1, OCR=1.
    y=initial(mat,p=100+20*(1-z))
    v0=float(y[3])
    h=30/n
    for _ in range(n):
        dy=rate(y,'k0',mat)
        dsv=dy[0]+2*dy[1]/3
        if dsv <= 0:
            raise ValueError('nonpositive constrained tangent')
        mid=y+.5*h*dy/dsv
        mid[3]=v0*math.exp(-mid[4])
        dy=rate(mid,'k0',mat)
        dsv=dy[0]+2*dy[1]/3
        if dsv <= 0:
            raise ValueError('nonpositive constrained midpoint tangent')
        y=y+h*dy/dsv
        y[3]=v0*math.exp(-y[4])
    return y


def unload_point(y, pressure_change, mat=Material()):
    if pressure_change > 0:
        raise ValueError('elastic unloading requires nonpositive stress change')
    p,q,lpc,v,ev=map(float,y)
    dp=pressure_change/(1+4*mat.g_over_k/3)
    pn=p+dp
    if pn <= 0:
        raise ValueError('p<=0 during unloading')
    vn=v-mat.kap*math.log(pn/p)
    out=np.array([pn,q+2*mat.g_over_k*dp,lpc,vn,ev+math.log(v/vn)])
    if yield_value(out,mat)>1e-8:
        raise ValueError('elastic unload leaves yield surface; unsupported branch')
    return out


def point_history(z,n=1024):
    m=Material()
    seed=initial(m,p=100+20*(1-z))
    load=load_point(z,n,m)
    unload=unload_point(load,-30.,m)
    cut=unload_point(unload,-10.,m) if z <= .5 else unload.copy()
    return [seed,load,unload,cut]


def integrated_stages(z, n=1024, quadrature=16, cells=None):
    """All stages on remaining 0<=z<=.5 domain, integrating each history once."""
    if not 0 <= z <= .5:
        raise ValueError('outside common remaining domain')
    if z==0:
        return [0.,0.,0.]
    if cells is None:
        x,w=np.polynomial.legendre.leggauss(quadrature)
        points=(x+1)*z/2;weights=w*z/2
    else:
        points=(np.arange(cells)+.5)*z/cells
        weights=np.full(cells,z/cells)
    values=np.zeros(3)
    for p,w in zip(points,weights):
        hs=point_history(float(p),n)
        values-=float(w)*np.array([y[4] for y in hs[1:]])
    return values.tolist()


def diagnose(material_levels=(256,512,1024),spatial_cells=(8,16,32),spatial_material_n=2048):
    stage_names=('initial','load','unload','cut')
    report=dict(source_ids=['S1','S3','S8','S9','S10'],accepted=None,tolerance='未定義',
                geometry_m=dict(width=1.,depth=1.,height=1.,layer_interface=.5),
                initial=dict(surcharge_kpa=100.,unit_weight_kn_m3=20.,k0=1.,pc_over_p=1.),
                pressure_kpa=[100.,130.,100.,100.],cut_height_m=.5,
                kinematics='fixed domain; small strain; u=-integral(eps_v dz); logarithmic v evolution',
                material_levels=list(material_levels),spatial_cells=list(spatial_cells),
                spatial_fixed_material_n=spatial_material_n)
    # Compare remaining lower-column top at z=.5 for all three loaded stages.
    increments=[integrated_stages(.5,n=n,quadrature=16) for n in material_levels]
    spatial=[integrated_stages(.5,n=spatial_material_n,cells=c) for c in spatial_cells]
    # 16->32 Gauss at fixed n estimates continuum integration error separately.
    gauss16=integrated_stages(.5,n=spatial_material_n,quadrature=16)
    gauss32=integrated_stages(.5,n=spatial_material_n,quadrature=32)
    report.update(material_increment_values_m=increments,spatial_values_m=spatial,
                  material_increment_diagnostic=richardson(increments),
                  spatial_diagnostic=richardson(spatial),
                  gauss16_to_32_max_abs_m=float(np.max(np.abs(np.array(gauss32)-gauss16))),
                  continuum_quadrature_values_m=gauss32,
                  stage_names=stage_names[1:],
                  original_mixed_mcc_failure_preserved=dict(pc_max_abs_kpa=1.773e-3,threshold_kpa=1e-4,
                      later_eight_layer_pc_abs_kpa=5.334e-4,accepted=False,source='S8'),
                  general_3d_spatial_order=None)
    # Explicit states at fixed heights; upper points freeze on removal, no history filling.
    report['point_histories']=[dict(z_m=z,active=[True,True,True,z<=.5],
        rows=[dict(stage=stage,**invariants(y)) for stage,y in zip(stage_names,point_history(z,spatial_material_n))])
        for z in (.125,.375,.625,.875)]
    return report
