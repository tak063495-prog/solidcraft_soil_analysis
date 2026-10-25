"""Optional Matplotlib figure of raw Tet10 nodal uz on excavation surfaces.

The geometry is undeformed. Four display triangles per Tri6 retain the
exported edge-node values without smoothing; colours are displacement in mm.
"""
import argparse
import hashlib
import json
import shutil
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import fem_reference as ref
from verify_acceptance_packet import verify


def plot(run,output):
    integrity=verify(run)
    output.mkdir(parents=True,exist_ok=False)
    fig=plt.figure(figsize=(13,6.2),layout='constrained')
    # Include every raw nodal extremum: the floor mean is substantially below
    # the peak heave, so fixed limits chosen from the mean clip the result.
    extrema=[]
    for name in ('flat','curved'):
        state=json.loads((run/name/'stage-results.json').read_text())['last_committed']
        extrema.extend(1000*np.asarray(state['displacement_m'])[:,2])
    minimum=float(np.min(extrema));maximum=float(np.max(extrema))
    if not minimum<0<maximum:raise ValueError('this comparison plot requires settlement and heave')
    norm=TwoSlopeNorm(vmin=minimum,vcenter=0,vmax=maximum);cmap=plt.get_cmap('coolwarm');inputs={}
    for index,(name,label) in enumerate((('flat','Box excavation'),('curved','Cylindrical excavation'))):
        mp=run/name/'ground/mesh.json';rp=run/name/'stage-results.json'
        mesh=json.loads(mp.read_text());data=json.loads(rp.read_text());state=data['last_committed']
        xyz=np.asarray(mesh['nodes']);uz=1000*np.asarray(state['displacement_m'])[:,2]
        active=np.asarray(state['active_elements']);owners=ref.topology(mesh);faces=[];values=[]
        for face,adjacent in owners.items():
            live=[e for e in adjacent if active[e]]
            if len(live)!=1:continue
            nodes=ref.face_nodes(mesh,live[0],face);corners=xyz[nodes[:3]]
            center=corners.mean(axis=0);normal=np.cross(corners[1]-corners[0],corners[2]-corners[0])
            # Original top, excavation floor and walls; omit outer sides/base.
            inner=not any(abs(center[a]-edge)<1e-8 for a in (0,1) for edge in (0,10))
            if normal[2]<=1e-12 and not (inner and len(adjacent)==2):continue
            for local in ((0,3,5),(3,1,4),(5,4,2),(3,4,5)) if len(nodes)==6 else ((0,1,2),):
                ids=np.asarray(nodes)[list(local)];faces.append(xyz[ids]);values.append(uz[ids].mean())
        ax=fig.add_subplot(1,2,index+1,projection='3d')
        surface=Poly3DCollection(faces,facecolors=cmap(norm(values)),edgecolors=(0,0,0,.12),linewidths=.25)
        ax.add_collection3d(surface);ax.set(xlim=(0,10),ylim=(0,10),zlim=(5.5,8.5),xlabel='X [m]',ylabel='Y [m]',title=label)
        ax.text2D(.85,.72,'Z [m]',transform=ax.transAxes,fontsize=10)
        ax.set_box_aspect((10,10,4));ax.view_init(elev=35,azim=-55)
        m=data['frames'][-1]['monitors'];floor=m['excavation_floor']['mean_displacement_z_m']*1000;settle=m['surrounding_surface']['minimum_displacement_z_m']*1000
        ax.set_title(f'{label}\nFloor mean: +{floor:.3f} mm | Surface min: {settle:.3f} mm',fontsize=11,pad=8)
        for path in (mp,rp):inputs[str(path.resolve())]=hashlib.sha256(path.read_bytes()).hexdigest()
    cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),ax=fig.axes,shrink=.7,pad=.03,label='Vertical displacement uz [mm] (+ heave / - settlement)')
    cb.set_ticks([minimum,0,maximum/2,maximum])
    fig.suptitle('SolidCraft geometry → Rust 3D MCC excavation\nUndeformed geometry; nodal displacement colour; fixed pore pressure',fontsize=13)
    for suffix in ('png','svg'):fig.savefig(output/f'local-excavation.{suffix}',dpi=180)
    plt.close(fig)
    shutil.copyfile(Path(__file__),output/'plot-source.py')
    inputs[str(Path(__file__).resolve())]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    report={'schema':'soil-fem-local-excavation-figure/1','parent_integrity':integrity,
        'numpy_version':np.__version__,'matplotlib_version':matplotlib.__version__,
        'inputs_sha256':inputs,'colour_limits_mm':[minimum,maximum],
        'display':'Undeformed top/floor/wall geometry; four display triangles per Tri6, colour mean of raw nodal uz. No field smoothing or colour clipping. Floor mean is nodal, not area weighted.'}
    with (output/'provenance.json').open('x',encoding='utf-8') as stream:json.dump(report,stream,indent=2)
    print(output.resolve())


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();plot(a.run.resolve(),a.output.resolve())
