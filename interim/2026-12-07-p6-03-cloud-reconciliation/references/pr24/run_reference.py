"""Generate small independent material references and original MC column diagnosis."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import platform
import time
import numpy as np
import scipy
from mcc_reference import Material, integrate, isotropic, mc_column, richardson
from compare_history import META, compare
from column_reference import diagnose, integrated_stages

BASE = 'f377e8eda41777c84c3bb849af556bf61683b4cb'
SOURCES = {'S1': 'docs/MCC_RATE_FEM.ja.md', 'S2': 'docs/MCC_TRIAXIAL_VERIFICATION.ja.md',
           'S3': 'docs/MCC_AFFINE_SEED.ja.md', 'S4': 'docs/MCC_MIDPOINT.ja.md',
           'S5': 'docs/PLANE_STRAIN_2D.ja.md', 'S6': 'docs/2D_COMPARISON.ja.md',
           'S7': 'interim/2026-11-16-issue-triage/todo.ja.md',
           'S8': 'docs/MCC_MIXED.ja.md', 'S9': 'tools/prepare_mixed_mcc_reference.py',
           'S10': 'examples/ground-mixed-mcc.json'}
FIELDS = ('p', 'q', 'pc', 'e', 'eps_v')


def write(path, value):
    with Path(path).open('x', encoding='utf-8') as f:
        json.dump(value, f, indent=2, ensure_ascii=False, allow_nan=False)
        f.write('\n')


def history(rows, mode, ocr, material, source_hash):
    physics = dict(mode=mode, ocr=ocr, material=vars(material), p0=100.,
                   end=rows[-1]['axial_strain'], strain_measure='log_volume')
    digest = hashlib.sha256(json.dumps(physics, sort_keys=True).encode()).hexdigest()
    return dict(metadata=dict(META, case_id=f'{mode}-ocr{ocr:g}', kind='gp',
                              physics_sha256=digest, source_sha256=source_hash,
                              ids=['gp0'], load_coordinates=[r['axial_strain'] for r in rows]),
                rows=[dict(step=r['step'], id='gp0', load_coordinate=r['axial_strain'],
                           **{k:r[k] for k in FIELDS}) for r in rows])


def generate(output, repo, levels=(128, 256, 512), column_material_levels=(256,512,1024),
             column_spatial_cells=(8,16,32), column_spatial_material_n=2048):
    started = time.monotonic()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    if len(levels) != 3 or levels[1] != 2*levels[0] or levels[2] != 2*levels[1] or levels[0] < 1:
        raise ValueError('three positive doubled n-per-output levels required')
    here = Path(__file__).resolve().parent
    hashes = {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(here.glob('*.py'))}
    source_hash = hashes['mcc_reference.py']
    sources = {ident:dict(path=path, sha256=hashlib.sha256((repo/path).read_bytes()).hexdigest())
               for ident, path in SOURCES.items()}
    mat = Material()
    report = dict(base_commit=BASE, status='diagnostic_only', accepted=None, tolerance='未定義',
                  p6_03_complete=False, p6_complete=False, nc_01_complete=False,
                  source_documents=sources, implementation_sha256=hashes,
                  environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
                  parameters=dict(vars(mat), p0=100., pressure_unit='kPa',
                                  source_ids=['S1','S2','S3']), material_paths={})
    iso = isotropic([100.,130.,100.])
    report['isotropic'] = dict(rows=iso, solution='closed_form; no load-step truncation',
                              observed_order=None, accepted=None, tolerance='未定義')
    write(output/'isotropic.json', report['isotropic'])
    for mode in ('k0','undrained','drained'):
        for ocr in (1.,4.):
            name = f'{mode}-ocr{ocr:g}'
            all_rows = []
            for n in levels:
                rows, event = integrate(mode, n, ocr=ocr, mat=mat)
                all_rows.append(rows)
                write(output/f'{name}-n{n}.json', dict(mode=mode, ocr=ocr, event_axial_strain=event,
                      substeps_per_output_interval=n, output_intervals=16, rows=rows))
            rows = all_rows[-1]
            comparisons = {k:richardson([[r[k] for r in rr] for rr in all_rows]) for k in FIELDS}
            # S3 invariant. Signed f retained; never project onto the yield surface.
            drifts = []
            for rr in all_rows:
                plastic = [r for r in rr if event is not None and r['axial_strain'] >= event]
                drifts.append(dict(
                    max_normalized_yield_abs=max((abs((r['q']/r['pc'])**2+mat.M**2*r['p']/r['pc']*(r['p']/r['pc']-1)) for r in plastic), default=0.),
                    max_ncl_volume_abs=max(abs(1+r['e']-(mat.N-mat.lam*math.log(r['pc']/mat.pref)+mat.kap*math.log(r['pc']/r['p']))) for r in rr),
                    max_volume_strain_identity_abs=max(abs(math.log((1+rr[0]['e'])/(1+r['e']))-r['eps_v']) for r in rr)))
            item = dict(accepted=None, tolerance='未定義', source_ids=['S1','S2','S3'],
                        initial=rows[0], final=rows[-1], event_axial_strain=event,
                        axial_strain_end=.08, substeps_per_output_interval=list(levels),
                        increment_series=comparisons, invariant_drift_by_level=drifts,
                        pc_regression=any(b['pc']<a['pc'] for a,b in zip(rows,rows[1:])),
                        spatial_order=None, spatial_reason='material point; no spatial discretization')
            if mode in ('undrained','drained'):
                pcs = 100*(ocr/2)**((mat.lam-mat.kap)/mat.lam) if mode=='undrained' else 100/(1-mat.M/3)
                item['critical_state_limit'] = dict(p=pcs, q=mat.M*pcs, pc=2*pcs,
                    source='S2; asymptotic limit, NOT finite axial strain endpoint reference')
            report['material_paths'][name] = item
            h = history(rows,mode,ocr,mat,source_hash)
            if item['pc_regression']:
                # Do not lie about the exchange contract; these files are raw model results only.
                item['exchange_contract'] = 'REJECTED: physical OC dilative pc decrease; monotone-pc exchange unsupported'
            else:
                write(output/f'{name}-gp-reference.json', h)
                with (output/f'{name}-gp-reference.csv').open('x',encoding='utf-8',newline='') as f:
                    writer=csv.DictWriter(f,fieldnames=['step','id','load_coordinate',*FIELDS])
                    writer.writeheader();writer.writerows(h['rows'])
                write(output/f'{name}-gp-reference.metadata.json',h['metadata'])
                # Deliberately self comparison: format/validator smoke only, never body evidence.
                write(output/f'{name}-format-self-comparison.json',compare(h,h))
    # Original column is MC, including p'=0 at the free top in gravity. Never initialize MCC there.
    z = np.linspace(0,.2,65)
    column = {stage:mc_column(z,stage) for stage in ('gravity','load','unload')}
    physics_sha=hashlib.sha256(b'S5 S6 original MC column .2 x .2 x 1 E1000 nu.25 gamma20 c2 phi10 psi0 P12').hexdigest()
    node_meta=dict(META,kind='node',case_id='original-mc-column',physics_sha256=physics_sha,
                   source_sha256=source_hash,ids=[f'n{i}' for i in range(len(z))],
                   load_coordinates=[0.,1.,2.,3.],
                   node_coordinates_m={f'n{i}':[.1,.5,float(zz)] for i,zz in enumerate(z)})
    node_rows=[]
    for step,stage in enumerate((None,'gravity','load','unload')):
        for i in range(len(z)):
            node_rows.append(dict(step=step,id=f'n{i}',load_coordinate=float(step),ux=0.,uy=0.,
                                  uz=0. if stage is None else column[stage]['uz'][i]))
    write(output/'original-mc-column-node-reference.json',dict(metadata=node_meta,rows=node_rows))
    write(output/'original-mc-column.json',dict(source_ids=['S5','S6'],material='nonassociated_MC',stages=column))
    # Exact quadratic interpolation error at 8 interior Gauss points in EACH fine cell.
    # This is 1D P1 field representation, not a new Tet4 solve or MCC spatial certification.
    gx,gw=np.polynomial.legendre.leggauss(8)
    common=np.concatenate([(.2/32)*(i+(gx+1)/2) for i in range(32)])
    w=np.tile(gw*.2/(2*32),32)
    exact=np.array(mc_column(common,'gravity')['uz'])
    representations=[]; errors=[]
    for n in (4,8,16):
        nodes=np.linspace(0,.2,n+1)
        value=np.interp(common,nodes,mc_column(nodes,'gravity')['uz'])
        representations.append(value*np.sqrt(w))
        errors.append(float(np.linalg.norm((value-exact)*np.sqrt(w))))
    order=richardson(representations)
    order['fine_vs_exact_l2_m_sqrt_m']=errors[-1]
    order['representation_error_orders']=[math.log(a/b,2) for a,b in zip(errors,errors[1:])]
    order['description']='independent 1D MC gravity P1 interpolation only; NOT Tet4 or MCC FE results'
    order['h_m']=[.05,.025,.0125]
    # Samples/weighted extrapolation are not displacement node histories.
    write(output/'one-dimensional-spatial-diagnostic.json',order)
    report['original_column']=dict(material='MC; MCC mapping 未定義', source_ids=['S5','S6'],
        top_displacement_m={stage:column[stage]['uz'][-1] for stage in column},
        independent_1d_representation_order=order['observed_order'],
        original_tet4_fine_relative_l2=.068282,original_tet4_fine_accepted=False,
        original_continuum_relative_l2=.010336,original_continuum_accepted=True,
        original_reported_order=1.197,original_threshold_relative_l2=.05,
        original_discretization_component=.057946,
        original_stress_extrusion_accepted=None,original_finite_width_stress_accepted=None,
        mcc_column_same_conditions='未確認: original .2m MC column has no documented MCC mapping; separate 1m MCC column S8--S10 identified',
        mcc_general_spatial_order=None)
    for seq in (column_material_levels,column_spatial_cells):
        if len(seq)!=3 or seq[0]<1 or seq[1]!=2*seq[0] or seq[2]!=2*seq[1]:
            raise ValueError('column series needs three positive doubled levels')
    if column_spatial_material_n<1:
        raise ValueError('column material n must be positive')
    mcc_column=diagnose(column_material_levels,column_spatial_cells,column_spatial_material_n)
    write(output/'documented-mcc-column.json',mcc_column)
    report['documented_mcc_column']=mcc_column
    # Exchange only the common surviving lower domain; upper removed histories
    # remain explicit raw point histories in documented-mcc-column.json.
    col_physics=hashlib.sha256(b'S8-S10 fixed 1m MCC column lower-domain histories').hexdigest()
    col_source=hashlib.sha256((hashes['mcc_reference.py']+hashes['column_reference.py']).encode()).hexdigest()
    gp_points=mcc_column['point_histories'][:2]
    col_meta=dict(META,kind='gp',case_id='documented-mcc-column-lower-domain',
        physics_sha256=col_physics,source_sha256=col_source,
        ids=['gp0','gp1'],load_coordinates=[0.,1.,2.,3.],
        gp_coordinates_m={'gp0':[.5,.5,.125],'gp1':[.5,.5,.375]},
        scope='only common surviving LowerClay domain; no upper history interpolation')
    col_rows=[]
    for step in range(4):
        for i,pt in enumerate(gp_points):
            r=pt['rows'][step]
            col_rows.append(dict(step=step,id=f'gp{i}',load_coordinate=float(step),**{k:r[k] for k in FIELDS}))
    col_gp=dict(metadata=col_meta,rows=col_rows)
    compare(col_gp,col_gp)  # strict validation only
    write(output/'documented-mcc-column-gp-reference.json',col_gp)
    col_node_meta=dict(col_meta,kind='node',ids=['base','middle','floor'])
    col_node_meta.pop('gp_coordinates_m')
    col_node_meta['node_coordinates_m']={'base':[.5,.5,0.],'middle':[.5,.5,.25],'floor':[.5,.5,.5]}
    nodes=[('base',[0.,0.,0.]),('middle',integrated_stages(.25,n=column_spatial_material_n,quadrature=16)),
           ('floor',mcc_column['continuum_quadrature_values_m'])]
    col_node=dict(metadata=col_node_meta,rows=[dict(step=step,id=ident,load_coordinate=float(step),ux=0.,uy=0.,
        uz=0. if step==0 else vals[step-1]) for step in range(4) for ident,vals in nodes])
    compare(col_node,col_node)
    write(output/'documented-mcc-column-node-reference.json',col_node)

    common=dict(physics_sha256=col_physics,observable='lower_top_uz',unit='m',
                sample_ids=['load','unload','cut'],depth_m=1.,boundary_id='S9-all-horizontal-fixed-bottom-z-fixed',
                load_path_id='S8-100-130-100-cut',load_increment=30.,
                mesh_sha256=hashlib.sha256(b'1D continuum GL16 fixed rule').hexdigest())
    material_series=dict(axis='material_increment',parameter_unit='kPa',run_status='completed',
        levels=[dict(common,material_increment=30/n,values=values)
            for n,values in zip(column_material_levels,mcc_column['material_increment_values_m'])])
    spatial_series=dict(axis='spatial',parameter_unit='m',run_status='completed',nested_h_refinement=True,
        levels=[dict(common,h_m=.5/c,material_increment=30/column_spatial_material_n,values=values,
                    mesh_sha256=hashlib.sha256(f'1D midpoint cells {c}'.encode()).hexdigest())
            for c,values in zip(column_spatial_cells,mcc_column['spatial_values_m'])])
    write(output/'column-material-series-input.json',material_series)
    write(output/'column-spatial-series-input.json',spatial_series)
    report['elapsed_seconds']=time.monotonic()-started
    write(output/'summary.json',report)
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[2])
    parser.add_argument('--levels',nargs=3,type=int,default=[128,256,512])
    parser.add_argument('--column-material-levels',nargs=3,type=int,default=[256,512,1024])
    parser.add_argument('--column-spatial-cells',nargs=3,type=int,default=[8,16,32])
    parser.add_argument('--column-spatial-material-n',type=int,default=2048)
    args=parser.parse_args()
    result=generate(args.output,args.repo,args.levels,args.column_material_levels,
                    args.column_spatial_cells,args.column_spatial_material_n)
    print(json.dumps(dict(status=result['status'],accepted=None,elapsed_seconds=result['elapsed_seconds']),ensure_ascii=False))


if __name__=='__main__':
    main()
