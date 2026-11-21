"""Small dry-MC excavation h series through the unchanged production CLI.

Tet4 stress is elementwise constant. Intersection-volume norms therefore compare
the complete displacement and stress fields without GP history interpolation.
The default refines x/z at fixed depth. A separate depth series fixes x/z.
The physical domain and cut stay fixed in both series.
"""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
import numpy as np
from tet10_overlap import Overlay

ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT/'target/release'/('soil-fem-3d.exe' if os.name == 'nt' else 'soil-fem-3d')
EDGES = ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3))


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def write(path, data):
    # Canonical LF keeps input fingerprints reproducible on Windows as well.
    with Path(path).open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(data, stream, indent=2, allow_nan=False)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def richardson(values):
    """Signed monotone scalar sequence, no external reference or fitted truth."""
    values = np.asarray(values, float)
    result = dict(values=values.tolist(), order=None, estimate=None,
                  finest_discretization_component=None, reason='insufficient_levels')
    if not np.isfinite(values).all():
        raise ValueError('nonfinite convergence values')
    if len(values) < 3:
        return result
    a, b = np.diff(values[-3:])
    if max(abs(a),abs(b)) <= 1e-8*max(np.max(abs(values)),np.finfo(float).tiny):
        result['reason'] = 'changes_at_force_tolerance_scale_no_order'
        return result
    if a*b <= 0:
        result['reason'] = 'zero_or_oscillating_differences'
        return result
    if abs(b) >= .7*abs(a):
        result['reason'] = 'differences_do_not_decrease_by_0.7'
        return result
    order = float(np.log2(abs(a/b)))
    if not .5 <= order <= 8:
        result['reason'] = 'order_outside_existing_0.5_to_8_guard'
        return result
    correction = float(b/(2**order-1))
    result.update(order=order, estimate=float(values[-1]+correction),
                  finest_discretization_component=abs(correction), reason='conditional_asymptotic_estimate')
    return result


def validate_divisions(n, ny, finite_width):
    if n < 2 or n > 64 or n & (n-1):
        raise ValueError('in-plane divisions must be a power of two in 2..64')
    if ny < 1 or ny > 64 or ny & (ny-1) or (finite_width and ny < 2):
        raise ValueError('depth divisions must be a power of two in 1..64 (finite width requires >=2)')


def inputs(n, finite_width=False, increment=1/16, material_substep=1e-4, depth_divisions=None):
    ny = depth_divisions if depth_divisions is not None else (2 if finite_width else 1)
    validate_divisions(n, ny, finite_width)
    nodes = [[i/n, j/ny, k/n] for k in range(n+1) for j in range(ny+1) for i in range(n+1)]
    index = lambda i, j, k: i+(n+1)*(j+(ny+1)*k)
    cells, regions = [], []
    for k in range(n):
        for j in range(ny):
            for i in range(n):
                c = [index(i,j,k), index(i+1,j,k), index(i+1,j+1,k), index(i,j+1,k),
                     index(i,j,k+1), index(i+1,j,k+1), index(i+1,j+1,k+1), index(i,j+1,k+1)]
                cut = i < n//2 and k >= n//2 and (not finite_width or j < ny//2)
                for local in ((0,1,2,6),(0,2,3,6),(0,3,7,6),(0,7,4,6),(0,4,5,6),(0,5,1,6)):
                    tet = [c[q] for q in local]
                    xyz = np.asarray([nodes[q] for q in tet])
                    if np.linalg.det((xyz[1:]-xyz[0]).T) < 0:
                        tet[1], tet[2] = tet[2], tet[1]
                    cells.append(tet)
                    regions.append(2 if cut else 1)
    xyz = np.asarray(nodes)
    sets = dict(all=list(range(len(nodes))), bottom=np.flatnonzero(xyz[:,2] == 0).tolist(),
                sides=np.flatnonzero((xyz[:,0] == 0) | (xyz[:,0] == 1)).tolist(),
                ends=np.flatnonzero((xyz[:,1] == 0) | (xyz[:,1] == 1)).tolist())
    cut_volume = .125 if finite_width else .25
    mesh = dict(schema='soil-fem-mesh/2', units='m-kN-kPa', axis='z_up', crs='LOCAL', origin_m=[0.,0.,0.],
                nodes=nodes, tetrahedra=cells, material_ids=[1]*len(cells), region_ids=regions,
                node_sets=sets, face_sets={}, excavation_sets={}, regions=[
                    dict(id=q, body=name, material_id=1, bounds_m=[[0.,0.,0.],[1.,1.,1.]],
                         volume_m3=volume, source_sha256=hashlib.sha256(b'dry-MC-unit-notch-v1').hexdigest())
                    for q,name,volume in ((1,'Soil',1-cut_volume),(2,'Cut',cut_volume))])
    zero = dict(effective_stress_compression_positive_kpa=[0.]*6, plastic_strain_compression_positive=[0.]*6,
                void_ratio=.7, equivalent_plastic_strain=0., plastic_work_kpa=0., last_return='elastic')
    ground = dict(schema='soil-fem-ground/1', units='m-kN-kPa',
                  initialization=dict(method='equilibrate', surface_surcharge_kpa=0., water_level_m=None,
                      water_unit_weight_kn_m3=9.80665, relative_equilibrium_tolerance=1e-8,
                      geostatic=dict(seed=dict(kind='body_states', states={'Soil':zero,'Cut':zero}),
                          initial_boundary_factor=0., maximum_initial_engineering_strain=.05,
                          relative_volume_tolerance=1e-8)),
                  materials=[dict(id=1, name='dry nonassociated MC notch', unit_weight_kn_m3=20.,
                      mohr_coulomb=dict(young_kpa=10000., poisson=.25, cohesion_kpa=2., friction_deg=20.,
                          dilation_deg=0., tension_cutoff_kpa=None, flow='non_associated'))],
                  bodies=[dict(body=name,material_id=1) for name in ('Soil','Cut')], layers=[])
    case = read(ROOT/'examples/stages-mc-column.json')
    case['settings'].update(kinematics='explicit', relative_force_tolerance=1e-8,
        relative_linear_tolerance=1e-10, maximum_newton_iterations=35, initial_increment=increment,
        minimum_increment=1/4096, material_maximum_substep=material_substep)
    case['boundary'] = dict(constraints=[dict(node_set=name, direction=np.eye(3)[axis].tolist(),value_m=0.)
        for name,axis in [('sides',0),('ends' if finite_width else 'all',1),('bottom',2)]])
    case['monitors'] = {}; case['monitor_sets'] = {}; case['verification_reject_once_stages'] = []
    case['stages'] = [dict(kind='excavate',name='cut',bodies=['Cut'],remaining_surface_pressure_kpa=0.,
                          initial_increment=increment)]
    return mesh, ground, case


def lift(mesh, displacement):
    """Exact representation of a linear field for the existing Tet10 overlay."""
    result = copy.deepcopy(mesh)
    nodes = list(mesh['nodes']); u = list(displacement); edges = {}; cells = []
    for tet in mesh['tetrahedra']:
        mids = []
        for a,b in EDGES:
            key = tuple(sorted((tet[a],tet[b])))
            if key not in edges:
                edges[key] = len(nodes)
                nodes.append(((np.asarray(nodes[key[0]])+nodes[key[1]])/2).tolist())
                u.append(((np.asarray(u[key[0]])+u[key[1]])/2).tolist())
            mids.append(edges[key])
        cells.append(tet+mids)
    result.update(schema='soil-fem-mesh/3', nodes=nodes, tetrahedra10=cells, quadrature_order=3)
    return result, u


def compare(mesh_a, run_a, mesh_b, run_b, diagnostics=False):
    ma, _ = lift(mesh_a, run_a['frames'][0]['state']['displacement_m'])
    mb, _ = lift(mesh_b, run_b['frames'][0]['state']['displacement_m'])
    overlay = Overlay(ma, mb)
    reports = []
    for fa,fb in zip(run_a['frames'],run_b['frames'],strict=True):
        if fa['name'] != fb['name']:
            raise ValueError('frame sequence differs')
        a,b = fa['state'],fb['state']
        _,ua = lift(mesh_a,a['displacement_m']); _,ub = lift(mesh_b,b['displacement_m'])
        active_a,active_b = a['active_elements'],b['active_elements']
        selected,coverage = overlay.coverage(active_a,active_b)
        # Tet4 has exactly one material point per cell: do not silently accept Tet10 GP arrays.
        if len(a['material_states']) != len(mesh_a['tetrahedra']) or len(b['material_states']) != len(mesh_b['tetrahedra']):
            raise ValueError('element-constant Tet4 stress required')
        sa = np.asarray([s['effective_stress_compression_positive_kpa'] for s in a['material_states']])
        sb = np.asarray([s['effective_stress_compression_positive_kpa'] for s in b['material_states']])
        delta = sa[overlay.ea[selected]]-sb[overlay.eb[selected]]
        weights = overlay.volume[selected]
        # Tensor Frobenius norm: engineering stress storage has three distinct shears.
        metric = np.array([1,1,1,2,2,2])
        error = float(np.sum(weights*np.sum(delta**2*metric,axis=1)))
        ref = float(np.sum(weights*np.sum(sb[overlay.eb[selected]]**2*metric,axis=1)))
        displacement = overlay.compare(ua,ub,active_a,active_b)
        initial_a = np.asarray(run_a['frames'][0]['state']['displacement_m'])
        initial_b = np.asarray(run_b['frames'][0]['state']['displacement_m'])
        _,da = lift(mesh_a,(np.asarray(a['displacement_m'])-initial_a).tolist())
        _,db = lift(mesh_b,(np.asarray(b['displacement_m'])-initial_b).tolist())
        reports.append(dict(stage=fa['name'], displacement=displacement,
            construction_displacement=overlay.compare(da,db,active_a,active_b),
            stress=dict(**coverage,rms_difference_kpa=float(np.sqrt(error/coverage['active_volume_m3'])),
                        relative_l2=float(np.sqrt(error/ref)) if ref > 0 else None)))
        if diagnostics:
            components = np.sum(weights[:,None]*delta**2*metric,axis=0)
            pressure = np.mean(delta[:,:3],axis=1)
            pressure_error = float(np.sum(weights*3*pressure**2))
            deviator = delta.copy(); deviator[:,:3] -= pressure[:,None]
            deviator_error = float(np.sum(weights*np.sum(deviator**2*metric,axis=1)))
            pa = np.asarray([s['equivalent_plastic_strain'] > 0 for s in a['material_states']])[overlay.ea[selected]]
            pb = np.asarray([s['equivalent_plastic_strain'] > 0 for s in b['material_states']])[overlay.eb[selected]]
            density = np.sum(delta**2*metric,axis=1)
            reports[-1]['stress_diagnostic'] = dict(
                component_order=['xx','yy','zz','xy','yz','xz'],
                tensor_weighted_component_rms_kpa=np.sqrt(components/coverage['active_volume_m3']).tolist(),
                component_squared_integrals_kpa2_m3=components.tolist(),
                hydrostatic_squared_integral_kpa2_m3=pressure_error,
                deviatoric_squared_integral_kpa2_m3=deviator_error,
                plastic_partition=[dict(region=name,volume_m3=float(weights[mask].sum()),
                    squared_integral_kpa2_m3=float(np.sum(weights[mask]*density[mask])))
                    for name,mask in [('both_plastic',pa & pb),('one_plastic',pa ^ pb),('both_elastic',~pa & ~pb)]])
    return dict(geometry=overlay.metadata, stages=reports)


def summarize(mesh, run):
    frames = []
    for f in run['frames']:
        state = f['state']; volume=0.; mean=np.zeros(6); plastic_volume=0.
        for point,material in zip(run['integration_points'],state['material_states'],strict=True):
            if state['active_elements'][point['element_index']]:
                w=point['weight_m3']; volume+=w
                mean+=w*np.asarray(material['effective_stress_compression_positive_kpa'])
                if material['equivalent_plastic_strain'] > 0:
                    plastic_volume+=w
        # Same retained physical node on every mesh; never nearest-node sampling.
        node = mesh['nodes'].index([.5,0.,.5])
        frames.append(dict(stage=f['name'], monitor_uz_m=state['displacement_m'][node][2],
            construction_monitor_uz_m=state['displacement_m'][node][2]-run['frames'][0]['state']['displacement_m'][node][2],
            active_volume_m3=volume, active_plastic_volume_m3=plastic_volume,
            mean_stress_kpa=(mean/volume).tolist(), relative_force_residual=f['relative_force_residual'],
            maximum_normalized_yield=f['maximum_normalized_yield'],
            maximum_engineering_strain=f['maximum_absolute_engineering_strain_from_seed'],
            global_force_relative_error=float(np.linalg.norm(f['force_balance_error_kn'])/
                max(np.linalg.norm(f['applied_force_kn']),np.finfo(float).tiny))))
    return dict(complete=run['complete'],failure_reason=run['failure_reason'],nodes=len(mesh['nodes']),
                elements=len(mesh['tetrahedra']),frames=frames)


def run_series(folder, binary, divisions=(2,4,8), finite_width=False, increment=1/16, timeout=120,
               material_substep=1e-4, axis='in_plane', depth_divisions=None, in_plane_divisions=2):
    if axis not in ('in_plane','depth'):
        raise ValueError('axis must be in_plane or depth')
    if len(divisions)<3 or any(b!=2*a for a,b in zip(divisions,divisions[1:])):
        raise ValueError('at least three levels with ratio two required')
    if axis == 'depth' and depth_divisions is not None:
        raise ValueError('fixed depth divisions cannot be combined with a depth series')
    # Validate every level before creating an output directory or starting any solve.
    for n in divisions:
        validate_divisions(n if axis == 'in_plane' else in_plane_divisions,
                           (depth_divisions if depth_divisions is not None else (2 if finite_width else 1))
                           if axis == 'in_plane' else n, finite_width)
    if not 0 < increment <= 1 or not 0 < material_substep <= 1e-4 or not np.isfinite(timeout) or timeout <= 0:
        raise ValueError('positive timeout, increment in (0,1] and material substep in (0,1e-4] required')
    folder = Path(folder); folder.mkdir(parents=True,exist_ok=False)
    levels=[]; solved=[]; failures=[]
    for n in divisions:
        sub=folder/(f'n{n}' if axis == 'in_plane' else f'ny{n}'); sub.mkdir()
        nx = n if axis == 'in_plane' else in_plane_divisions
        ny = (depth_divisions if depth_divisions is not None else (2 if finite_width else 1)) if axis == 'in_plane' else n
        mesh,ground,case=inputs(nx,finite_width,increment,material_substep,ny)
        for name,value in [('mesh',mesh),('ground',ground),('case',case)]:
            write(sub/f'{name}.json',value)
        command=[str(Path(binary).resolve()),'solve-stages',str(sub/'mesh.json'),str(sub/'ground.json'),
                 str(sub/'case.json'),str(sub/'result.json')]
        started=time.monotonic()
        try:
            process=subprocess.run(command,capture_output=True,text=True,timeout=timeout,check=False)
        except subprocess.TimeoutExpired:
            write(sub/'failure.json',dict(reason='wall_time_limit',timeout_seconds=timeout))
            failures.append(dict(divisions=n,reason='wall_time_limit',timeout_seconds=timeout))
            break
        (sub/'solver.log').write_text(process.stdout+process.stderr,encoding='utf-8')
        if not (sub/'result.json').exists():
            write(sub/'failure.json',dict(reason='no_result',exit_code=process.returncode))
            failures.append(dict(divisions=n,reason='no_result',exit_code=process.returncode))
            break
        result=read(sub/'result.json')
        summary=summarize(mesh,result)
        summary.update(divisions=n,h_xz_m=1/nx,h_y_m=1/ny,in_plane_divisions=nx,depth_divisions=ny,
                       solver_elapsed_seconds=time.monotonic()-started,
                       files_sha256={name:sha(sub/name) for name in ('mesh.json','ground.json','case.json','result.json')})
        levels.append(summary)
        write(sub/'summary.json',summary)
        if not result['complete']:
            failures.append(dict(divisions=n,reason=result['failure_reason']))
            break
        if process.returncode != 0 or len(summary['frames']) != 2 or any(
            f['relative_force_residual']>1e-8 or f['global_force_relative_error']>1e-8
            or f['maximum_normalized_yield']>1e-7 or f['maximum_engineering_strain'] is None
            or f['maximum_engineering_strain']>.05 for f in summary['frames']):
            failures.append(dict(divisions=n,reason='original_physical_audit_gate_failed'))
            break
        solved.append((mesh,result))
    pairs=[compare(*a,*b,diagnostics=True) for a,b in zip(solved,solved[1:])]
    quantities=[]
    if solved:
        for stage_index,frame in enumerate(levels[0]['frames']):
            for quantity,extract in [('monitor_uz_m',lambda f:f['monitor_uz_m']),
                                    ('construction_monitor_uz_m',lambda f:f['construction_monitor_uz_m']),
                                    ('mean_sigma_xx_kpa',lambda f:f['mean_stress_kpa'][0]),
                                    ('mean_sigma_zz_kpa',lambda f:f['mean_stress_kpa'][2])]:
                values=[extract(level['frames'][stage_index]) for level in levels[:len(solved)]]
                quantities.append(dict(stage=frame['stage'],quantity=quantity,**richardson(values)))
    # Norm-ratio orders are recorded, but a vector field has no scalar signed Richardson value.
    norm_orders=[]
    if len(pairs)>=2:
        for stage in range(len(pairs[-1]['stages'])):
            for field,key in [('displacement','rms_difference_m'),('construction_displacement','rms_difference_m'),
                              ('stress','rms_difference_kpa')]:
                a,b=[p['stages'][stage][field][key] for p in pairs[-2:]]
                order=float(np.log2(a/b)) if a>0 and b>0 and b<.7*a else None
                observable=order is not None and .5<=order<=8
                norm_orders.append(dict(stage=pairs[-1]['stages'][stage]['stage'],field=field,
                    difference_norms=[a,b],order=order if observable else None,
                    finest_error_norm_estimate=b/(2**order-1) if observable else None,
                    conditional_estimate=True, reason='conditional_same_leading_error_field' if observable else 'indeterminate'))
    report=dict(schema='soil-fem-plastic-spatial-convergence/1',p6_complete=False,
        general_3d_spatial_convergence_demonstrated=False,griffiths_slope_five_percent_closed=False,
        refinement_axis=axis,finite_width=finite_width,requested_divisions=list(divisions),levels=levels,pairs=pairs,
        scalar_richardson=quantities,field_difference_orders=norm_orders,
        solver_sha256=sha(binary),driver_sha256=sha(__file__),
        complete=len(solved)==len(divisions),failures=failures, fs=dict(status='not_attempted',
            reason='no independently excluded upper factor for this nonuniform cut; nonconvergence is not collapse'),
        limitations=[('Only x/z are refined; depth discretization is fixed and needs a separate study.' if axis == 'in_plane'
                      else 'Only depth is refined; x/z discretization is fixed and needs a separate study.'),
            'Three levels give conditional orders, not a proof of the asymptotic regime.',
            'Load and material increments are fixed; run the half-increment control separately.',
            'The original Griffiths 661-node 2D displacement data is absent; its 5% comparison stays open.',
            'Existing MCC excavations, Griffiths slope, water, structures, JOINT and P6 remain unaccepted.'])
    write(folder/'report.json',report)
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    parser.add_argument('--binary',type=Path,default=BINARY)
    parser.add_argument('--divisions',type=int,nargs='+',default=[2,4,8])
    parser.add_argument('--finite-width',action='store_true')
    parser.add_argument('--axis',choices=['in_plane','depth'],default='in_plane')
    parser.add_argument('--depth-divisions',type=int,help='fixed depth for an in-plane series; use a new output')
    parser.add_argument('--in-plane-divisions',type=int,default=2,help='fixed x/z for a depth series')
    parser.add_argument('--increment',type=float,default=1/16)
    parser.add_argument('--material-substep',type=float,default=1e-4)
    parser.add_argument('--control-baseline',type=Path,help='compare each level with an existing same-mesh baseline')
    parser.add_argument('--timeout',type=float,default=120,help='per solve; workstation users may explicitly increase')
    args=parser.parse_args()
    if len(args.divisions)<3 or any(b!=2*a for a,b in zip(args.divisions,args.divisions[1:])):
        parser.error('at least three levels with ratio two required')
    report=run_series(args.output.resolve(),args.binary,args.divisions,args.finite_width,args.increment,args.timeout,
                      args.material_substep,args.axis,args.depth_divisions,args.in_plane_divisions)
    if args.control_baseline and report['complete']:
        comparisons=[]
        for n in args.divisions:
            name=f'n{n}' if args.axis == 'in_plane' else f'ny{n}'
            a=args.control_baseline/name; b=args.output/name
            if read(a/'mesh.json')!=read(b/'mesh.json') or read(a/'ground.json')!=read(b/'ground.json'):
                raise ValueError('control requires identical mesh and ground')
            comparisons.append(dict(divisions=n,**compare(read(a/'mesh.json'),read(a/'result.json'),
                                                         read(b/'mesh.json'),read(b/'result.json'))))
        write(args.output/'control-comparison.json',comparisons)
    print(json.dumps(dict(complete=report['complete'],report=str(args.output/'report.json'),p6_complete=False)))
    return 0 if report['complete'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
