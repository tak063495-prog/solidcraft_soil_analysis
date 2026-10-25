"""NumPy isoparametric basis/integration used only for independent verification.

The independent rule uses tensor Gauss-Legendre with a Duffy transform,
whereas Rust uses a Gauss-Jacobi product rule.
"""
import numpy as np

EDGES=((0,1),(1,2),(2,0),(0,3),(1,3),(2,3))
FACES=((1,2,3),(0,3,2),(0,1,3),(0,2,1))


def normalize_conditions(value):
    value.setdefault('nodal_loads',[]);value.setdefault('surface_loads',[])
    for c in value['constraints']:
        c.setdefault('gradient',[0.,0.,0.]);c.setdefault('reference_point_m',[0.,0.,0.])
    for load in value['surface_loads']:
        load.setdefault('gradient_kpa_per_m',[[0.]*3 for _ in range(3)] if load['kind']=='traction' else [0.]*3)
        load.setdefault('reference_point_m',[0.,0.,0.])
    return value


def normalize_stage_case(value):
    value.setdefault('verification_reject_once_stages',[])
    for field in ('monitor_sets','pressure_on_original_surface_only','initially_inactive_bodies'):
        if not value.get(field):value.pop(field,None)
    if value.get('boundary') is None:value.pop('boundary',None)
    else:normalize_conditions(value['boundary'])
    if value.get('joint_system') is None:value.pop('joint_system',None)
    else:
        joint=value['joint_system']
        if not joint.get('structure_bodies'):joint.pop('structure_bodies',None)
        for d in joint['joints']:
            if d.get('water') is None:d.pop('water',None)
            elif d['water']['kind']=='prescribed':d['water'].setdefault('level_m',None)
    anchor=value.get('anchor_system')
    if anchor is None:value.pop('anchor_system',None)
    else:
        for a in anchor['anchors']:
            for point in [a['head']]+a['bonded_path']:point.setdefault('body',None)
    beam=value.get('beam_system')
    if beam is None:value.pop('beam_system',None)
    else:
        for key in ('rotation_constraints','translation_constraints','nodal_loads'):beam.setdefault(key,[])
        for node in beam['nodes']:node.setdefault('soil_node_set',None);node.setdefault('match_soil_node',False)
        for member in beam['members']:
            member.setdefault('initial_axial_force_kn',0.)
            for key in ('constant_line_load_kn_m','boundary_line_load_kn_m'):member.setdefault(key,[0.]*3)
            member.setdefault('rotational_releases',[[False]*3 for _ in range(2)])
        for load in beam['nodal_loads']:
            for key in ('constant_force_kn','boundary_force_kn','constant_moment_kn_m','boundary_moment_kn_m'):load.setdefault(key,[0.]*3)
    if value['settings'].get('linear_solver',{}).get('method')=='dense_lu':value['settings'].pop('linear_solver')
    for stage in value['stages']:
        if stage['kind']!='reset_display':stage.setdefault('initial_increment',None)
        if stage['kind']=='steady_seepage':
            flow=stage['flow'];flow.setdefault('fluxes',[])
            flow.setdefault('settings',{'relative_linear_tolerance':1e-11,'pressure_head_tolerance_m':1e-9,'maximum_free_heads':300000,
                'linear_solver':{'method':'sparse_gmres_ilu0','options':{'restart':60,'maximum_iterations':3000,'pivot_tolerance':1e-14}}})
            for material in flow['materials']:material.setdefault('volume_source_per_s',0.)
            for head in flow['heads']:head.setdefault('gradient',[0.,0.,0.]);head.setdefault('reference_point_m',[0.,0.,0.])
            for flux in flow['fluxes']:
                if flux['kind']=='vector':flux.setdefault('gradient_per_s',[[0.]*3 for _ in range(3)]);flux.setdefault('reference_point_m',[0.,0.,0.])
    return value


def normalize_ground(value):
    for key in ('layers','bodies','excavation_bodies','construction_material_ids'):
        if not value.get(key):value.pop(key,None)
    cfg=value.get('initialization',{}).get('geostatic')
    if cfg is not None:
        base=cfg.get('base_loads')
        if base is not None:
            base.setdefault('nodal_loads',[]);base.setdefault('surface_loads',[])
            normalized=normalize_conditions({'constraints':[],**base});normalized.pop('constraints')
            cfg['base_loads']=normalized
            if not normalized['nodal_loads'] and not normalized['surface_loads']:cfg.pop('base_loads')
        for key in ('hydrostatic_face_sets','verification_reject_once'):
            if not cfg.get(key):cfg.pop(key,None)
    return value


def cells(mesh):
    return np.asarray(mesh.get('tetrahedra10',mesh['tetrahedra']),dtype=int)


def jacobi(n,alpha):
    j=np.zeros((n,n))
    for k in range(n):
        t=2*k+alpha;j[k,k]=0 if t==0 else -alpha**2/(t*(t+2))
        if k:
            j[k,k-1]=j[k-1,k]=2/t*np.sqrt(k*k*(k+alpha)**2/((t-1)*(t+1)))
    x,v=np.linalg.eigh(j)
    return (x+1)/2,v[0]**2/(alpha+1)


def volume_rule(mesh,independent=False,independent_order=5):
    if not mesh.get('tetrahedra10'):return np.array([[.25]*3]),np.array([1/6])
    if independent:
        x,w=np.polynomial.legendre.leggauss(independent_order);x=(x+1)/2;w=w/2
        axis=[(x,w*(1-x)**2),(x,w*(1-x)),(x,w)]
    else:axis=[jacobi(mesh.get('quadrature_order',3),a) for a in (2,1,0)]
    r=[];weights=[]
    for u,wu in zip(*axis[0]):
        for v,wv in zip(*axis[1]):
            for z,wz in zip(*axis[2]):
                r.append([u,(1-u)*v,(1-u)*(1-v)*z]);weights.append(wu*wv*wz)
    return np.array(r),np.array(weights)


def basis(natural,quadratic=True):
    r=np.asarray(natural);l=np.column_stack((1-r.sum(axis=1),r));dl=np.vstack((-np.ones(3),np.eye(3)))
    if not quadratic:return l,np.broadcast_to(dl,(len(r),4,3)).copy()
    shape=np.empty((len(r),10));gradient=np.empty((len(r),10,3))
    shape[:,:4]=l*(2*l-1);gradient[:,:4]=(4*l-1)[:,:,None]*dl
    for k,(i,j) in enumerate(EDGES):
        shape[:,k+4]=4*l[:,i]*l[:,j]
        gradient[:,k+4]=4*(l[:,i,None]*dl[j]+l[:,j,None]*dl[i])
    return shape,gradient


def volume(mesh,ei,independent=False,independent_order=5):
    xyz=np.asarray(mesh['nodes']);cell=cells(mesh)[ei];p=xyz[cell];r,w=volume_rule(mesh,independent,independent_order)
    shape,gradient=basis(r,len(cell)==10);j=np.einsum('ai,qaj->qij',p-p[0],gradient)
    determinants=np.linalg.det(j);assert np.all(determinants>0)
    g=np.einsum('qai,qij->qaj',gradient,np.linalg.inv(j));b=np.zeros((len(r),6,3*len(cell)))
    for i in range(len(cell)):
        gx,gy,gz=g[:,i].T;k=3*i
        b[:,0,k]=gx;b[:,1,k+1]=gy;b[:,2,k+2]=gz
        b[:,3,k]=gy;b[:,3,k+1]=gx;b[:,4,k+1]=gz;b[:,4,k+2]=gy;b[:,5,k]=gz;b[:,5,k+2]=gx
    return {'natural':r,'position':p[0]+shape@(p-p[0]),'weight':w*determinants,'shape':shape,'b':b,'nodes':cell}


def topology(mesh):
    result={}
    for ei,tet in enumerate(mesh['tetrahedra']):
        for f in FACES:
            key=tuple(sorted(tet[i] for i in f));result.setdefault(key,[]).append(ei)
    return result


def face_nodes(mesh,owner,face):
    cell=cells(mesh)[owner]
    for f in FACES:
        if set(cell[i] for i in f)!=set(face):continue
        nodes=[cell[i] for i in f]
        if len(cell)==10:
            for a,b in ((0,1),(1,2),(2,0)):
                edge=next(k for k,e in enumerate(EDGES) if set(e)=={f[a],f[b]});nodes.append(cell[edge+4])
        return np.asarray(nodes)
    raise ValueError('owner/face mismatch')


def face(mesh,owner,key):
    ids=face_nodes(mesh,owner,key);p=np.asarray(mesh['nodes'])[ids]
    x,w=np.polynomial.legendre.leggauss(5);x=(x+1)/2;w=w/2
    r=np.array([[u,(1-u)*v] for u in x for v in x])
    weights=np.array([wu*wv*(1-u) for u,wu in zip(x,w) for wv in w])
    l=np.column_stack((1-r.sum(axis=1),r));dl=np.array([[-1,-1],[1,0],[0,1]])
    if len(ids)==3:shape=l;gradient=np.broadcast_to(dl,(len(r),3,2))
    else:
        shape=np.empty((len(r),6));gradient=np.empty((len(r),6,2))
        shape[:,:3]=l*(2*l-1);gradient[:,:3]=(4*l-1)[:,:,None]*dl
        for k,(i,j) in enumerate(((0,1),(1,2),(2,0))):
            shape[:,k+3]=4*l[:,i]*l[:,j];gradient[:,k+3]=4*(l[:,i,None]*dl[j]+l[:,j,None]*dl[i])
    j=np.einsum('ai,qaj->qij',p-p[0],gradient)
    normal=np.cross(j[:,:,0],j[:,:,1])*weights[:,None]
    return {'nodes':ids,'shape':shape,'position':p[0]+shape@(p-p[0]),'area_vector':normal}


def elastic_tensor(e,nu):
    g=e/(2*(1+nu));lam=e*nu/((1+nu)*(1-2*nu));d=np.zeros((6,6));d[:3,:3]=lam
    for i in range(3):d[i,i]+=2*g;d[i+3,i+3]=g
    return d


def conditions_load(mesh,conditions,active,skip_removed_surfaces=False):
    xyz=np.asarray(mesh['nodes']);f=np.zeros_like(xyz);owners=topology(mesh)
    used=np.unique(cells(mesh)[active])
    for load in conditions.get('nodal_loads',[]):
        ids=np.array(mesh['node_sets'][load['node_set']]);assert np.isin(ids,used).all()
        f[ids]+=np.asarray(load['total_force_kn'])/len(ids)
    for load in conditions.get('surface_loads',[]):
        for key in mesh['face_sets'][load['face_set']]:
            live=[ei for ei in owners[tuple(sorted(key))] if active[ei]]
            if skip_removed_surfaces and not live:continue
            assert len(live)==1
            q=face(mesh,live[0],key);dx=q['position']-load.get('reference_point_m',[0,0,0]);area=q['area_vector']
            if load['kind']=='pressure':
                p=load['pressure_kpa']+dx@load.get('gradient_kpa_per_m',[0,0,0]);force=-p[:,None]*area
            else:
                t=load['traction_kpa']+dx@np.asarray(load.get('gradient_kpa_per_m',np.zeros((3,3)))).T
                force=t*np.linalg.norm(area,axis=1)[:,None]
            f[q['nodes']]+=q['shape'].T@force
    return f


def staged_constraint_space(mesh,case,state):
    xyz=np.array(mesh['nodes']);active=np.array(state['active_elements']);used=np.unique(cells(mesh)[active])
    kinematics=case['settings']['kinematics'];factor=state.get('boundary_factor',0);eq=[[] for _ in xyz];rhs=[[] for _ in xyz]
    if kinematics=='side_rollers':
        for n in used:
            for a in range(3):
                fixed=abs(xyz[n,2]-xyz[:,2].min())<1e-8 if a==2 else min(abs(xyz[n,a]-xyz[:,a].min()),abs(xyz[n,a]-xyz[:,a].max()))<1e-8
                if fixed:eq[n].append(np.eye(3)[a]);rhs[n].append(0.)
    else:assert kinematics=='explicit'
    for c in case['boundary']['constraints']:
        d=np.asarray(c['direction']);norm=np.linalg.norm(d)
        for n in mesh['node_sets'][c['node_set']]:
            if n not in used:continue
            eq[n].append(d/norm);rhs[n].append(factor*(c['value_m']+np.asarray(c.get('gradient',[0,0,0]))@(xyz[n]-c.get('reference_point_m',[0,0,0])))/norm)
    basis={};constraint_error=0.;u=np.asarray(state['displacement_m'])
    for n in used:
        if eq[n]:
            a=np.array(eq[n]);b=np.array(rhs[n]);_,s,vh=np.linalg.svd(a,full_matrices=True)
            rank=int((s>1e-10*max(s.max(),1)).sum());basis[n]=vh[rank:].T
            constraint_error=max(constraint_error,float(np.max(np.abs(a@u[n]-b))))
        else:basis[n]=np.eye(3)
    return basis,constraint_error


def initial_fields(mesh,ground,ei,position):
    regions={r['id']:r for r in mesh['regions']};layers={l['body']:l for l in ground['layers']};materials={m['id']:m for m in ground['materials']}
    region=regions[mesh['region_ids'][ei]];layer=layers[region['body']];material=materials[region['material_id']]
    z=np.asarray(position)[:,2];total=np.full(len(z),ground['initialization']['surface_surcharge_kpa'])
    for r in mesh['regions']:
        total+=np.maximum(r['bounds_m'][1][2]-np.maximum(z,r['bounds_m'][0][2]),0)*materials[r['material_id']]['unit_weight_kn_m3']
    h=ground['initialization']['water_level_m'];pore=np.zeros(len(z)) if h is None else np.maximum(h-z,0)*ground['initialization']['water_unit_weight_kn_m3']
    s=np.zeros((len(z),6));s[:,:3]=(total-pore)[:,None]*np.array([layer['k0_x'],layer['k0_y'],1])
    return s,pore,material,layer


def verify_gmsh_basis():
    import gmsh
    gmsh.initialize()
    try:
        local=np.asarray(gmsh.model.mesh.getElementProperties(11)[4]).reshape(-1,3)
        permutation=list(range(4))+[int(np.argmin(np.linalg.norm(local-(local[i]+local[j])/2,axis=1))) for i,j in EDGES]
        points=np.array([[.13,.21,.17],[.4,.2,.1],[.08,.04,.4]])
        shape,grad=basis(points)
        _,values,_=gmsh.model.mesh.getBasisFunctions(11,points.reshape(-1).tolist(),'Lagrange')
        _,derivatives,_=gmsh.model.mesh.getBasisFunctions(11,points.reshape(-1).tolist(),'GradLagrange')
        a=np.asarray(values).reshape(-1,10)[:,permutation];b=np.asarray(derivatives).reshape(-1,10,3)[:,permutation]
        assert np.max(np.abs(a-shape))<1e-13 and np.max(np.abs(b-grad))<1e-13
        return {'gmsh_shape_maximum_error':float(np.max(np.abs(a-shape))),'gmsh_gradient_maximum_error':float(np.max(np.abs(b-grad)))}
    finally:gmsh.finalize()
