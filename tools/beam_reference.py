"""Independent 3D frame: section flexibility by virtual work, explicit releases.

Unlike Rust's closed nodal matrix and Schur complement, this implementation
keeps all frame translations/rotations and released end rotations in a dense
system. Soil bindings are imposed by a kinematic transformation.
"""
import numpy as np


def constraint_space(constraints,node):
    eq=[c for c in constraints if c['node_id']==node]
    if not eq:return np.zeros(3),np.eye(3),np.zeros((3,0))
    a=np.asarray([c['direction'] for c in eq],float);b=np.asarray([c['value'] for c in eq]);_,s,vh=np.linalg.svd(a,full_matrices=True)
    rank=int((s>1e-10*max(s.max(),1.)).sum());offset=np.linalg.lstsq(a,b,rcond=None)[0]
    assert np.max(np.abs(a@offset-b))<1e-10
    return offset,vh[rank:].T,vh[:rank].T


def element(a,b,local_y,section):
    dx=b-a;l=np.linalg.norm(dx);ex=dx/l;ey=np.asarray(local_y,float);ey-=ex*(ey@ex);ey/=np.linalg.norm(ey);r=np.array([ex,ey,np.cross(ex,ey)])
    e=section['young_kpa'];g=e/(2*(1+section['poisson']));form=section['formulation']
    comp=[1/(e*section['area_m2']),1/(g*section['torsion_j_m4']),1/(e*section['iy_m4']),1/(e*section['iz_m4']),0.,0.]
    if form['kind']=='timoshenko':comp[4:]=[1/(g*form['shear_area_y_m2']),1/(g*form['shear_area_z_m2'])]
    f=np.zeros((6,6));xi,wi=np.polynomial.legendre.leggauss(5)
    for x,w in zip((xi+1)*l/2,wi*l/2):
        h=np.zeros((6,6));h[0,0]=h[1,1]=1.;h[2,2]=h[3,4]=-(1-x/l);h[2,3]=h[3,5]=x/l
        h[4,4:6]=-1/l;h[5,2:4]=1/l;f+=h.T@np.diag(comp)@h*w
    basic=np.zeros((6,12));basic[0,[0,6]]=[-1,1];basic[1,[3,9]]=[-1,1]
    basic[2,[2,8,4]]=[-1/l,1/l,1];basic[3,[2,8,10]]=[-1/l,1/l,1]
    basic[4,[1,7,5]]=[1/l,-1/l,1];basic[5,[1,7,11]]=[1/l,-1/l,1]
    k=basic.T@np.linalg.solve(f,basic);t=np.kron(np.eye(4),r)
    return l,r,t,k


def uniform_load(l,r,global_line):
    q=r@global_line;f=np.zeros(12);xi,wi=np.polynomial.legendre.leggauss(5)
    for z,w in zip((xi+1)/2,wi*l/2):
        f[[0,6]]+=q[0]*np.array([1-z,z])*w
        for ids,component,sign in (([1,5,7,11],q[1],1),([2,4,8,10],q[2],-1)):
            h=np.array([1-3*z*z+2*z**3,sign*l*(z-2*z*z+z**3),3*z*z-2*z**3,sign*l*(-z*z+z**3)])
            f[ids]+=component*h*w
    return f


class System:
    def __init__(self,mesh,cfg):
        self.mesh=mesh;self.cfg=cfg;self.nodes=cfg['nodes'];self.members=cfg['members'];self.node_ids={n['id']:i for i,n in enumerate(self.nodes)}
        self.attach=[];self.positions=[]
        for n in self.nodes:
            if n.get('match_soil_node'):
                found=np.flatnonzero(np.linalg.norm(np.asarray(mesh['nodes'])-n['position_m'],axis=1)<=1e-8);assert len(found)==1;node=int(found[0])
            elif n.get('soil_node_set') is not None:
                selected=mesh['node_sets'][n['soil_node_set']];assert len(selected)==1;node=selected[0]
            else:node=None
            self.attach.append(node);self.positions.append(n['position_m'] if node is None else mesh['nodes'][node])
        self.positions=np.asarray(self.positions);self.data=[];self.released={};self.size=6*len(self.nodes);sections={s['id']:s for s in cfg['sections']}
        for mi,m in enumerate(self.members):
            ends=[self.node_ids[n] for n in m['nodes']];s=sections[m['section_id']];l,r,t,k=element(*self.positions[ends],m['local_y'],s)
            for end in range(2):
                for a in range(3):
                    if m.get('rotational_releases',[[False]*3]*2)[end][a]:self.released[mi,end,a]=self.size;self.size+=1
            self.data.append((ends,s,l,r,t,k))

    def initial_state(self):
        return {'active_members':[m['initially_active'] for m in self.members],'reference_global_kinematics':[[0.]*12 for _ in self.members],
            'node_displacement_m':[[0.]*3 for _ in self.nodes],'node_rotation_rad':[[0.]*3 for _ in self.nodes]}

    def active_nodes(self,state):
        active=np.zeros(len(self.nodes),bool)
        for mi,(ends,*_) in enumerate(self.data):
            if state['active_members'][mi]:active[ends]=True
        return active

    def change(self,state,soil_u,stage):
        import copy
        result=copy.deepcopy(state);install=stage['kind']=='install_beams';old=self.active_nodes(state)
        for name in stage['members']:
            mi=next(i for i,m in enumerate(self.members) if m['id']==name);assert state['active_members'][mi]!=install;result['active_members'][mi]=install
            if install:
                ends=self.data[mi][0]
                for end,n in enumerate(ends):
                    if self.attach[n] is not None:result['node_displacement_m'][n]=list(soil_u[self.attach[n]])
                    elif stage['placement']=='cad_new_nodes' and not old[n]:result['node_displacement_m'][n]=[0.]*3
                    if stage['placement']=='cad_new_nodes' and not old[n]:result['node_rotation_rad'][n]=[0.]*3
                    result['reference_global_kinematics'][mi][6*end:6*end+6]=result['node_displacement_m'][n]+result['node_rotation_rad'][n]
        return result

    def assemble(self,state,factor):
        k=np.zeros((self.size,self.size));rest=np.zeros(self.size);load=np.zeros(self.size);effective=np.zeros((len(self.nodes),3));info={};external_f=np.zeros(3);external_m=np.zeros(3)
        for mi,(ends,s,l,r,t,local_k) in enumerate(self.data):
            if not state['active_members'][mi]:continue
            member=self.members[mi];p=np.zeros((12,self.size))
            for end,n in enumerate(ends):
                p[6*end:6*end+3,6*n:6*n+3]=r
                for a in range(3):
                    if (mi,end,a) in self.released:p[6*end+3+a,self.released[mi,end,a]]=1.
                    else:p[6*end+3+a,6*n+3:6*n+6]=r[a]
            line=np.asarray(member.get('constant_line_load_kn_m',[0.]*3))+factor*np.asarray(member.get('boundary_line_load_kn_m',[0.]*3));line=line.astype(float);line[2]-=s['unit_weight_kn_m3']*s['area_m2']
            f=uniform_load(l,r,line);reference=t@state['reference_global_kinematics'][mi];pre=np.zeros(12);pre[0]=-member.get('initial_axial_force_kn',0.);pre[6]=-pre[0]
            k+=p.T@local_k@p;rest+=p.T@(pre-local_k@reference);load+=p.T@f
            # Independently recover the load convention used by the condensed
            # soil interface, while preserving all release DOFs in our solve.
            released=[6*end+3+a for end in range(2) for a in range(3) if (mi,end,a) in self.released];retained=[i for i in range(12) if i not in released];eff=f.copy()
            if released:
                eff[:]=0;eff[retained]=f[retained]-local_k[np.ix_(retained,released)]@np.linalg.solve(local_k[np.ix_(released,released)],f[released])
            global_f=t.T@f;effective_f=t.T@eff
            for end,n in enumerate(ends):
                effective[n]+=effective_f[6*end:6*end+3];external_f+=global_f[6*end:6*end+3]
                external_m+=np.cross(self.positions[n],global_f[6*end:6*end+3])+global_f[6*end+3:6*end+6]
            info[mi]=(p,local_k,reference,pre,f,l,r,line)
        active=self.active_nodes(state)
        for d in self.cfg.get('nodal_loads',[]):
            n=self.node_ids[d['node_id']];f=np.asarray(d.get('constant_force_kn',[0.]*3))+factor*np.asarray(d.get('boundary_force_kn',[0.]*3));m=np.asarray(d.get('constant_moment_kn_m',[0.]*3))+factor*np.asarray(d.get('boundary_moment_kn_m',[0.]*3));assert active[n] or not np.r_[f,m].any()
            load[6*n:6*n+3]+=f;load[6*n+3:6*n+6]+=m;effective[n]+=f;external_f+=f;external_m+=np.cross(self.positions[n],f)+m
        return k,rest,load,effective,external_f,external_m,info

    def space(self,state,soil_u,soil_basis=None):
        active=self.active_nodes(state);columns=[];offset=np.zeros(self.size);normals=[];binding=np.zeros((self.size,0 if soil_basis is None else soil_basis.shape[1]))
        for i,n in enumerate(self.nodes):
            if not active[i]:offset[6*i:6*i+3]=state['node_displacement_m'][i];offset[6*i+3:6*i+6]=state['node_rotation_rad'][i];continue
            for part,cs in ((0,self.cfg.get('translation_constraints',[])),(1,self.cfg.get('rotation_constraints',[]))):
                if part==0 and self.attach[i] is not None:
                    node=self.attach[i]
                    if soil_basis is None:offset[6*i:6*i+3]=soil_u[node]
                    else:binding[6*i:6*i+3]=soil_basis[3*node:3*node+3]
                    continue
                u0,b,fixed=constraint_space(cs,n['id']);offset[6*i+3*part:6*i+3*part+3]=u0;normals.append((6*i+3*part,fixed))
                for col in b.T:
                    v=np.zeros(self.size);v[6*i+3*part:6*i+3*part+3]=col;columns.append(v)
        for (mi,end,a),index in self.released.items():
            if state['active_members'][mi]:v=np.zeros(self.size);v[index]=1.;columns.append(v)
        q=np.column_stack(columns) if columns else np.zeros((self.size,0))
        return offset,binding,q,normals

    def solve_interface(self,state,soil_u,factor):
        k,rest,load,eff,extf,extm,info=self.assemble(state,factor);offset,_,q,normals=self.space(state,soil_u)
        d=offset+q@np.linalg.solve(q.T@k@q,q.T@(load-rest-k@offset)) if q.shape[1] else offset
        return self.results(state,d,k,rest,load,eff,extf,extm,info,normals)

    def results(self,state,d,k,rest,load,eff,extf,extm,info,normals):
        import copy
        result=copy.deepcopy(state);active=self.active_nodes(state);residual=k@d+rest-load;support_f=np.zeros(3);support_m=np.zeros(3);interface={};attached_load={};members={}
        for i in np.flatnonzero(active):
            result['node_displacement_m'][i]=d[6*i:6*i+3].tolist();result['node_rotation_rad'][i]=d[6*i+3:6*i+6].tolist()
            if self.attach[i] is not None:
                node=self.attach[i];interface[node]=interface.get(node,np.zeros(3))+residual[6*i:6*i+3];attached_load[node]=attached_load.get(node,np.zeros(3))+eff[i]
        for start,fixed in normals:
            r=fixed@(fixed.T@residual[start:start+3])
            if start%6==0:support_f+=r;support_m+=np.cross(self.positions[start//6],r)
            else:support_m+=r
        for mi,(p,local_k,reference,pre,f,l,r,line) in info.items():
            relative=p@d-reference;action=local_k@relative+pre-f
            samples=[];distributed=r@line
            for x in (0.,l/2,l):
                section_force=-action[:3]-distributed*x
                section_moment=-action[3:6]+np.cross([x,0.,0.],action[:3]+distributed*x/2)
                samples.append(np.r_[section_force,section_moment])
            members[self.members[mi]['id']]={'local_end_actions':action,'local_end_relative_rotation_rad':np.array([relative[3:6],relative[9:12]]),
                                           'section_forces_at_start_mid_end':np.asarray(samples),'local_axes':r,'length_m':l}
        return {'state':result,'interface_residual':interface,'attached_load':attached_load,'external_force':extf,'external_moment':extm,'support_force':support_f,'support_moment':support_m,'members':members,'residual':residual,'positions':self.positions,'node_active':active}


def verify_state(system,native,report,soil_u,factor):
    reference=system.solve_interface(native,soil_u,factor)
    assert np.allclose(reference['state']['node_displacement_m'],native['node_displacement_m'],atol=1e-10,rtol=1e-8)
    assert np.allclose(reference['state']['node_rotation_rad'],native['node_rotation_rad'],atol=1e-10,rtol=1e-8)
    for name,key in (('external_force_kn','external_force'),('external_moment_origin_kn_m','external_moment'),('support_force_kn','support_force'),('support_moment_origin_kn_m','support_moment')):assert np.allclose(report[name],reference[key],atol=1e-6,rtol=1e-8),(name,report[name],reference[key])
    assert {m['id'] for m in report['members']}==set(reference['members'])
    assert np.array_equal(report['node_active'],reference['node_active'])
    assert np.allclose(report['node_positions_m'],reference['positions'],atol=1e-12,rtol=1e-12)
    for member in report['members']:
        expected=reference['members'][member['id']]
        for key in ('local_end_actions','local_end_relative_rotation_rad','section_forces_at_start_mid_end','local_axes','length_m'):assert np.allclose(member[key],expected[key],atol=1e-7,rtol=1e-8),(member['id'],key)
    return reference


def coupled_forces(system,state,soil_internal,soil_external,report=None):
    """Match the soil-interface force convention, retaining explicit frame DOFs."""
    u=np.asarray(state['displacement_m']);factor=state.get('boundary_factor',0.)
    result=(system.solve_interface(state['beams'],u,factor) if report is None else
            verify_state(system,state['beams'],report,u,factor))
    internal=soil_internal.copy();external=soil_external.copy()
    for node,load in result['attached_load'].items():
        external[node]+=load;internal[node]+=load+result['interface_residual'][node]
    return internal,external,result
