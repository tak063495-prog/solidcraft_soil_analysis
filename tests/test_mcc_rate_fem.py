"""Physics checks for the independent continuum-rate verification reference."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
import numpy as np
from mcc_rate_fem import System, seed, solve, terms
from mcc_rate_midpoint import advance_cached


class RatePhysics(unittest.TestCase):
    def setUp(self):
        self.parameters=np.array([[1.2**2,.2,.04,.3]])
        self.state=np.zeros((1,14));self.state[0,:3]=100.;self.state[0,6:8]=np.log([100.,2.])

    def test_isotropic_nc_bulk_follows_lambda_and_elastic_bulk_follows_kappa(self):
        _,d,n,dn,den=terms(self.state,self.parameters)
        ep=d[0]-np.outer(dn[0],dn[0])/den[0];unit=np.r_[np.ones(3)/3,np.zeros(3)]
        self.assertAlmostEqual(float((d[0]@unit)[:3].mean()),2*100/.04,places=9)
        self.assertAlmostEqual(float((ep@unit)[:3].mean()),2*100/.2,places=9)
        np.testing.assert_allclose(ep,ep.T,atol=1e-10,rtol=1e-13)

    def test_shear_gradient_uses_engineering_strain_conjugacy(self):
        self.state[0,:6]=[120.,90.,90.,12.,-7.,4.]
        _,_,normal,_,_=terms(self.state,self.parameters)
        numerical=[];pc=np.exp(self.state[0,6]);m2=self.parameters[0,0]
        for axis in range(6):
            plus=self.state.copy();minus=self.state.copy();plus[0,axis]+=1e-3;minus[0,axis]-=1e-3
            numerical.append((terms(plus,self.parameters)[0][0]-terms(minus,self.parameters)[0][0])*m2*pc**2/(2e-3))
        np.testing.assert_allclose(normal[0],numerical,atol=1e-8,rtol=1e-9)

    def test_ocr_two_isotropic_elastic_centre_has_zero_yield_gradient(self):
        self.state[0,6]=np.log(200.)
        yf,d,n,dn,den=terms(self.state,self.parameters)
        self.assertAlmostEqual(yf[0],-.25)
        np.testing.assert_allclose(n,0.,atol=1e-12)
        np.testing.assert_allclose(dn,0.,atol=1e-8)
        self.assertGreater(np.linalg.eigvalsh(d[0]).min(),0.)

    def test_consistency_rate_cancels_yield_derivative_for_general_shear(self):
        # A non-coaxial 3D strain rate, exactly on an MCC ellipse.
        self.state[0,:6]=[130.,85.,85.,8.,-5.,3.]
        stress=self.state[0,:6];p=stress[:3].mean();dev=stress.copy();dev[:3]-=p
        q2=1.5*(dev[:3]@dev[:3]+2*(dev[3:]@dev[3:]))
        pc=p+q2/(self.parameters[0,0]*p);self.state[0,6]=np.log(pc)
        yf,d,n,dn,den=terms(self.state,self.parameters);self.assertLess(abs(yf[0]),1e-14)
        de=np.array([.002,.001,.003,.004,-.002,.001]);mult=float(dn[0]@de/den[0]);self.assertGreater(mult,0.)
        flow=mult*n[0];ds=d[0]@(de-flow);dpc=pc*2/(.2-.04)*flow[:3].sum()
        self.assertAlmostEqual(float(n[0]@ds-self.parameters[0,0]*p*dpc),0.,places=10)


class GlobalReference(unittest.TestCase):
    def specimen(self,rotate=False,quadratic=False):
        points=np.vstack((np.zeros(3),np.eye(3)));q=np.eye(3)
        if quadratic:points=np.vstack((points,[(points[a]+points[b])/2 for a,b in ((0,1),(1,2),(2,0),(0,3),(1,3),(2,3))]))
        if rotate:
            a=.43;b=.29
            q=np.array([[np.cos(a),-np.sin(a),0],[np.sin(a),np.cos(a),0],[0,0,1.]])@np.array([[1.,0,0],[0,np.cos(b),-np.sin(b)],[0,np.sin(b),np.cos(b)]])
        mesh=dict(nodes=(points@q.T).tolist(),tetrahedra=[[0,1,2,3]],material_ids=[1],node_sets={f'n{i}':[i] for i in range(4)})
        if quadratic:mesh.update(tetrahedra10=[list(range(10))],quadrature_order=3)
        pars=dict(critical_state_m=1.2,**{'lambda':.2,'kappa':.04},poisson=.3)
        ground=dict(materials=[dict(id=1,mcc=pars)])
        constraints=[]
        for node,axes in [(0,[0,1,2]),(1,[1,2]),(2,[2])]:
            for axis in axes:constraints.append(dict(node_set=f'n{node}',direction=q[:,axis].tolist(),value_m=0.))
        case=dict(settings=dict(kinematics='explicit'),boundary=dict(constraints=constraints))
        return mesh,ground,case,q

    def test_elastic_tetrahedron_with_rotated_directional_supports_has_analytic_affine_solution(self):
        for rotated,quadratic in ((False,False),(True,False),(False,True),(True,True)):
            mesh,ground,case,q=self.specimen(rotated,quadratic);system=System(mesh,ground,[True],case)
            state=np.zeros((system.b.shape[1],14));state[:,:3]=100.;state[:,6:8]=np.log([200.,2.])
            ds=np.array([[1.,.4,.2],[.4,2.,.3],[.2,.3,3.]])
            bulk=5000.;shear=3*bulk*(1-.6)/(2*1.3);young=2*shear*1.3
            de=(1.3*ds-.3*np.trace(ds)*np.eye(3))/young
            affine=-np.array([[de[0,0],2*de[0,1],2*de[0,2]],[0,de[1,1],2*de[1,2]],[0,0,de[2,2]]])
            expected=((np.asarray(mesh['nodes'])@q)@affine.T)@q.T
            rotated_stress=q@ds@q.T;vector=rotated_stress[[0,1,2,0,1,0],[0,1,2,1,2,2]]
            candidate=state.copy();candidate[:,:6]+=vector
            force=system.force(candidate)-system.force(state)
            du,dy,_,_=system.rate(state,force,np.zeros(len(state),dtype=bool))
            np.testing.assert_allclose(du.reshape(-1,3),expected,rtol=1e-10,atol=1e-13)
            np.testing.assert_allclose(dy[:,:6],np.broadcast_to(vector,dy[:,:6].shape),rtol=1e-10,atol=1e-12)
            self.assertLess(np.linalg.norm(system.reduce(force-system.force(candidate)+system.force(state))),1e-10)

    def test_sparse_linear_solution_and_indefinite_fallback_match_dense_independent_lu(self):
        for matrix in (np.array([[3.,1.],[1.,2.]]),np.array([[1.,2.],[2.,1.]])):
            system=System.__new__(System);system.free=np.arange(2);system.sparse_rows=np.array([0,0,1,1]);system.sparse_cols=np.array([0,1,0,1])
            system.inverse=np.arange(4);system.mask=np.ones(4,dtype=bool);system.linear_iterations=0;system.dense_fallbacks=0
            rhs=np.array([1.,0.]);actual=system.linear_solve(matrix[None],rhs)
            np.testing.assert_allclose(actual,np.linalg.solve(matrix,rhs),rtol=1e-12,atol=1e-12)
            self.assertEqual(system.dense_fallbacks,0 if matrix[0,0]==3 else 1)

    def test_seed_point_positions_are_checked_against_independent_geometry(self):
        mesh,ground,_,_=self.specimen();mesh.update(regions=[dict(id=1,body='soil')],region_ids=[1])
        row=dict(effective_stress_compression_positive_kpa=[100.,100.,100.,0.,0.,0.],preconsolidation_pressure_kpa=200.,void_ratio=1.,plastic_strain_compression_positive=[0.]*6)
        point=dict(element_index=0,local_index=0,material_id=1,position_m=[.25]*3,state=row)
        ground['initialization']=dict(geostatic=dict(seed=dict(kind='integration_points',mesh_sha256='a'*64,points=[point])))
        identity=dict(schema='soil-fem-mesh-identity/1',mesh_sha256='a'*64,integration_points=[dict(element_index=0,local_index=0,position_m=[.25]*3,weight_m3=1/6)])
        self.assertEqual(seed(mesh,ground,identity).shape,(1,1,14))
        point['position_m']=[.3,.25,.25]
        with self.assertRaisesRegex(ValueError,'position'):seed(mesh,ground,identity)

    def test_general_geostatic_equilibrium_and_boundary_cycle_follow_isotropic_closed_form(self):
        mesh,ground,case,_=self.specimen();mesh.update(regions=[dict(id=1,body='soil')],region_ids=[1],
            face_sets={'all':[[1,2,3],[0,3,2],[0,1,3],[0,2,1]]})
        ground['materials'][0]['unit_weight_kn_m3']=0.
        row=dict(effective_stress_compression_positive_kpa=[100.,100.,100.,0.,0.,0.],preconsolidation_pressure_kpa=200.,void_ratio=1.,plastic_strain_compression_positive=[0.]*6)
        ground['initialization']=dict(method='equilibrate',surface_surcharge_kpa=0.,water_level_m=None,
            geostatic=dict(seed=dict(kind='body_states',states={'soil':row}),initial_boundary_factor=0.,
                maximum_initial_engineering_strain=.05,base_loads=dict(surface_loads=[dict(kind='pressure',face_set='all',pressure_kpa=100.)])))
        case['boundary']['surface_loads']=[dict(kind='pressure',face_set='all',pressure_kpa=10.)]
        case['stages']=[dict(kind='boundary',name='load',target_factor=1.),dict(kind='boundary',name='unload',target_factor=0.)]
        result=solve(mesh,ground,case,100)
        np.testing.assert_allclose(result['reference_displacement_m'],0.,atol=1e-14)
        loaded=np.array(result['frames'][1]['material_state']).reshape(-1,14)[0]
        expected_v=2.-.04*np.log(1.1);expected_strain=np.log(2./expected_v)/3
        np.testing.assert_allclose(loaded[:3],110.,rtol=1e-11,atol=1e-10)
        self.assertAlmostEqual(np.exp(loaded[7]),expected_v,places=9)
        np.testing.assert_allclose(result['frames'][1]['displacement_m'],-np.asarray(mesh['nodes'])*expected_strain,rtol=1e-6,atol=1e-10)
        np.testing.assert_allclose(result['frames'][-1]['displacement_m'],0.,atol=2e-10)

    def test_cached_event_start_preserves_every_state_and_displacement_bit(self):
        for quadratic in (False,True):
            for ratio in (1.,1.1):
                mesh,ground,case,_=self.specimen(True,quadratic)
                old=System(mesh,ground,[True],case);cached=System(mesh,ground,[True],case)
                state=np.zeros((old.b.shape[1],14));state[:,:3]=100.;state[:,6:8]=np.log([100*ratio,2.])
                trial=state.copy();trial[:,:3]+=20.
                load=old.force(trial)-old.force(state);zero=np.zeros(np.asarray(mesh['nodes']).size)
                a,u,da=old.advance_uncached(state,zero,load,25);b,v,db=cached.advance(state,zero,load,25)
                np.testing.assert_array_equal(a,b);np.testing.assert_array_equal(u,v)
                self.assertEqual(a.tobytes(),b.tobytes());self.assertEqual(u.tobytes(),v.tobytes())
                self.assertEqual(da['accepted_steps'],db['accepted_steps']);self.assertEqual(da['events'],db['events'])
                if ratio>1:self.assertLess(db['rate_assemblies'],da['rate_assemblies'])
                a,u,da=old.advance_uncached(a,u,-load,25);b,v,db=cached.advance(b,v,-load,25)
                np.testing.assert_array_equal(a,b);np.testing.assert_array_equal(u,v)
                self.assertEqual(a.tobytes(),b.tobytes());self.assertEqual(u.tobytes(),v.tobytes())


if __name__=='__main__':unittest.main()
