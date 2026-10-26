"""Check both uncondensed MCC/member rate forms against a full elastic solve."""
import os
for name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[name]='1'
from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from test_mcc_factor_reference import specimen
from beam_reference import System as Members
from mcc_structural_reference import StructuralSystem
import mcc_rate_fem as rate


def configuration(mesh):
    nodes=[dict(id=f'n{i}',position_m=mesh['nodes'][i],soil_node_set=f'n{i}') for i in range(3)]
    nodes.append(dict(id='free',position_m=(np.asarray(mesh['nodes'][1])+[.3,.4,.8]).tolist(),soil_node_set=None))
    return dict(schema='soil-fem-structures/1',nodes=nodes,
        sections=[dict(id='b',young_kpa=2e6,poisson=.2,area_m2=.03,iy_m4=1e-4,iz_m4=2e-4,
            torsion_j_m4=1e-4,unit_weight_kn_m3=24.,formulation=dict(kind='euler_bernoulli'))],
        members=[dict(id='beam',nodes=['n1','n2'],section_id='b',local_y=[0.,0.,1.],initially_active=True,
            rotational_releases=[[False,False,True],[False,False,False]],boundary_line_load_kn_m=[.1,-.2,-.3])],
        plate_sections=[dict(id='p',young_kpa=3e6,poisson=.2,thickness_m=.08,unit_weight_kn_m3=24.,drilling_factor=1e-4)],
        plates=[dict(id='plate',nodes=['n0','n1','free'],section_id='p',initially_active=True,boundary_traction_kpa=[.2,-.1,-.4])],
        rotation_constraints=[dict(node_id='n0',direction=np.eye(3)[a].tolist(),value=0.) for a in range(3)],
        translation_constraints=[],nodal_loads=[],rotation_length_scale_m=1.,maximum_condensed_dofs=100)


class StructuralRate(unittest.TestCase):
    def check_form(self,projected):
        mesh,ground,case=specimen();system=StructuralSystem(mesh,ground,[True],case,projected=projected)
        members=Members(mesh,configuration(mesh));history=members.initial_state();system.configure_structures(members,history,.7)
        state=np.zeros((system.weight.size,14));state[:,:3]=100.;state[:,6:8]=np.log([200.,2.])
        elastic=rate.terms(state,system.parameters)[1]
        force=np.random.default_rng(60261026).normal(size=system.xyz.size)
        du,de=system.linear_rate(elastic,force)
        identity=np.eye(system.xyz.size)
        transform=np.stack([system.lift(identity[i]) for i in system.free],axis=1)
        ksoil=np.zeros((system.xyz.size,system.xyz.size))
        for e,dofs in enumerate(system.dofs):
            local=np.einsum('qai,qab,qbj,q->ij',system.b[e],elastic.reshape(*system.weight.shape,6,6)[e],system.b[e],system.weight[e])
            ksoil[np.ix_(dofs,dofs)]+=local
        kb,_,f0,*_=members.assemble(history,0.);_,_,f1,*_=members.assemble(history,1.)
        _,binding,q,_=members.space(history,np.zeros_like(system.xyz),transform)
        matrix=np.block([[transform.T@ksoil@transform+binding.T@kb@binding,binding.T@kb@q],
                         [q.T@kb@binding,q.T@kb@q]])
        load=(f1-f0)*.7
        rhs=np.r_[transform.T@force+binding.T@load,q.T@load]
        direct=np.linalg.solve(matrix,rhs);expected=transform@direct[:len(system.free)]
        np.testing.assert_allclose(du,expected,rtol=1e-9,atol=1e-11)
        strain=-np.einsum('eqai,ei->eqa',system.b,expected[system.dofs]).reshape(-1,6)
        np.testing.assert_allclose(de,strain,rtol=1e-9,atol=1e-11)
        self.assertLess(max(system.maximum_block_residuals.values()),1e-10)
        self.assertGreater(system.retained_structural_dofs,0)

    def test_projected_rates_match_full_soil_plate_frame_released_rotation_system(self):self.check_form(True)
    def test_gp_volume_rates_match_full_soil_plate_frame_released_rotation_system(self):self.check_form(False)


if __name__=='__main__':unittest.main()
