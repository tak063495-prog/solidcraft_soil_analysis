"""Reject nonfinite external GP data before any error reduction (Issue #14)."""
import copy,json,subprocess,sys,tempfile,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
import fem_reference as fem
from joint_water_reference import Water
from verify_soil_source_column import audit

def fixture():
    read=lambda n:json.loads((ROOT/'examples'/n).read_text())
    mesh=read('joint-water-column-tet10-mesh.json')
    ground=fem.normalize_ground(read('ground-wet-plate-beam.json'))
    case=fem.normalize_stage_case(read('stages-joint-water-soil.json'))
    stage=next(s for s in case['stages'] if s['kind']=='steady_seepage')
    stage['flow']['materials'][0]['volume_source_per_s']=2*stage['flow']['materials'][0]['hydraulic_conductivity_m_s'][2][2]
    water=Water(mesh,ground,case);xyz=np.asarray(mesh['nodes']);used=np.unique(fem.cells(mesh)[water.soil])
    points=[]
    for i in np.flatnonzero(water.soil):
        q=fem.volume(mesh,i)
        for pos,w in zip(q['position'],q['weight'],strict=True):
            z=pos[2];points.append(dict(element_index=int(i),position_m=pos.tolist(),weight_m3=float(w),head_m=2+z-z*z,pressure_head_m=2-z*z,darcy_flux_m_s=[0.,0.,2e-4*(z-.5)]))
    total=2e-4*sum(p['weight_m3'] for p in points)
    balance={k:total for k in ('source_m3_s','absolute_source_m3_s','dirichlet_reaction_outflow_m3_s','gradient_boundary_outflow_m3_s')}
    balance.update(neumann_outflow_m3_s=0.,relative_nodal_balance_error=0.)
    solution=dict(complete=True,active_elements=water.soil.tolist(),points=points,balances=[balance],boundary_fluxes=[dict(boundary_kind='head',gradient_outflow_m3_s=total),dict(boundary_kind='no_flow',gradient_outflow_m3_s=0.)])
    field=dict(defined_nodes=np.isin(np.arange(len(xyz)),used).tolist(),head_m=(2+xyz[:,2]-xyz[:,2]**2).tolist())
    index=case['stages'].index(stage);gamma=ground['initialization']['water_unit_weight_kn_m3']
    frames=[{} for _ in range(index+2)];frames[-1]=dict(state=dict(joint_water=dict(physical_pressure_kpa=[gamma*1.75]),pore_pressure_kpa=[gamma*(2-p['position_m'][2]**2) for p in points]))
    native=dict(complete=True,flow_audits=[dict(stage_index=index,solution=solution,joint_flow=dict(original_field=field))],frames=frames,integration_points=copy.deepcopy(points))
    return mesh,ground,case,native

class SoilGP(unittest.TestCase):
    def call(self,values,output):
        mesh,ground,case,native=values;folder=output.parent
        for name,value in zip(('mesh','ground','stages','stage-results'),values,strict=True):(folder/(name+'.json')).write_text(json.dumps(value))
        return subprocess.run([sys.executable,str(ROOT/'tools/verify_soil_source_column.py'),'--evidence',str(folder),'--output',str(output)],capture_output=True,text=True)

    def test_normal_api_and_cli(self):
        values=fixture();self.assertEqual(len(values[-1]['flow_audits'][0]['solution']['points']),162)
        self.assertTrue(audit(*values)['passed'])
        with tempfile.TemporaryDirectory() as d:
            output=Path(d)/'report.json';p=self.call(values,output)
            self.assertEqual(p.returncode,0,p.stderr);self.assertTrue(json.loads(output.read_text())['passed'])

    def test_nonfinite_mismatch_and_shape_api_and_cli(self):
        base=fixture()
        for key in ('head_m','pressure_head_m','darcy_flux_m_s','position_m'):
            for count in (1,162):
                for value in (float('nan'),float('inf'),-float('inf'),100.):
                    with self.subTest(key=key,count=count,value=value),tempfile.TemporaryDirectory() as d:
                        values=copy.deepcopy(base)
                        for point in values[-1]['flow_audits'][0]['solution']['points'][:count]:point[key]=[value]*3 if key in ('position_m','darcy_flux_m_s') else value
                        with self.assertRaises((ValueError,AssertionError)):audit(*values)
                        output=Path(d)/'report.json';p=self.call(values,output)
                        self.assertNotEqual(p.returncode,0,p.stdout);self.assertFalse(output.exists())
            for shape in ([],[0.],[[0.,0.,0.]]):
                with self.subTest(key=key,shape=shape),tempfile.TemporaryDirectory() as d:
                    values=copy.deepcopy(base);values[-1]['flow_audits'][0]['solution']['points'][0][key]=shape
                    with self.assertRaises(ValueError):audit(*values)
                    output=Path(d)/'report.json';self.assertNotEqual(self.call(values,output).returncode,0);self.assertFalse(output.exists())

if __name__=='__main__':unittest.main()
