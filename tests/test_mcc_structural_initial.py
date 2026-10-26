"""General CAD seed acceptance requires independently balanced original input."""
import contextlib
import io
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from prepare_mcc_structures import inputs
from run_mcc_structural_reference import solve


class GeneralStructuralInitial(unittest.TestCase):
    def run_seed(self,folder,pressure):
        root=Path(__file__).resolve().parents[1];mesh,ground,case=inputs(root)
        states={}
        for material in ground['materials']:
            material['unit_weight_kn_m3']=0.;p=material['mcc'];body=next(x['body'] for x in ground['layers'] if x['material_id']==material['id'])
            states[body]=dict(effective_stress_compression_positive_kpa=[100.]*3+[0.]*3,
                preconsolidation_pressure_kpa=200.,plastic_strain_compression_positive=[0.]*6,
                void_ratio=p['specific_volume_ncl_at_reference']-p['lambda']*math.log(200.)+p['kappa']*math.log(2.)-1)
        ground['bodies']=[dict(body=x['body'],material_id=x['material_id']) for x in ground.pop('layers')]
        ground['initialization'].update(method='equilibrate',surface_surcharge_kpa=pressure,
            geostatic=dict(seed=dict(kind='body_states',states=states),initial_boundary_factor=0.,
                maximum_initial_engineering_strain=.05,relative_volume_tolerance=1e-8,base_loads=dict(surface_loads=[])))
        case['stages']=[];case['verification_reject_once_stages']=[]
        paths=[folder/n for n in ('mesh.json','ground.json','case.json')]
        for path,value in zip(paths,(mesh,ground,case),strict=True):path.write_text(json.dumps(value),encoding='utf-8')
        output=folder/'reference.json'
        with contextlib.redirect_stdout(io.StringIO()):solve(*paths,output,2)
        return json.loads(output.read_text()),output

    def test_equilibrated_original_body_seed_stays_zero_without_importing_rust_history(self):
        with tempfile.TemporaryDirectory() as tmp:
            result,_=self.run_seed(Path(tmp),100.)
            self.assertEqual(result['frames'][0]['name'],'initial_equilibrated')
            self.assertLess(result['frames'][0]['relative_force_residual'],1e-10)
            self.assertFalse(any(any(row) for row in result['frames'][0]['displacement_m']))

    def test_unbalanced_general_seed_is_rejected_before_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError,'already equilibrated'):self.run_seed(Path(tmp),99.)
            self.assertFalse((Path(tmp)/'reference.json').exists())


if __name__=='__main__':unittest.main()
