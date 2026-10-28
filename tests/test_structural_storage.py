import copy,json
from pathlib import Path
import sys,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from plate_surfaces import attach
import test_plate_surfaces as fixtures
from verify_mcc_structural_study import physical_case

class StructuralStorage(unittest.TestCase):
    def test_opt_in_budget_preserves_every_geometry_attachment_and_section(self):
        mesh,case,spec=fixtures.PlateSurfaces().fixture();a,b,bind=attach(mesh,case,spec)
        self.assertNotIn('storage_budget',b['beam_system'])
        proposed=copy.deepcopy(spec);proposed['storage_budget']=dict(maximum_nodes=750,maximum_compilation_mib=1024,maximum_cached_mib=512)
        x,y,other=attach(mesh,case,proposed)
        self.assertEqual(y['beam_system'].pop('storage_budget'),proposed['storage_budget'])
        self.assertEqual((x,y,other),(a,b,bind))

    def test_invalid_conflicting_and_small_budgets_reject_without_input_mutation(self):
        good=dict(maximum_nodes=750,maximum_compilation_mib=1024,maximum_cached_mib=512)
        variants=[None,{**good,'maximum_nodes':True},{**good,'maximum_nodes':1001},{**good,'maximum_compilation_mib':0},{**good,'maximum_compilation_mib':2049},{**good,'maximum_cached_mib':1025},{**good,'unknown':1},{**good,'maximum_nodes':1}]
        for budget in variants:
            mesh,case,spec=fixtures.PlateSurfaces().fixture();spec['storage_budget']=budget;before=copy.deepcopy((mesh,case,spec))
            with self.assertRaises(ValueError):attach(mesh,case,spec)
            self.assertEqual((mesh,case,spec),before)

    def test_storage_is_numerical_while_section_material_remains_physical(self):
        root=Path(__file__).resolve().parents[1];case=json.loads((root/'examples/stages-beam-slope.json').read_text(encoding='utf-8'))
        other=copy.deepcopy(case);other['beam_system']['storage_budget']=dict(maximum_nodes=750,maximum_compilation_mib=1024,maximum_cached_mib=512)
        self.assertEqual(physical_case(case),physical_case(other))
        other['beam_system']['sections'][0]['young_kpa']*=2
        self.assertNotEqual(physical_case(case),physical_case(other))

if __name__=='__main__':unittest.main()
