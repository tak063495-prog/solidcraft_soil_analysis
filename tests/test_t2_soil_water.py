import copy
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
from t1_water_balance import LIMITS, audit


def state(masses):
    return dict(hydraulic=dict(time_s=0, relative_cell_head_m=[None if v is None else 1. for v in masses],
                               cell_storage_offset_m3=[None if v is None else 0. for v in masses]),
                fluid=dict(hydraulic_points=[dict(element_index=i, weight_m3=1., water_equivalent_mass_per_volume=v)
                                              for i, v in enumerate(masses) if v is not None]))


def tanks():
    old, new = state([1., None]), state([.75, .25])
    cells = [dict(element_index=i, before=dict(water_equivalent_volume_m3=1.) if i == 0 else None,
                  after=dict(water_equivalent_volume_m3=v)) for i, v in enumerate([.75,.25])]
    record = dict(name='fill', time_s=0, water=dict(changed_cells=cells, water_added_m3=0., water_removed_m3=0.,
                  internal_transfers=[dict(donor_element=0, receiver_element=1, volume_m3=.25)]))
    frame = dict(phase='material_transition', name='fill', state=new)
    return dict(complete=True, case=dict(one_way=dict(flow=dict(materials=[]))),
                initialization=dict(last_committed=old), frames=[frame, copy.deepcopy(frame)],
                transitions=[record], excavations=[])


class SoilWaterTest(unittest.TestCase):
    def test_required_null_defers_both_profiles(self):
        report = tanks()
        report['frames'][0]['state']['fluid']['hydraulic_points'][0]['water_equivalent_mass_per_volume'] = None
        missing_state = tanks()
        missing_state['frames'][0]['state']['hydraulic'] = None
        missing_endpoint = tanks()
        del missing_endpoint['transitions'][0]['water']['changed_cells'][1]['after']
        for profile in LIMITS:
            self.assertIsNone(audit({}, missing_state, profile)['passed'])
            self.assertIsNone(audit({}, missing_endpoint, profile)['passed'])
            self.assertIsNone(audit({}, report, profile)['passed'])
            self.assertIsNone(audit({}, {'complete':True, 'case':{'stages':[]}}, profile)['passed'])

    def test_internal_pair_recorded_once_during_release(self):
        for profile in LIMITS:
            result = audit({}, tanks(), profile)
            self.assertTrue(result['passed'])
            global_rows = [r for r in result['rows'] if r['kind']=='global']
            self.assertEqual(len(global_rows), 1)
            self.assertEqual([global_rows[0][key] for key in ('P','N','e','D','E','U')], [0.]*6)
            self.assertEqual([r['D'] for r in result['rows'] if r['kind']=='event_pair'], [.25,.25])

    def test_local_leak_and_double_supply_detected(self):
        report = tanks()
        report['transitions'][0]['water']['internal_transfers'][0]['volume_m3'] = .3
        for profile in LIMITS:
            result = audit({}, report, profile)
            self.assertFalse(result['passed'])
            self.assertTrue([r for r in result['rows'] if r['kind']=='global'][0]['passed'])
        report['transitions'][0]['water']['water_added_m3'] = .25
        with self.assertRaises(ValueError):
            audit({}, report, 'verification')
