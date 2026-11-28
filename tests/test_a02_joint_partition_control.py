"""Default identity and bounded diagnostic policy, no heavy histories in CI."""
from pathlib import Path
import sys
import unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'interim/2026-12-21-a02-joint-birth-cause-diagnosis'))
from a02cause_reference import controlled_run
from a02cause_binding import reference,study


class A02JointPartitionControl(unittest.TestCase):
    def test_default_is_the_unchanged_original_function(self):
        self.assertIs(controlled_run(),reference.run)
        self.assertIs(controlled_run(1),reference.run)
        for value in (0,17,True,1.5):
            with self.assertRaises(ValueError):controlled_run(value)

    def test_diagnostic_minimum_changes_partition_observation_only(self):
        values=study.fixture('stress_free',2)
        baseline=[];observed=[]
        a=reference.run(reference.Model(*values),2,lambda r,f:baseline.append((r,f)))
        b=controlled_run(2)(reference.Model(*values),2,lambda r,f:observed.append((r,f)))
        self.assertEqual(a['accepted_states'],b['accepted_states'])
        self.assertTrue(all(min(r['contact_partitions'])>=2 for r,_ in observed))
        for (ra,fa),(rb,fb) in zip(baseline,observed):
            self.assertEqual((ra['stage'],ra['progress']),(rb['stage'],rb['progress']))
            for key in ('generation','joint_active','archived_generation_count'):
                self.assertEqual(fa[key].tolist(),fb[key].tolist())


if __name__=='__main__':unittest.main()
