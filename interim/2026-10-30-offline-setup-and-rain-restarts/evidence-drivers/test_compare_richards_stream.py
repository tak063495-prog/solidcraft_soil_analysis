from pathlib import Path
import copy,json,tempfile,unittest
from unittest.mock import patch
from compare_richards_stream import compare,sha

class StreamingRichardsComparison(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.folder=Path(self.temp.name)
        self.report=dict(schema='soil-fem-richards-result/1',complete=True,
            case=dict(settings=dict(elapsed_seconds=2.)),
            frames=[dict(name='rain',state=dict(time_s=100.,relative_cell_head_m=[.1,.2]))],
            attempts=[dict(accepted=True,iterations=[dict(linear=dict(elapsed_seconds=.2))])],
            last_committed=dict(time_s=100.,relative_cell_head_m=[.1,.2]))
    def write(self,name,value):
        q=self.folder/name;q.write_text(json.dumps(value,allow_nan=False),'utf-8');return q
    def check(self,value):
        a=self.write('left.json',self.report);b=self.write('right.json',value)
        return compare(a,b,expected_left_sha256=sha(a))
    def test_only_exact_solver_timer_is_normalized(self):
        changed=copy.deepcopy(self.report);changed['attempts'][0]['iterations'][0]['linear']['elapsed_seconds']=8.
        result=self.check(changed)
        self.assertEqual(result['ignored_solver_timer_values'],1)
        self.assertEqual((result['frames'],result['last_committed_time_s']),(1,100.))
    def test_physical_state_change_rejected(self):
        changed=copy.deepcopy(self.report);changed['last_committed']['relative_cell_head_m'][1]+=.001
        with self.assertRaisesRegex(ValueError,'non-time report field differs'):self.check(changed)
    def test_timer_named_field_elsewhere_is_not_ignored(self):
        changed=copy.deepcopy(self.report);changed['case']['settings']['elapsed_seconds']=3.
        with self.assertRaisesRegex(ValueError,'non-time report field differs'):self.check(changed)
    def test_extra_state_and_truncated_reports_rejected(self):
        for mutation in ('extra','truncated'):
            with self.subTest(mutation=mutation):
                a=self.write('left.json',self.report);b=self.write('right.json',self.report)
                if mutation=='extra':
                    text=b.read_text('utf-8');b.write_text(text[:-1]+',"unexpected":1}','utf-8')
                else:b.write_bytes(b.read_bytes()[:-12])
                with self.assertRaises(ValueError):compare(a,b,expected_left_sha256=sha(a))
    def test_large_arrays_use_streams_and_reference_hash_is_bound(self):
        self.report['last_committed']['relative_cell_head_m']=[.125]*100000
        a=self.write('left.json',self.report);b=self.write('right.json',self.report);expected=sha(a)
        with patch.object(Path,'read_text',side_effect=AssertionError('whole read forbidden')),patch.object(Path,'read_bytes',side_effect=AssertionError('whole read forbidden')):
            result=compare(a,b,expected_left_sha256=expected)
        self.assertGreater(result['numeric_values'],100000)
        with self.assertRaisesRegex(ValueError,'reference hash differs'):compare(a,b,expected_left_sha256='0'*64)
