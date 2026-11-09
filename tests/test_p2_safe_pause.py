"""Pure snapshot data only: this test never opens or terminates a process."""
import ast,copy,datetime,unittest
from pathlib import Path
PATH=Path(__file__).resolve().parents[1]/'interim/2026-10-30-pause19-durable-SRM-mixed-elements/pause/scripts/safe-pause-latest102.py'
# Avoid importing the interim Windows entry point, which inventories real CIM processes.
tree=ast.parse(PATH.read_text());namespace={'datetime':datetime}
exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ('process_identity','valid_edge','expand_owned')],type_ignores=[]),str(PATH),'exec'),namespace)
expand=namespace['expand_owned']
def record(pid,parent,created):return dict(ProcessId=pid,ParentProcessId=parent,Created=created,ExecutablePath='C:/isolated-fixture/python.exe')
class Ownership(unittest.TestCase):
    def test_reused_parent_pid_does_not_claim_older_child(self):
        records={100:record(100,1,'2026-10-30T23:00:00+00:00'),200:record(200,100,'2026-10-30T22:00:00+00:00')}
        self.assertEqual(expand(records,{100},set()),{100})
    def test_ordered_children_grandchildren_and_excluded_are_respected(self):
        records={100:record(100,1,'2026-10-30T22:00:00+00:00'),200:record(200,100,'2026-10-30T22:01:00+00:00'),300:record(300,200,'2026-10-30T22:02:00+00:00'),400:record(400,100,'2026-10-30T22:03:00+00:00')}
        self.assertEqual(expand(records,{100},{400}),{100,200,300})
        self.assertEqual(expand(records,{100},{200}),{100,400})
    def test_unknown_missing_unparseable_and_identity_mismatch_fail_closed(self):
        base={100:record(100,1,'2026-10-30T22:00:00+00:00'),200:record(200,100,'2026-10-30T22:01:00+00:00')}
        for pid in (100,200):
            for key,value in (('Created',None),('Created','unknown'),('Created','2026-10-30T22:01:00'),('ProcessId',999),('ExecutablePath',None),('ExecutablePath','')):
                with self.subTest(pid=pid,key=key,value=value):
                    records=copy.deepcopy(base);records[pid][key]=value
                    self.assertEqual(expand(records,{100},set()),set() if pid==100 else {100})
            records=copy.deepcopy(base);records[pid].pop('Created');self.assertEqual(expand(records,{100},set()),set() if pid==100 else {100})
        self.assertEqual(expand({}, {100},set()),set())
    def test_live_handle_identity_guards_are_retained(self):
        source=PATH.read_text();self.assertIn('kernel.QueryFullProcessImageNameW(h,0,image,ctypes.byref(size))',source)
        self.assertIn("norm(image.value)==norm(r['ExecutablePath'])",source)
        self.assertIn('kernel.GetProcessTimes(h,',source);self.assertIn("abs(actual-expected)<10000,'PID was reused'",source)
if __name__=='__main__':unittest.main()
