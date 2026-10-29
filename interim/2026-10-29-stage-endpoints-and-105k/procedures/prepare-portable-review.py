"""Portable data/theory review of the completed105456 nonlinear MCC run."""
from pathlib import Path
import datetime
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
import zipfile

p = Path(__file__).resolve().parent
root = p.parents[1]
source = root/'runs/mcc-parallel-100k-durable-20261029-v1'
case = source/'n26'
package = p/'review-package'
read = lambda q: json.loads(q.read_text('utf-8-sig'))


def sha(q):
    h = hashlib.sha256()
    with q.open('rb') as f:
        for block in iter(lambda:f.read(4<<20),b''):
            h.update(block)
    return h.hexdigest()


def write(path,value):
    with path.open('x',encoding='utf-8') as f:
        json.dump(value,f,indent=2,allow_nan=False)


def copy(a,b):
    if (p/'pause-requested').exists():
        raise SystemExit(2)
    b.parent.mkdir(parents=True,exist_ok=True)
    digest = sha(a)
    shutil.copy2(a,b)
    assert sha(a)==sha(b)==digest


write(p/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
audit = read(case/'analytical-verification.json')
resources = read(case/'process-terminal.json')
assert read(case/'terminal.json')['passed'] and audit['passed']
assert audit['tetrahedra']==105456 and audit['all_nodes_and_integration_points_checked']
assert resources['exit_code']==0 and resources['held_handle_creation_and_image_verified']
assert resources['solver_sha256']=='d55ffe42c3b61e0891fec513489c9f5881830733ac9e003a01a0cf69f495905d'
package.mkdir()
for name in ['mesh.json','ground.json','stages.json','preparation.json']:
    copy(case/name,package/'n26'/name)
copy(case/'solver-job/result.json',package/'n26/solver-job/result.json')
copy(source/'verify.py',package/'verify.py')
for name in ['analytical-verification.json','process-terminal.json','process-observations.jsonl','terminal.json']:
    copy(case/name,package/'recorded-reports'/name)
copy(source/'original-input-equivalence.json',package/'recorded-reports/original-input-equivalence.json')
frozen = root/'runs/richards-local-rain-3d-20261029-v1'
for name in ['source-snapshot.json','source-snapshot.zip']:
    copy(frozen/name,package/'solver-sources'/name)
assert len(read(package/'solver-sources/source-snapshot.json'))==611
assert sha(package/'n26/solver-job/result.json')==audit['result_sha256']
text = '''# 105456要素の非線形MCC：全場・理論照合用データ

本体3870c2fの数値ソースに対応する、最大4スレッドの材料更新を使った一様静水圧パッチ試験です。105456 Tet4・19683節点・105456積分点の初期平衡、塑性載荷、除荷を含みます。重力は既知の節点荷重で相殺しています。CAD形状、掘削、水理、構造部材を含む10万要素解析の受入ではありません。

元の全結果JSON、未解析入力、独立した理論・幾何監査、資源記録、当該数値ソースを含みます。実行exeと保存点は含めません。保存点を使う場合は、それを生成した実行exeの完全一致が必要です。

展開後、NumPyを用意したPythonで次を実行できます。元の記録はrecorded-reportsに保持し、新しい監査はn26/analytical-verification.jsonへ作成します。

```powershell
python verify.py 26
```

全節点・全積分点の変位、応力、pc、間隙比、塑性ひずみ、独立した体積・境界面積・相殺荷重と、載荷中の棄却復元を検査します。この操作は保存済み全場の理論監査で、非線形求解の再実行ではありません。新しい実行exeで再計算するには、solver-sourcesを展開してビルドし、n26の入力をprepare-durable-stages-jobへ指定します。新exeで生成した保存点を使用してください。

経過6543.93秒、CPU23261.39秒、観測peak working set660246528 bytesです。他の検証を同時実行した同一PCの観測であり、専有環境の速度保証や直列版との速度比ではありません。除荷後の最大変位差5.984e-12 m未満、応力差2.098e-11 kPa未満です。時間・履歴の全場照合と新規配置での再監査の記録を区別します。P0～P6全体は未完です。
'''
(package/'README.ja.md').write_text(text,'utf-8',newline='\n')
write(package/'metadata.json',dict(schema='soil-fem-105k-portable-review/1',passed_source_audit=True,
    full_result_included=True,full_result_bytes=(case/'solver-job/result.json').stat().st_size,
    full_result_sha256=audit['result_sha256'],tetrahedra=105456,nodes=19683,integration_points=105456,
    solver_sha256=resources['solver_sha256'],numerical_source_commit='3870c2fdea261b5c4c6a740a283f147db81848c7',
    solver_source_files=611,resources=resources,independent_audit_replay_does_not_rerun_nonlinear_solver=True,
    loaded_snapshot_recovery_still_separate=True,new_pc_or_gpu_not_accepted=True,p6_complete=False))
files = {q.relative_to(package).as_posix():sha(q) for q in package.rglob('*') if q.is_file()}
write(package/'files-sha256.json',files)
archive = p/'mcc-105456-full-field-review.zip'
with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for q in sorted(package.rglob('*')):
        if q.is_file():
            assert q.suffix.lower() not in ['.exe','.dll','.pdb']
            z.write(q,q.relative_to(package).as_posix())
assert archive.stat().st_size<95_000_000
extracted = p/'Portable105456 review'
extracted.mkdir()
with zipfile.ZipFile(archive) as z:
    z.extractall(extracted)
assert {name:sha(extracted/name) for name in files}==files
assert not (extracted/'n26/analytical-verification.json').exists()
if (p/'pause-requested').exists():
    raise SystemExit(2)
started = time.monotonic()
env = dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1')
with (p/'replay.log').open('x',encoding='utf-8') as log:
    child = subprocess.Popen([sys.executable,'-X','utf8',str(extracted/'verify.py'),'26'],cwd=extracted,env=env,
        stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
    write(p/'replay-child.json',dict(pid=child.pid,parent_pid=os.getpid(),executable=sys.executable,script=str(extracted/'verify.py'),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    code = child.wait()
write(p/'replay-terminal.json',dict(exit_code=code,elapsed_seconds=time.monotonic()-started))
assert code==0
actual = read(extracted/'n26/analytical-verification.json')
assert actual==audit
assert files=={name:sha(package/name) for name in files}
assert sha(case/'solver-job/result.json')==audit['result_sha256']
write(p/'terminal.json',dict(passed=True,archive=archive.name,archive_bytes=archive.stat().st_size,archive_sha256=sha(archive),
    full105456_result_included=True,full_data_file_count=len(files),portable_extracted_all_source_and_data_hashes_exact=True,
    actual_fresh_folder_with_spaces_independent_all_node_point_theory_audit_exact=True,
    no_solver_execution_or_new_pc_acceptance_inferred=True,p6_complete=False))
print('PASS105456 portable full-field archive and exact fresh-folder theoretical re-audit',archive.stat().st_size,flush=True)
