"""Match held process evidence to one exact staged job; saved PID alone is insufficient."""
from pathlib import Path

def matches_stage_job_process(job,status,owned,records,parse_argv):
    pid=status.get('process_id')
    if status.get('state')!='running' or type(pid) is not int or pid not in owned or pid not in records:
        return False
    record=records[pid];expected=(Path(job).resolve()/'solver/soil-fem-3d.exe').resolve()
    normalize=lambda q:str(Path(q).resolve()).casefold()
    image=record.get('ExecutablePath')
    if not image or not Path(image).is_absolute() or normalize(image)!=normalize(expected):
        return False
    try:args=parse_argv(record.get('CommandLine') or '')
    except (OSError,ValueError,TypeError):return False
    if len(args)!=3 or args[1]!='run-stages-job':return False
    if not Path(args[0]).is_absolute() or normalize(args[0])!=normalize(expected):return False
    if not Path(args[2]).is_absolute() or normalize(args[2])!=normalize(job):return False
    return True
