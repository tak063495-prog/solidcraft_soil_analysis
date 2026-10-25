"""Preserve native stdout/stderr even where PowerShell transcripts omit them."""
from pathlib import Path
import subprocess
import sys


def main():
    path=Path(sys.argv[1]);command=sys.argv[2:]
    with path.open('x',encoding='utf-8') as log:
        log.write('Command: '+repr(command)+'\n');log.flush()
        process=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,encoding='utf-8',errors='replace')
        for line in process.stdout:
            log.write(line);log.flush();sys.stdout.write(line);sys.stdout.flush()
        status=process.wait();log.write(f'\nExit code: {status}\n')
    return status


if __name__=='__main__':sys.exit(main())
