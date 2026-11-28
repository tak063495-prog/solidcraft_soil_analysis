"""Re-solve legacy inputs; compare committed bytes with the immutable prior index."""
import argparse
import bisect
import filecmp
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import zipfile
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('c01_path_compat_study',HERE/'study.py')
st=importlib.util.module_from_spec(spec);spec.loader.exec_module(st)

def canonical(value):
    if isinstance(value,dict):
        return {k:canonical(v) for k,v in value.items() if k not in ('elapsed_seconds','solver_binary_sha256')}
    if isinstance(value,list): return [canonical(v) for v in value]
    return value

def canonical_bytes(value):
    return json.dumps(canonical(value),ensure_ascii=False,sort_keys=True,separators=(',',':'),allow_nan=False).encode('utf-8')

class Parts(io.RawIOBase):
    """Seek across sealed parts without making a second 600 MB ZIP copy."""
    def __init__(self,paths):
        self.handles=[p.open('rb') for p in paths];self.starts=[0];self.position=0
        for p in paths: self.starts.append(self.starts[-1]+p.stat().st_size)
    def seekable(self): return True
    def tell(self): return self.position
    def seek(self,offset,whence=0):
        position=offset+(self.position if whence==1 else self.starts[-1] if whence==2 else 0)
        if position<0: raise ValueError('negative part seek')
        self.position=position;return position
    def read(self,size=-1):
        size=self.starts[-1]-self.position if size<0 else min(size,self.starts[-1]-self.position)
        chunks=[]
        while size>0:
            i=bisect.bisect_right(self.starts,self.position)-1
            n=min(size,self.starts[i+1]-self.position);f=self.handles[i]
            f.seek(self.position-self.starts[i]);data=f.read(n)
            if len(data)!=n: raise ValueError('truncated sealed part')
            chunks.append(data);self.position+=n;size-=n
        return b''.join(chunks)
    def close(self):
        for f in self.handles: f.close()
        super().close()

def compare_prior_runs(out,report=None):
    prior=st.ROOT/'interim/2026-12-11-c01-global-step-control'
    index=st.read(prior/'archive-index.json')
    if report is None: report=st.read(out/'compatibility.json')
    assert report['prior_index_sha256']==st.sha(prior/'archive-index.json')
    paths=[]
    for row in index['parts']:
        p=prior/row['file'];assert p.stat().st_size==row['bytes'] and st.sha(p)==row['sha256'];paths.append(p)
    with Parts(paths) as parts:
        h=hashlib.sha256()
        while chunk:=parts.read(1024*1024): h.update(chunk)
        assert h.hexdigest()==index['sha256']
        parts.seek(0)
        entries={e['path']:e for e in index['entries']}
        with zipfile.ZipFile(parts) as z:
            for row in report['rows']:
                if row['mode']!='switch': continue
                n=row['steps'];name=f'paths/g{n}/native-raw/run.json';data=z.read(name)
                assert hashlib.sha256(data).hexdigest()==entries[name]['sha256']
                target=out/'prior-run'/f'g{n}'/'run.json';target.parent.mkdir(parents=True,exist_ok=True)
                if target.exists(): assert target.read_bytes()==data
                else: target.write_bytes(data)
                assert canonical_bytes(json.loads(data))==canonical_bytes(st.read(out/'switch'/f'g{n}'/'raw/run.json'))
                row['previous_switch_run_exact_except_clock_and_executable']=True
    report.update(prior_raw_all_parts_verified=True,prior_archive_entries=len(index['entries']),prior_archive_parts=len(index['parts']))
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--solver',type=Path,required=True)
    p.add_argument('--baseline-default-raw',type=Path,required=True)
    p.add_argument('--prior-index',type=Path,default=st.ROOT/'interim/2026-12-11-c01-global-step-control/archive-index.json')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); out=a.output.resolve();st.initialize(out,a.solver.resolve())
    index=st.read(a.prior_index); entries={e['path']:e for e in index['entries']}
    case=st.read(out/'stages.json');rows=[]
    binary=out/('probe.exe' if (out/'probe.exe').exists() else 'probe')
    for mode,levels in [('default',(16,)),('switch',st.STEPS)]:
        for n in levels:
            folder=out/mode/f'g{n}';folder.mkdir(parents=True)
            cfg=st.prior.configuration(case,n,5e-5)
            if mode=='switch':
                cfg['joint_system']['global_step_control']={'switch_bisections':4}
                cfg['settings']['minimum_increment']=min(cfg['settings']['minimum_increment'],1/(n*32))
            st.write(folder/'stages.json',cfg)
            cmd=[str(binary),str(out/'mesh-h0.json'),str(out/'ground.json'),str(folder/'stages.json'),str(folder/'raw')]
            with (folder/'process.log').open('x') as log:
                result=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=1800)
            assert result.returncode==0,(mode,n,result.returncode)
            raw=folder/'raw';run=st.read(raw/'run.json')
            assert run['complete']
            digest=st.sha(raw/'accepted.jsonl')
            if mode=='default':
                baseline=a.baseline_default_raw.resolve()
                assert filecmp.cmp(raw/'accepted.jsonl',baseline/'accepted.jsonl',shallow=False)
                assert canonical_bytes(run)==canonical_bytes(st.read(baseline/'run.json'))
                expected=st.sha(baseline/'accepted.jsonl')
                old=st.read(baseline/'run.json')['provenance']
                assert {k:v for k,v in run['provenance'].items() if k!='solver_binary_sha256'}=={k:v for k,v in old.items() if k!='solver_binary_sha256'}
            else:
                expected=entries[f'paths/g{n}/native-raw/accepted.jsonl']['sha256']
                assert digest==expected,('prior switch bytes changed',n,digest,expected)
            rows.append(dict(mode=mode,steps=n,complete=True,accepted_states=sum(1 for _ in (raw/'accepted.jsonl').open()),
                             accepted_bytes_exact=True,accepted_sha256=digest,expected_sha256=expected,
                             provenance=run['provenance'],default_run_exact_except_clock_and_executable=(True if mode=='default' else None)))
            print('legacy byte comparison complete',mode,n,flush=True)
    report=dict(starting_commit=st.BASE,legacy_default_unchanged=True,
             previous_switch_unchanged=True,prior_index_sha256=st.sha(a.prior_index),
             old_checkpoint_executable_identity_rewritten=False,rows=rows)
    st.write(out/'compatibility.json',compare_prior_runs(out,report))

if __name__=='__main__': main()
