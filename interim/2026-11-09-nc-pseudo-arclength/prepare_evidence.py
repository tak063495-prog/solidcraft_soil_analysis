"""Verify the prior archive and reproduce the unchanged native observation crate."""
import argparse,hashlib,json,pathlib,shutil,subprocess,tarfile,time
parser=argparse.ArgumentParser();parser.add_argument('--output',type=pathlib.Path,required=True);parser.add_argument('--source-only',action='store_true');args=parser.parse_args()
root=pathlib.Path(__file__).resolve().parents[2];w=args.output.resolve();w.mkdir(parents=True,exist_ok=False);prior=root/'interim/2026-11-09-nc-turning-path-diagnosis'
base='e21778c1e0aeeb28181322808dd0800a9747a507'
assert subprocess.check_output(['git','rev-parse',base+'^'],cwd=root).decode().strip()=='06bff71dcc54f7060932aff57194f04ffe0593f8'
assert not subprocess.check_output(['git','diff',base,'--','src','Cargo.toml','Cargo.lock','rust-toolchain.toml'],cwd=root)

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
if not args.source_only:
    manifest=json.loads((prior/'archive-sha256.json').read_text())['archives']['evidence.tar.xz']
    assert sha(prior/'evidence.tar.xz')==manifest['archive_sha256']
    t=time.monotonic();count=0
    with tarfile.open(prior/'evidence.tar.xz','r:xz') as z:
        assert set(z.getnames())==set(manifest['entries'])
        for m in z.getmembers():
            assert m.isfile(); b=z.extractfile(m).read();assert hashlib.sha256(b).hexdigest()==manifest['entries'][m.name];count+=1
            if m.name in ('field-records.json','assembled-matrices.json') or m.name.startswith('prior/fresh-head-midpoint/input/') or m.name=='prior/fresh-head-midpoint/job/recovery-checkpoint-0.json':
                p=w/'prior'/m.name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b)
    (w/'prior-binding.json').write_text(json.dumps(dict(archive_sha256=manifest['archive_sha256'],verified_entries=count,elapsed_seconds=time.monotonic()-t),indent=2))
s=w/'isolated';s.mkdir();shutil.copytree(root/'src',s/'src')
for n in ('Cargo.toml','Cargo.lock','rust-toolchain.toml'):shutil.copyfile(root/n,s/n)
for n,a in [('mcc.rs','mcc-diagnostic-addition.rs'),('mcc_midpoint.rs','midpoint-diagnostic-addition.rs')]:
    p=s/'src'/n;p.write_bytes(p.read_bytes().replace(b'\r\n',b'\n')+(prior/a).read_bytes())
(s/'examples').mkdir();shutil.copyfile(prior/'material_diagnostic.rs',s/'examples/material_diagnostic.rs')
subprocess.run(['cargo','fmt','--manifest-path',str(s/'Cargo.toml')],check=True)
(w/'source-sha256.json').write_text(json.dumps({str(p.relative_to(s)):sha(p) for p in sorted(s.rglob('*')) if p.is_file()},indent=2))
if args.source_only:raise SystemExit(0)
subprocess.run(['cargo','build','--locked','--offline','--release','--manifest-path',str(s/'Cargo.toml'),'--target-dir',str(w/'target'),'--example','material_diagnostic'],check=True)
(w/'binary-sha256.json').write_text(json.dumps({'material_diagnostic':sha(w/'target/release/examples/material_diagnostic')},indent=2))
