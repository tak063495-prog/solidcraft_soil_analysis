"""Verify the complete unsolved input archive before creating a fresh output folder."""
from pathlib import Path,PurePosixPath
import argparse,hashlib,json,zipfile

def unpack(output):
 base=Path(__file__).resolve().parent
 manifest=json.loads((base/'inputs-sha256.json').read_text(encoding='utf-8'));archive=base/'inputs.zip'
 digest=lambda data:hashlib.sha256(data).hexdigest()
 if digest(archive.read_bytes())!=manifest['zip_sha256']:raise ValueError('input archive SHA256 differs')
 files=manifest['files'];verified={}
 with zipfile.ZipFile(archive) as z:
  names=z.namelist()
  if len(names)!=len(set(names)) or set(names)!=set(files):raise ValueError('archive inventory differs')
  for name in names:
   q=PurePosixPath(name);info=z.getinfo(name)
   if q.is_absolute() or '..' in q.parts or '\\' in name or ':' in name or (info.external_attr>>16)&0o170000==0o120000:raise ValueError('unsafe archive path')
   if info.file_size!=files[name]['bytes']:raise ValueError('archive size differs')
   data=z.read(name)
   if digest(data)!=files[name]['sha256']:raise ValueError('archive member SHA256 differs')
   verified[name]=data
 output=output.resolve()
 if output.exists():raise FileExistsError('output must be a new folder')
 for name in verified:
  if not (output/name).resolve().is_relative_to(output):raise ValueError('output escapes destination')
 output.mkdir(parents=True)
 for name,data in verified.items():
  target=output/name;target.parent.mkdir(parents=True,exist_ok=True)
  with target.open('xb') as f:f.write(data)
 print(json.dumps(dict(verified=True,unsolved_inputs=len(verified),output=str(output))))

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('output',type=Path);unpack(parser.parse_args().output)
