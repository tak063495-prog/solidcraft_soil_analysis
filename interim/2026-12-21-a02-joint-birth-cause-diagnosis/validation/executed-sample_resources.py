import pathlib,time,json
pid=3216
out=pathlib.Path('/workspace/solidcraft_soil_analysis/work/a02-resource-samples.jsonl')
with out.open('x') as f:
 while pathlib.Path(f'/proc/{pid}').exists():
  selected={pid};info={}
  for path in pathlib.Path('/proc').glob('[0-9]*/status'):
   try:
    data={k:v.strip() for k,v in [line.split(':',1) for line in path.read_text().splitlines() if ':' in line]}
    info[int(path.parent.name)]=data
   except (OSError,ValueError):pass
  for _ in range(5):selected|={p for p,d in info.items() if int(d.get('PPid','0')) in selected}
  rows=[dict(pid=p,name=info[p].get('Name'),rss_kib=int(info[p].get('VmRSS','0 kB').split()[0]),hwm_kib=int(info[p].get('VmHWM','0 kB').split()[0])) for p in sorted(selected) if p in info]
  f.write(json.dumps(dict(monotonic=time.monotonic(),processes=rows,total_rss_kib=sum(x['rss_kib'] for x in rows)))+'\n');f.flush();time.sleep(1)
