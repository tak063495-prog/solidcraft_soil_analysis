"""Disk-backed JSON components for complete, source-bound construction audits.

This changes storage only. Existing mechanical audits still see every frame,
geostatic seed, transition and frozen material history through Mapping/Sequence.
The512MiB component/cache accounting is separate from an auditor's NumPy arrays
and references it retains. It is not a512MiB whole-process memory guarantee.
"""
from collections import OrderedDict
from collections.abc import Mapping,Sequence
from pathlib import Path
import hashlib,json,math,operator,re

MAX_SOURCE=2*1024**3
MAX_COMPONENT=256*1024**2
MAX_PAYLOAD=512*1024**2
MAX_CACHE_PAYLOAD=512*1024**2
SPLIT_ARRAYS={'frames','transition_audits'}
SPECIAL=re.compile(rb'["\\]')
NUMBER=re.compile(rb'[-+0-9.eE]+')

def sha(path):
 with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

class CheckedReader:
 def __init__(self,file,maximum):
  self.file=file;self.maximum=maximum;self.count=0;self.digest=hashlib.sha256()
  self.inside=False;self.escape=False;self.string_bytes=0;self.number_bytes=0
 def read(self,n=-1):
  if n<0:n=65536
  data=self.file.read(min(n,65536));self.count+=len(data)
  if self.count>self.maximum:raise ValueError('source exceeds encoded byte budget')
  self.digest.update(data);cursor=0
  while cursor<len(data):
   if not self.inside:
    quote=data.find(b'"',cursor);end=len(data) if quote<0 else quote
    segment=data[cursor:end]
    for match in NUMBER.finditer(segment):
     length=len(match.group())+(self.number_bytes if match.start()==0 else 0)
     if length>4096:raise ValueError('source numeric token exceeds budget')
     self.number_bytes=length if match.end()==len(segment) else 0
    if segment and not NUMBER.fullmatch(segment[-1:]):self.number_bytes=0
    if quote<0:break
    self.inside=True;self.string_bytes=0;self.number_bytes=0;cursor=quote+1
   else:
    if self.escape:
     self.string_bytes+=1;self.escape=False;cursor+=1
    else:
     match=SPECIAL.search(data,cursor);end=len(data) if match is None else match.start()
     self.string_bytes+=end-cursor;cursor=end
     if match is None:
      if self.string_bytes>1024**2:raise ValueError('source string token exceeds budget')
      break
     if data[cursor]==34:self.inside=False
     else:self.escape=True
     cursor+=1
    if self.string_bytes>1024**2:raise ValueError('source string token exceeds budget')
  return data

class _Store:
 def __init__(self,source,folder,index,cache_payload):
  self.source=source;self.folder=folder;self.index=index
  self.cache=OrderedDict();self.cache_payload=cache_payload;self.retained=0
 def load(self,descriptor):
  name=descriptor['file']
  if name in self.cache:
   value,cost=self.cache.pop(name);self.cache[name]=(value,cost);return value
  path=self.folder/name
  if path.resolve().parent!=self.folder.resolve():raise ValueError('component escapes cache')
  if path.stat().st_size!=descriptor['bytes']:raise ValueError('source component changed')
  cost=descriptor['payload']
  while self.cache and self.retained+cost>self.cache_payload:
   _,(_,used)=self.cache.popitem(last=False);self.retained-=used
  with path.open('rb') as f:data=f.read(self.index['component_limit']+1)
  if len(data)!=descriptor['bytes'] or hashlib.sha256(data).hexdigest()!=descriptor['sha256']:raise ValueError('source component changed')
  value=json.loads(data)
  if cost<=self.cache_payload:self.cache[name]=(value,cost);self.retained+=cost
  return value
 def value(self,descriptor):
  kind=descriptor['kind']
  if kind=='value':return self.load(descriptor)
  if kind=='mapping':return ChunkedRun(self,descriptor['children'])
  if kind=='sequence':return ChunkedList(self,descriptor['children'])
  raise ValueError('invalid component descriptor')
 def clear(self):self.cache.clear();self.retained=0
 def verify(self):
  if sha(self.source)!=self.index['source_sha256']:raise ValueError('original source changed')
  for descriptor in self.index['components']:
   path=self.folder/descriptor['file']
   if path.stat().st_size!=descriptor['bytes'] or sha(path)!=descriptor['sha256']:raise ValueError('source component changed')

class ChunkedRun(Mapping):
 def __init__(self,store,children):self._store=store;self._children=children
 def __len__(self):return len(self._children)
 def __iter__(self):return iter(self._children)
 def __getitem__(self,key):return self._store.value(self._children[key])
 def clear_cache(self):self._store.clear()
 def verify_source(self):self._store.verify()

class ChunkedList(Sequence):
 def __init__(self,store,children):self._store=store;self._children=children
 def __len__(self):return len(self._children)
 def __getitem__(self,key):
  if isinstance(key,slice):return ChunkedList(self._store,self._children[key])
  return self._store.value(self._children[operator.index(key)])

def index_source(source,folder,*,source_limit=MAX_SOURCE,component_limit=MAX_COMPONENT,payload_limit=MAX_PAYLOAD,cache_payload=MAX_CACHE_PAYLOAD):
 import ijson
 from ijson.common import ObjectBuilder,JSONError
 source=Path(source).resolve();folder=Path(folder).resolve()
 if min(source_limit,component_limit,payload_limit,cache_payload)<=0:raise ValueError('invalid source budget')
 if source.stat().st_size>source_limit:raise ValueError('source exceeds encoded byte budget')
 folder.mkdir(parents=True,exist_ok=False)
 children={};components=[];builder=None;target=None;prefix_root=None;payload=0;stack=[];root_seen=False;root_pending=None;geo_pending=None
 def save(value,cost):
  if len(components)>=2048:raise ValueError('too many source components')
  name='component-%04d.json'%len(components);path=folder/name
  written=0
  with path.open('xb') as f:
   for piece in json.JSONEncoder(ensure_ascii=False,allow_nan=False,separators=(',',':')).iterencode(value):
    data=piece.encode('utf-8');written+=len(data)
    if written>component_limit:raise ValueError('encoded source component exceeds budget')
    f.write(data)
  descriptor=dict(kind='value',file=name,bytes=path.stat().st_size,sha256=sha(path),payload=cost)
  components.append(descriptor);return descriptor
 with source.open('rb') as raw:
  reader=CheckedReader(raw,source_limit)
  def events():
   try:yield from ijson.parse(reader,use_float=True)
   except JSONError as error:raise ValueError('invalid source JSON: '+str(error)) from error
  for prefix,event,value in events():
   if not root_seen:
    if (prefix,event)!=('','start_map'):raise ValueError('source object required')
    root_seen=True
   if event=='map_key':
    if len(value)>1024 or len(stack[-1])>=4096 or value in stack[-1]:raise ValueError('duplicate or excessive source map key')
    stack[-1].add(value)
   elif event in ('start_map','start_array'):
    if len(stack)>=128:raise ValueError('source nesting exceeds budget')
    stack.append(set() if event=='start_map' else None)
   elif event in ('end_map','end_array'):stack.pop()
   elif event in ('number','double','integer') and (not math.isfinite(value) or abs(value)>1e100):raise ValueError('invalid source number')
   if prefix=='' and event=='map_key' and len(stack)==1:
    if len(children)>=256:raise ValueError('too many top-level source fields')
    root_pending=value
   if prefix=='geostatic' and event=='map_key' and len(stack)==2 and builder is None:geo_pending=value
   if builder is None:
    # Prefix strings alone are ambiguous when a literal map key contains a
    # dot, "item", or an empty string. Use the actual direct map key first.
    if root_pending is not None and prefix==root_pending and event!='map_key':
     key=root_pending;root_pending=None
     if key in SPLIT_ARRAYS and event=='start_array':
      children[key]=dict(kind='sequence',children=[]);continue
     if key=='geostatic' and event=='start_map':
      children[key]=dict(kind='mapping',children={});continue
     target=(children,key)
    elif geo_pending is not None and prefix=='geostatic.'+geo_pending and event!='map_key':
     target=(children['geostatic']['children'],geo_pending);geo_pending=None
    elif prefix in SPLIT_ARRAYS and event=='end_array' or prefix=='geostatic' and event=='end_map':continue
    elif prefix.endswith('.item') and prefix[:-5] in SPLIT_ARRAYS:
     key=prefix[:-5];sequence=children[key]['children']
     if len(sequence)>=(257 if key=='frames' else 512):raise ValueError('too many source frames/audits')
     target=sequence
    else:continue
    builder=ObjectBuilder();prefix_root=prefix;payload=0
   if event not in ('end_map','end_array'):
    payload+=64+(2*len(value) if event in ('string','map_key') else 0)
    if payload>payload_limit:raise ValueError('retained source component exceeds budget')
   builder.event(event,value)
   if prefix==prefix_root and event not in ('start_map','start_array','map_key'):
    descriptor=save(builder.value,payload)
    if isinstance(target,tuple):target[0][target[1]]=descriptor
    else:target.append(descriptor)
    builder=None
  if builder is not None or stack:raise ValueError('incomplete component source')
  digest=reader.digest.hexdigest();count=reader.count
 if sha(source)!=digest or source.stat().st_size!=count:raise ValueError('original source changed while indexing')
 index=dict(schema='soil-fem-source-components/1',complete=True,source_sha256=digest,source_bytes=count,
    component_limit=component_limit,payload_limit=payload_limit,cache_payload=cache_payload,children=children,components=components)
 with (folder/'index.json').open('x',encoding='utf-8') as f:json.dump(index,f,indent=2,allow_nan=False)
 return ChunkedRun(_Store(source,folder,index,cache_payload),children)
