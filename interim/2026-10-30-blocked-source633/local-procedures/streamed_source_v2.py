"""Experimental bounded blocks for large immutable construction result audits.

Maps remain real dictionaries for existing mechanical validators. Long arrays
are read-only sequences; every element and component remains hash checked.
The limits bound stored/cache accounting, not an auditor's total process RAM.
No production exporter uses this module until the complete real-data gates pass.
"""
from bisect import bisect_right
from collections.abc import Sequence
from pathlib import Path
import hashlib,json,math,operator
from streamed_stage_source import CheckedReader,_Store,ChunkedRun,MAX_COMPONENT,MAX_PAYLOAD,MAX_CACHE_PAYLOAD,sha

MAX_STREAM_SOURCE=32*1024**3
MAX_ITEMS=5_000_000
MAX_COMPONENTS=32768
BLOCK_ITEMS=4096

class Array(Sequence):
 def __init__(self,store,children,blocked=False):
  self.store=store;self.children=children;self.blocked=blocked
  self.ends=[];total=0
  for child in children:
   total+=child['count'] if blocked else 1;self.ends.append(total)
  self.length=total
 def __len__(self):return self.length
 def __iter__(self):
  for child in self.children:
   if self.blocked:yield from self.store.load(child)
   else:yield self.store.value(child)
 def __getitem__(self,key):
  if isinstance(key,slice):return [self[i] for i in range(*key.indices(len(self)))]
  index=operator.index(key)
  if index<0:index+=len(self)
  if index<0 or index>=len(self):raise IndexError('source array index out of range')
  block=bisect_right(self.ends,index)
  if not self.blocked:return self.store.value(self.children[block])
  start=0 if block==0 else self.ends[block-1]
  return self.store.load(self.children[block])[index-start]
 def __eq__(self,other):
  if not isinstance(other,Sequence) or isinstance(other,(str,bytes)) or len(self)!=len(other):return False
  return all(a==b for a,b in zip(self,other,strict=True))

class Store(_Store):
 def value(self,descriptor):
  kind=descriptor['kind']
  if kind=='mapping':return {k:self.value(v) for k,v in descriptor['children'].items()}
  if kind=='sequence':return Array(self,descriptor['children'])
  if kind=='blocked':
   children=descriptor['children']
   if not children:return []
   if len(children)==1:return self.load(children[0])
   return Array(self,children,True)
  return super().value(descriptor)

def index_source(source,folder,*,source_limit=MAX_STREAM_SOURCE,
 component_limit=MAX_COMPONENT,payload_limit=MAX_PAYLOAD,
 cache_payload=MAX_CACHE_PAYLOAD,block_items=BLOCK_ITEMS):
 import ijson
 from ijson.common import JSONError
 source=Path(source).resolve();folder=Path(folder).resolve()
 if not 1<=block_items<=BLOCK_ITEMS or min(source_limit,component_limit,payload_limit,cache_payload)<=0:raise ValueError('invalid source budgets')
 if source.stat().st_size>source_limit:raise ValueError('source exceeds encoded byte budget')
 folder.mkdir(parents=True,exist_ok=False);components=[]
 def save(value,cost,count=None):
  if len(components)>=MAX_COMPONENTS:raise ValueError('too many source components')
  name='component-%05d.json'%len(components);path=folder/name;written=0
  with path.open('xb') as f:
   for piece in json.JSONEncoder(ensure_ascii=False,allow_nan=False,separators=(',',':')).iterencode(value):
    data=piece.encode('utf-8');written+=len(data)
    if written>component_limit:raise ValueError('encoded source component exceeds budget')
    f.write(data)
  d=dict(kind='value',file=name,bytes=written,sha256=sha(path),payload=cost)
  if count is not None:d['count']=count
  components.append(d);return d
 def cost(event,value):return 64+(2*len(value) if event in ('string','map_key') else 0)
 with source.open('rb') as raw:
  reader=CheckedReader(raw,source_limit);events=iter(ijson.basic_parse(reader,use_float=True))
  def next_event():
   try:return next(events)
   except StopIteration as error:raise ValueError('incomplete source JSON') from error
  def scalar(event,value):
   if event not in ('string','number','boolean','null'):raise ValueError('unexpected source token')
   if event=='number' and (not math.isfinite(value) or abs(value)>1e100):raise ValueError('invalid source number')
   return value
  def keep(event,value,depth,used):
   if depth>=128:raise ValueError('source nesting exceeds budget')
   used[0]+=cost(event,value)
   if used[0]>payload_limit:raise ValueError('retained source entry exceeds budget')
   if event=='start_map':
    answer={}
    while True:
     e,v=next_event()
     if e=='end_map':return answer
     if e!='map_key' or len(v)>1024 or len(answer)>=4096 or v in answer:raise ValueError('duplicate or excessive source map key')
     used[0]+=cost(e,v)
     if used[0]>payload_limit:raise ValueError('retained source entry exceeds budget')
     answer[v]=keep(*next_event(),depth+1,used)
   if event=='start_array':
    answer=[]
    while True:
     e,v=next_event()
     if e=='end_array':return answer
     if len(answer)>=MAX_ITEMS:raise ValueError('too many source array entries')
     answer.append(keep(e,v,depth+1,used))
   return scalar(event,value)
  def split(event,value,depth=0,structured_array=False,root=False):
   if depth>=128:raise ValueError('source nesting exceeds budget')
   if event=='start_map':
    children={}
    while True:
     e,v=next_event()
     if e=='end_map':return dict(kind='mapping',children=children)
     if e!='map_key' or len(v)>1024 or len(children)>=(256 if root else 4096) or v in children:raise ValueError('duplicate or excessive source map key')
     limit=(257 if v=='frames' else 512) if root and v in ('frames','transition_audits') else 0
     children[v]=split(*next_event(),depth+1,limit)
   if event=='start_array':
    children=[];block=[];block_cost=0;count=0
    while True:
     e,v=next_event()
     if e=='end_array':break
     count+=1
     if count>(structured_array if structured_array else MAX_ITEMS):raise ValueError('too many source array entries')
     if structured_array:
      children.append(split(e,v,depth+1));continue
     used=[0];entry=keep(e,v,depth+1,used)
     if block and (len(block)>=block_items or block_cost+used[0]>min(payload_limit,32*1024**2)):
      children.append(save(block,block_cost,len(block)));block=[];block_cost=0
     block.append(entry);block_cost+=used[0]
    if structured_array:return dict(kind='sequence',children=children)
    if block:children.append(save(block,block_cost,len(block)))
    return dict(kind='blocked',children=children)
   return save(scalar(event,value),cost(event,value))
  try:
   first=next_event()
   if first!=('start_map',None):raise ValueError('source object required')
   descriptor=split(*first,root=True)
   if next(events,None) is not None:raise ValueError('trailing source JSON')
  except (JSONError,RecursionError) as error:raise ValueError('invalid source JSON: '+str(error)) from error
  digest=reader.digest.hexdigest();count=reader.count
 if sha(source)!=digest or source.stat().st_size!=count:raise ValueError('original source changed while indexing')
 index=dict(schema='soil-fem-source-components/2',complete=True,source_sha256=digest,source_bytes=count,
  source_limit=source_limit,component_limit=component_limit,payload_limit=payload_limit,
  cache_payload=cache_payload,children=descriptor['children'],components=components,block_items=block_items)
 encoded=json.dumps(index,indent=2,allow_nan=False).encode('utf-8')
 if len(encoded)>64*1024**2:raise ValueError('source component index exceeds budget')
 with (folder/'index.json').open('xb') as f:f.write(encoded)
 return ChunkedRun(Store(source,folder,index,cache_payload),index['children'])
