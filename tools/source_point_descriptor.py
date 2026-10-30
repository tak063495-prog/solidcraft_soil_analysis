"""Explicit raw-point references for a version4 Native view packet.

No GP value is averaged, replaced, or put into this descriptor. Consumers must
hash the complete original result, validate its selected frame and derive raw
point queries from that same authoritative layout/history.
"""
import json
import math
import re
from import_solidcraft_capture import equal

SCHEMA='soil-fem-view-source-points/1'
VIEW_SCHEMA='soil-fem-solidcraft-view/4'
RESULT_SCHEMAS={'soil-fem-projected-mcc-stage-result/1','soil-fem-projected-mcc-stage-result/2'}
COMPONENTS=['xx','yy','zz','xy','yz','xz']

def load(path,*,maximum=8192):
 """Reject ambiguous or oversized JSON before validating a reference contract."""
 def pairs(items):
  result={}
  if len(items)>4096:raise ValueError('source-point JSON object exceeds key budget')
  for key,value in items:
   if key in result:raise ValueError('duplicate source-point JSON key: '+key)
   result[key]=value
  return result
 def reject(value):raise ValueError('non-finite source-point JSON: '+value)
 def real(value):
  number=float(value)
  if not math.isfinite(number):reject(value)
  return number
 with path.open('rb') as stream:raw=stream.read(maximum+1)
 if len(raw)>maximum:raise ValueError('source-point JSON exceeds encoded budget')
 try:value=json.loads(raw.decode('utf-8'),object_pairs_hook=pairs,parse_constant=reject,parse_float=real)
 except RecursionError as error:raise ValueError('source-point JSON nesting exceeds decoder budget') from error
 if not isinstance(value,dict):raise ValueError('source-point JSON object required')
 return value


def make(run,index,result_sha256):
 if run.get('schema') not in RESULT_SCHEMAS or run.get('complete') is not True or run.get('failure_reason') is not None:
  raise ValueError('source-point descriptor requires complete projected MCC')
 if not isinstance(result_sha256,str) or re.fullmatch('[0-9a-f]{64}',result_sha256) is None:
  raise ValueError('invalid authoritative result hash')
 if type(index) is not int or not 0<=index<len(run['frames']):raise ValueError('invalid source frame index')
 frame=run['frames'][index];state=frame['state'];layout=run['integration_points'];active=state['active_elements']
 if not isinstance(frame.get('name'),str) or not frame['name']:raise ValueError('source frame name required')
 if not 1<=len(layout)<=5_000_000 or len(state['material_states'])!=len(layout):raise ValueError('source history/layout count mismatch')
 if not active or any(type(value) is not bool for value in active):raise ValueError('invalid source active mask')
 count=0
 for point in layout:
  owner=point['element_index']
  if type(owner) is not int or not 0<=owner<len(active):raise ValueError('source point owner outside mask')
  count+=active[owner]
 if count==0:raise ValueError('source frame has no active integration points')
 return dict(schema=SCHEMA,source_result_schema=run['schema'],source_result_sha256=result_sha256,
    frame_index=index,frame_name=frame['name'],integration_point_count=len(layout),active_integration_point_count=count,
    length_unit='m',stress_unit='kPa',stress_sign='compression_positive',stress_components=COMPONENTS.copy())

def validate(descriptor,run,index,result_sha256):
 expected=make(run,index,result_sha256)
 if not equal(descriptor,expected):raise ValueError('source-point descriptor differs from authoritative frame')
 return expected
