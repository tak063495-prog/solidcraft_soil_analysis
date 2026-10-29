"""Check whole-source semantic preservation and the previously independent MCC values."""
from pathlib import Path
from collections.abc import Mapping,Sequence
import datetime,hashlib,json,struct,sys,time
p=Path(__file__).resolve().parent
sys.path[:0]=[str(p/'runtime/tools'),str(p.parent/'native-typed-source-development-20261029-v1/deps')]
import ijson
from streamed_stage_source import index_source,sha
read=lambda q:json.loads(q.read_text('utf-8'))
source=p.parent/'mcc-next-fixed-wall-durable-pair-20261029-v1/mesh8061/workflow/analysis/stage-results.json'
expected=next(v for k,v in read(p.parent/'mcc-next-fixed-wall-durable-pair-20261029-v1/mesh8061/field-verification.json')['sha256'].items() if k.endswith('analysis\\stage-results.json'))
out=p/'producer-source-index-v2';out.mkdir(exist_ok=True)
assert sha(source)==expected
start=time.monotonic();run=index_source(source,out/'components')
print(json.dumps(dict(phase='index_complete',elapsed_seconds=time.monotonic()-start)),flush=True)
reference=p.parent/'native-typed-source-development-20261029-v1/real8061-v4'
rows=[]
for i,frame in enumerate(run['frames']):
 report=read(reference/f'frame{i}.json');retained=read(reference/f'frame{i}.retained.json')
 actual=dict(frame);actual['state']=dict(frame['state']);states=actual['state'].pop('material_states')
 assert actual==retained['frames'][i] and run['case']==retained['case']
 digest=hashlib.sha256()
 for s in states:
  digest.update(struct.pack('<14d',*s['effective_stress_compression_positive_kpa'],s['void_ratio'],s['preconsolidation_pressure_kpa'],*s['plastic_strain_compression_positive']))
 assert len(states)==217647 and digest.hexdigest()==report['material_f64_le_sha256']
 rows.append(dict(frame=i,points=len(states),all_MCC_values_and_other_frame_fields_equal=True))
 print(json.dumps(dict(phase='frame_checked',frame=i)),flush=True)
del frame,states,actual;run.clear_cache()
digest=hashlib.sha256();layout=run['integration_points']
for q in layout:digest.update(struct.pack('<QQ17d',q['element_index'],q['local_index'],*q['natural_coordinates'],*q['position_m'],q['weight_m3'],*q['shape_functions']))
assert len(layout)==217647 and digest.hexdigest()==report['layout_u64_f64_le_sha256']
del layout;run.clear_cache()
def add(h,event,value):
 if event=='map_key' or event=='string':
  raw=value.encode('utf-8');h.update(b'K' if event=='map_key' else b'S');h.update(struct.pack('<Q',len(raw)));h.update(raw)
 elif event in ('integer','double','number'):
  if type(value) is int:
   raw=str(value).encode('ascii');h.update(b'I'+struct.pack('<Q',len(raw))+raw)
  else:h.update(b'D'+struct.pack('<d',value))
 elif event=='boolean':h.update(b'T' if value else b'F')
 else:h.update({'start_map':b'{','end_map':b'}','start_array':b'[','end_array':b']','null':b'N'}[event])
raw_hash=hashlib.sha256();raw_events=0
with source.open('rb') as f:
 for _,event,value in ijson.parse(f,use_float=True):add(raw_hash,event,value);raw_events+=1
print(json.dumps(dict(phase='raw_semantic_digest_complete',events=raw_events)),flush=True)
cached_hash=hashlib.sha256();cached_events=0
def walk(value):
 global cached_events
 def emit(event,value=None):
  global cached_events
  add(cached_hash,event,value);cached_events+=1
 if isinstance(value,Mapping):
  emit('start_map')
  for key in value:emit('map_key',key);walk(value[key])
  emit('end_map')
 elif isinstance(value,Sequence) and not isinstance(value,(str,bytes)):
  emit('start_array')
  for v in value:walk(v)
  emit('end_array')
 elif value is None:emit('null')
 elif type(value) is bool:emit('boolean',value)
 elif type(value) in (int,float):emit('number',value)
 elif type(value) is str:emit('string',value)
 else:raise ValueError('unsupported cached type')
walk(run)
assert cached_events==raw_events and cached_hash.digest()==raw_hash.digest()
run.verify_source();assert sha(source)==expected
index=read(out/'components/index.json')
v=dict(passed=True,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_sha256=expected,
 source_bytes=source.stat().st_size,frames=rows,all6_layout_fields_equal_to_Rust=True,
 whole_raw_vs_cached_semantic_event_digest_equal=True,semantic_event_count=raw_events,semantic_sha256=raw_hash.hexdigest(),
 includes_every_geostatic_seed_committed_transition_removed_history_and_diagnostic_field=True,
 source_and_component_hashes_verified_before_completion=True,components=len(index['components']),
 maximum_component_accounted_bytes=max(d['payload'] for d in index['components']),
 maximum_component_encoded_bytes=max(d['bytes'] for d in index['components']),
 elapsed_seconds=time.monotonic()-start,complete_existing_mechanical_audit_and_exporter_not_yet_connected=True,
 this_checks_storage_not_constitutive_accuracy=True,p6_complete=False)
with (out/'verification.json').open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
print(json.dumps(v),flush=True)
