"""Byte-identical assembly of independently frozen process groups; no solving."""
import argparse
import os
import shutil
from a04follow_common import *

def copy(source,destination):
    if source.is_dir():shutil.copytree(source,destination,copy_function=os.link,ignore=shutil.ignore_patterns('__pycache__'))
    else:destination.parent.mkdir(parents=True,exist_ok=True);os.link(source,destination)

def main():
    p=argparse.ArgumentParser();p.add_argument('--soil',type=Path,required=True);p.add_argument('--order',type=Path,required=True);p.add_argument('--exchange',type=Path,required=True);p.add_argument('--surface',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--soil-source-root',type=Path,required=True);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False);metadata=a.output/'component-metadata';metadata.mkdir();bindings={}
    for kind,root in (('soil',a.soil),('order',a.order),('exchange',a.exchange)):
        copy(root/kind,a.output/kind);meta=metadata/kind;meta.mkdir()
        for src in root.iterdir():
            if src.is_file() and src.suffix in ('.json','.log'):copy(src,meta/src.name)
        copy(root/'source-snapshot',meta/'source-snapshot')
        binding=read(root/'source-binding.json')
        for path,digest in binding.items():
            if sha(root/'source-snapshot'/path)!=digest:raise ValueError(('captured source changed',kind,path))
            if kind=='soil' and sha(a.soil_source_root/path)!=digest:raise ValueError(('soil execution source changed',path))
        bindings[kind]=dict(immutable_snapshot_entries=len(binding),source_binding_sha256=sha(root/'source-binding.json'))
    ledger=[]
    for folder in sorted((a.output/'soil').iterdir()):
        for eng in ('native','independent'):
            fields=folder/(eng+'-fields.jsonl')
            if not fields.exists():continue
            result=read(folder/('native.json' if eng=='native' else 'independent-result.json'));assert result['complete']
            ledger.append(dict(kind='soil',path=str(folder.relative_to(a.output)),engine=eng,states=sum(1 for _ in rows(fields))))
    for kind in ('soil','dry','hydro'):
        for folder in sorted((a.output/'order'/kind).iterdir()):
            if not folder.is_dir():continue
            for eng in ('native','independent'):
                native_result='native.json' if kind=='soil' else 'native/run.json'
                result=read(folder/(native_result if eng=='native' else 'independent-result.json'))
                if kind=='hydro' and eng=='independent':assert result['accepted']==16 and len(result['rows'])==16
                else:assert result['complete'] is True
                ledger.append(dict(kind='order',path=str(folder.relative_to(a.output)),engine=eng,states=sum(1 for _ in rows(folder/(eng+'-fields.jsonl')))))
    for folder in sorted((a.output/'exchange').iterdir()):
        if not folder.is_dir():continue
        result=read(folder/('raw/run.json' if folder.name.startswith('native') else 'result.json'));assert result['complete']
        ledger.append(dict(kind='exchange',path=str(folder.relative_to(a.output)),engine=folder.name.split('-')[0],states=sum(1 for _ in rows(folder/'fields.jsonl'))))
    main_count=len(ledger);main_states=sum(v['states'] for v in ledger)
    for folder in sorted(a.surface.glob('e*f*q1')):
        result=read(folder/'independent-result.json');assert result['accepted']==16 and len(result['rows'])==16
        ledger.append(dict(kind='full_surface',path=str(folder),engine='independent',states=sum(1 for _ in rows(folder/'independent-fields.jsonl'))))
    assert main_count==92 and len(ledger)==96,(main_count,len(ledger))
    assert main_states==6048 and sum(v['states'] for v in ledger)==6112
    write(a.output/'solve-ledger.json',dict(complete=True,main_paths=main_count,main_states=main_states,fresh_formal_paths=len(ledger),fresh_formal_accepted_states=sum(v['states'] for v in ledger),entries=ledger,
       solver_sources_fixed_per_process_group=bindings,assembly_is_not_additional_solves=True))

if __name__=='__main__':main()
