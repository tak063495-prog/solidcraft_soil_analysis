"""Standalone Darcy CLI, independent reference, export and negative outcomes."""
import argparse
import copy
import json
from pathlib import Path
from verify_checkpoint import invoke
from darcy_reference import verify
from export_darcy import export
import fem_reference as fem


def run(solver,output):
    output.mkdir(parents=True,exist_ok=False);root=Path(__file__).resolve().parents[1]
    mesh_path=root/'examples/darcy-residual-mesh.json';flow=root/'examples/darcy-residual-flow.json';target=output/'solution.json'
    invoke([str(solver),'solve-darcy',str(mesh_path),str(flow),str(target)],output/'solve.log')
    mesh=json.loads(mesh_path.read_text());case=json.loads(flow.read_text());result=json.loads(target.read_text());assert result['refinement_solves']
    comparison=verify(mesh,case,result);export(mesh,result,output/'results')
    unit_mesh=root/'examples/seepage-column-tet10-mesh.json';unit=json.loads(unit_mesh.read_text())
    exchange=json.loads((root/'examples/flow-seepage-column.json').read_text())
    exchange['heads']=[{'face_set':'bottom','head_m':2.,'gradient':[.3,0.,0.]}]
    exchange['fluxes']=[{'kind':'outward','face_set':'xmin','flux_m_s':3e-7},{'kind':'outward','face_set':'xmax','flux_m_s':-3e-7}]
    for name in ('circulation','balanced-sources'):
        if name=='balanced-sources':
            exchange['heads'][0]['gradient']=[0.,0.,0.];exchange['fluxes']=[]
            for material,rate in zip(exchange['materials'],(1e-6,-1e-6)):material['volume_source_per_s']=rate
        input_path=output/(name+'-case.json');target=output/(name+'-result.json');input_path.write_text(json.dumps(exchange,indent=2))
        invoke([str(solver),'solve-darcy',str(unit_mesh),str(input_path),str(target)],output/(name+'.log'))
        verify(unit,exchange,json.loads(target.read_text()))
    shifted=copy.deepcopy(exchange)
    for material in shifted['materials']:material['volume_source_per_s']=0.
    shifted['heads'][0]={'face_set':'bottom','head_m':1e12,'gradient':[.3,0.,0.]}
    shifted['fluxes']=[{'kind':'outward','face_set':'xmin','flux_m_s':3e-7},{'kind':'outward','face_set':'xmax','flux_m_s':-3e-7}]
    path=output/'large-datum-case.json';target=output/'large-datum-result.json';path.write_text(json.dumps(shifted,indent=2))
    invoke([str(solver),'solve-darcy',str(unit_mesh),str(path),str(target)],output/'large-datum.log');verify(unit,shifted,json.loads(target.read_text()))
    two=copy.deepcopy(unit);n=len(unit['nodes']);ne=len(unit['tetrahedra'])
    two['nodes'] += [[x+3,y,z] for x,y,z in unit['nodes']]
    for key in ('tetrahedra','tetrahedra10'):two[key] += [[v+n for v in row] for row in unit[key]]
    two['material_ids']*=2;two['region_ids'] += [v+2 for v in unit['region_ids']]
    regions=copy.deepcopy(unit['regions'])
    for r in regions:r.update(id=r['id']+2,body=r['body']+'Second',bounds_m=[[b[0]+3,b[1],b[2]] for b in r['bounds_m']])
    two['regions']+=regions;two['node_sets']={}
    faces=[list(f) for f,o in fem.topology(unit).items() if len(o)==1]
    two['face_sets']={'first':faces,'second':[[v+n for v in row] for row in faces]}
    two_path=output/'two-domains-mesh.json';two_path.write_text(json.dumps(two,indent=2))
    different=json.loads((root/'examples/flow-seepage-column.json').read_text());different['heads']=[{'face_set':'first','head_m':2.},{'face_set':'second','head_m':3.}]
    path=output/'two-domains-case.json';target=output/'two-domains-result.json';path.write_text(json.dumps(different,indent=2))
    invoke([str(solver),'solve-darcy',str(two_path),str(path),str(target)],output/'two-domains.log');verify(two,different,json.loads(target.read_text()))
    owners=fem.topology(two)
    for group,prefix in ((False,'first'),(True,'second')):
        for level,suffix in ((0.,'bottom'),(1.,'top')):
            two['face_sets'][prefix+'_'+suffix]=[list(f) for f,o in owners.items() if len(o)==1 and all((i>=n)==group and abs(two['nodes'][i][2]-level)<1e-12 for i in f)]
    two['material_ids'][ne:]=[i+2 for i in two['material_ids'][ne:]]
    for r in two['regions'][2:]:r['material_id']+=2
    contrast=copy.deepcopy(different);extra=copy.deepcopy(contrast['materials'])
    for m in extra:m['id']+=2;m['hydraulic_conductivity_m_s']=[[1e-20 if i==j else 0. for j in range(3)] for i in range(3)]
    contrast['materials']+=extra;contrast['heads']=[{'face_set':f'{prefix}_{suffix}','head_m':head} for prefix,suffix,head in [('first','bottom',2.),('first','top',1.),('second','bottom',5.),('second','top',4.)]]
    contrast_mesh=output/'contrast-mesh.json';contrast_mesh.write_text(json.dumps(two,indent=2));path=output/'contrast-case.json';target=output/'contrast-result.json';path.write_text(json.dumps(contrast,indent=2))
    invoke([str(solver),'solve-darcy',str(contrast_mesh),str(path),str(target)],output/'contrast.log');result_contrast=json.loads(target.read_text());verify(two,contrast,result_contrast)
    for point in result_contrast['points']:
        expected=1e-6 if point['element_index']<ne else 1e-20
        assert abs(point['darcy_flux_m_s'][2]-expected)/expected<1e-9
    dry=json.loads(json.dumps(case))
    for bc in dry['heads']:bc['head_m']=0.
    path=output/'unsaturated.json';path.write_text(json.dumps(dry));candidate=output/'incomplete.json'
    invoke([str(solver),'solve-darcy',str(mesh_path),str(path),str(candidate)],output/'incomplete.log',1)
    assert not json.loads(candidate.read_text())['complete']
    bad=json.loads(json.dumps(case));bad['heads'][0]['face_set']='missing'
    path=output/'bad-boundary.json';path.write_text(json.dumps(bad));target=output/'rejected.json'
    invoke([str(solver),'solve-darcy',str(mesh_path),str(path),str(target)],output/'bad-boundary.log',1);assert not target.exists()
    with (output/'verification.json').open('x') as stream:json.dump({'schema':'soil-fem-seepage-cli/1','passed':True,
        'original_row_refinement_exercised':True,'absolute_flow_normalization_exercised':True,'independent_domains_and_datum_exercised':True,
        'unsaturated_candidate_marked_incomplete':True,'bad_boundary_rejected':True,**comparison},stream,indent=2)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args();run(a.solver.resolve(),a.output_dir.resolve())
