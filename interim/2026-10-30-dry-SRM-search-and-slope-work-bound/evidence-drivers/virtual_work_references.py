"""Bind the independent dense-FEM arrays to actual Rust CLI input and output."""
from pathlib import Path
import copy,hashlib,json,math,re,struct
import numpy as np

def read(q):return json.loads(Path(q).read_text('utf-8-sig'))
def dump(q,v):
    with Path(q).open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
def prepare(runtime,out,reference):
    terminal=read(reference/'terminal.json');assert terminal['passed'] and terminal['actual_cases']==6
    base=read(runtime/'examples/mc-virtual-work-mesh.json')
    prototype=read(runtime/'examples/mc-virtual-work-problem.json')
    rows=[]
    for index,row in enumerate(terminal['rows']):
        order=row['order'];label=row['case'];factor=row['factor']
        source=reference/f'full-fields-Tet{order}-{label}.json';fields=read(source)
        mesh=copy.deepcopy(base);mesh['nodes']=fields['nodes']
        if order==10:
            mesh['schema']='soil-fem-mesh/3';mesh['tetrahedra10']=fields['cells'];mesh['quadrature_order']=3
        else:
            mesh['tetrahedra']=fields['cells'];mesh.pop('tetrahedra10',None);mesh.pop('quadrature_order',None)
        mesh['node_sets'].update(all=list(range(len(mesh['nodes']))),
            left=[i for i,x in enumerate(mesh['nodes']) if x[0]==0.],
            bottom=[i for i,x in enumerate(mesh['nodes']) if x[2]==0.])
        p=copy.deepcopy(prototype);params=p['materials'][0]['parameters']
        params['cohesion_kpa']=(10. if label=='uniaxial' else 100.)/factor
        phi=math.atan(math.tan(math.radians(30.))/factor);params['friction_deg']=math.degrees(phi)
        ratio=(1.+math.sin(phi))/(1.-math.sin(phi))
        beta=0. if label=='uniaxial' else .02
        p['physical_nodal_force_kn']=fields['force_kn']
        p['virtual_displacement_m']=[[(ratio+.2)*x[0]+beta*x[0]*x[2],0.,-x[2]] for x in mesh['nodes']]
        prefix=out/f'independent-case-{index}'
        meshfile=prefix.with_suffix('.mesh.json');problem=prefix.with_suffix('.problem.json')
        dump(meshfile,mesh);dump(problem,p)
        rows.append(dict(index=index,mesh=str(meshfile),problem=str(problem),
            output=str(prefix.with_suffix('.result.json')),reference=str(source),
            reference_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),expected=row))
    dump(out/'independent-reference-binding.json',dict(passed=True,cases=rows,
        source_terminal_sha256=hashlib.sha256((reference/'terminal.json').read_bytes()).hexdigest(),
        no_reference_tolerance_or_data_changed=True))
    return rows

def compact_outside_strings(text):
    quoted=escaped=False;out=[]
    for c in text:
        if quoted:
            out.append(c)
            if escaped:escaped=False
            elif c=='\\':escaped=True
            elif c=='"':quoted=False
        elif c=='"':quoted=True;out.append(c)
        elif not c.isspace():out.append(c)
    assert not quoted and not escaped
    return ''.join(out).encode('utf-8')

def point_stream_hash(text):
    markers=list(re.finditer(r'"point_checks"\s*:\s*\[',text));assert len(markers)==1
    pos=markers[0].end();decoder=json.JSONDecoder();digest=hashlib.sha256();count=0
    while True:
        while text[pos].isspace():pos+=1
        if text[pos]==']':break
        assert text[pos]=='{'
        _,end=decoder.raw_decode(text,pos)
        raw=compact_outside_strings(text[pos:end])
        digest.update(struct.pack('<Q',len(raw)));digest.update(raw);count+=1;pos=end
        while text[pos].isspace():pos+=1
        if text[pos]==',':pos+=1
        else:assert text[pos]==']'
    return digest.hexdigest(),count

def verify(row):
    text=Path(row['output']).read_text('utf-8');result=json.loads(text)
    reference=read(row['reference']);expected=row['expected'];problem=read(row['problem'])
    assert hashlib.sha256(Path(row['reference']).read_bytes()).hexdigest()==row['reference_sha256']
    assert result['complete_point_audit'] and result['assessment']['status']=='evaluated'
    c=result['assessment']['certificate'];points=result['point_checks'];actual_hash,count=point_stream_hash(text)
    assert actual_hash==c['point_checks_sha256'] and count==c['integration_points']==len(points)==len(reference['points'])
    assert c['excludes_admissible_discrete_equilibrium']==expected['excluded']
    assert c['all_constrained_components_exactly_zero'] and not c['continuum_collapse_or_stability_certified']
    assert abs(c['external_virtual_work_kn_m']-expected['external_work'])<=2e-12*max(1.,abs(expected['external_work']))
    assert abs(c['maximum_admissible_internal_virtual_work_kn_m']-expected['integrated_capacity'])<=2e-12*max(1.,abs(expected['integrated_capacity']))
    phi=math.radians(problem['materials'][0]['parameters']['friction_deg'])
    a=problem['materials'][0]['parameters']['cohesion_kpa']/math.tan(phi);ratio=(1.+math.sin(phi))/(1.-math.sin(phi))
    maximum=0.
    def near(actual,expected):
        nonlocal maximum
        error=abs(actual-expected)/max(1.,abs(expected));maximum=max(maximum,error)
        assert error<=2e-12,(actual,expected,error)
    key='virtual_engineering_strain_factor'+str(expected['factor'])
    per_element=1 if expected['order']==4 else 27
    for index,(actual,source) in enumerate(zip(points,reference['points'])):
        assert actual['element_index']==source['element']==index//per_element and actual['local_index']==index%per_element
        assert actual['material_id']==1 and actual['accepted_finite_support']
        near(actual['weight_m3'],source['weight_m3'])
        for av,ev in zip(actual['position_m'],source['position_m']):near(av,ev)
        eng=np.asarray(source[key]);D=np.array([[eng[0],eng[3]/2.,eng[5]/2.],[eng[3]/2.,eng[1],eng[4]/2.],[eng[5]/2.,eng[4]/2.,eng[2]]])
        scale=float(np.max(np.abs(D)));normal=D/scale;eigen=np.linalg.eigvalsh(normal)[::-1]
        near(actual['strain_tensor_scale'],scale)
        voigt=[normal[0,0],normal[1,1],normal[2,2],normal[0,1],normal[1,2],normal[0,2]]
        for av,ev in zip(actual['normalized_strain_tensor_tension_positive'],voigt):near(av,float(ev))
        for av,ev in zip(actual['normalized_principal_strains_descending'],eigen):near(av,float(ev))
        margin=float(np.trace(normal)+(ratio-1.)*np.minimum(eigen,0.).sum())/((1.+ratio)*float(np.linalg.norm(normal)))
        near(actual['relative_dual_margin'],margin)
        near(actual['support_density_kpa'],a*float(np.trace(D)))
    return dict(passed=True,index=row['index'],order=expected['order'],case=expected['case'],factor=expected['factor'],
        all_GP_checked=count,point_hash_framing_verified=True,maximum_scaled_difference=maximum,
        expected_excluded=expected['excluded'],actual_excluded=c['excludes_admissible_discrete_equilibrium'])
