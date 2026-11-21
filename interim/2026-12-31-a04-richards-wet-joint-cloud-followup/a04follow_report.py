"""Fixed-order measures of newly solved fields; no acceptance extrapolation."""
import argparse
from collections import Counter
import itertools
from a04follow_common import *

def dot(a,b):
    x=np.asarray(a,dtype=float).ravel(order='C');y=np.asarray(b,dtype=float).ravel(order='C')
    if x.shape!=y.shape:raise ValueError('dot shape')
    return math.fsum(float(v)*float(w) for v,w in zip(x,y))

def coordinate(row):
    return tuple((k,row[k]) for k in ('phase','name','stage','time_s','progress') if k in row)

def soil_observations(folder,engine,common=True):
    mesh=read(folder/'mesh.json');layout=read(folder/(engine+'-layout.json'))
    data=[];raw=[]
    for item in rows(folder/(engine+'-fields.jsonl')):
        row=item['row'];f=item['fields']
        if common and (row['phase']!='transient' or row['time_s']*4 != round(row['time_s']*4)):continue
        hm=layout['full' if all(f['active_elements']) else 'cut']
        pos=np.asarray(f['hydraulic_positions']);w=np.asarray(f['hydraulic_weights']);gp=len(w)//len(f['cell_ids'])
        hh=np.repeat(np.asarray(f['cell_head_m']),gp)
        moments=lambda x:[dot(w,np.asarray(x)*pos[:,2]**k) for k in range(4)]
        v=dict(head_moments=moments(hh),chip_moments=moments(f['hydraulic_effective_pore_stress_kpa']),
               water_density_moments=moments(f['hydraulic_water_equivalent_mass_per_volume']),
               water_m3=f['water_m3'],removed_water_m3=f['removed_water_m3'],net_inflow_m3=f['net_inflow_m3'],mass_defect_m3=f['mass_defect_m3'])
        raw.append(dict(row=row,**v));data.append(v)
    if common and len(data)!=16:raise ValueError(('common physical times',len(data)))
    result={}
    for k in data[0]:
        if k.endswith('_moments'):
            for power in range(4):result[k+'_z'+str(power)]=[v[k][power] for v in data]
        else:result[k]=[v[k] for v in data]
    return result,raw

def rates(values):
    a,b,c=map(lambda v:np.asarray(v,dtype=float),values);d1=a-b;d2=b-c;n1=norm(d1);n2=norm(d2)
    return dict(coarse_medium=compare(a,b),medium_fine=compare(b,c),observed_log2_rate=math.log2(n1/n2) if n1 and n2 else None,
                difference_cosine=dot(d1,d2)/(n1*n2) if n1 and n2 else None,extrapolation=False)

def full_soil_vectors(folder,engine):
    selected=[]
    for item in rows(folder/(engine+'-fields.jsonl')):
        row=item['row']
        if row['phase']=='transient' and row['time_s']*4==round(row['time_s']*4):selected.append(item['fields'])
    assert len(selected)==16
    return {k:np.concatenate([np.asarray(f[k],float).ravel(order='C') for f in selected]) for k in ('cell_head_m','face_head_m','fluid_internal','water_body','u','reaction')}

def soil_report(root):
    table={};raw={};agreement={};full={};full_agreement={};whole={};native_table={};native_full={}
    for h,n,q in itertools.product((2,4,8),(4,8,16),(3,5,7)):
        folder=root/'soil'/f'h{h}-n{n}-q{q}'
        key=f'h{h}-n{n}-q{q}';v,r=soil_observations(folder,'independent');table[(h,n,q)]=v;raw[key]=r
        full[h,n,q]=full_soil_vectors(folder,'independent')
        all_rows=list(rows(folder/'independent-fields.jsonl'))
        whole[key]=dict(accepted_states=len(all_rows),maximum_mass_defect_m3=max(abs(v['fields']['mass_defect_m3']) for v in all_rows))
        if q==3:
            nf=full_soil_vectors(folder,'native');full_agreement[key]={k:compare(nf[k],full[h,n,q][k]) for k in nf}
            native_rows=list(rows(folder/'native-fields.jsonl'));whole[key]['native_maximum_mass_defect_m3']=max(abs(v['fields']['mass_defect_m3']) for v in native_rows)
            nv,nr=soil_observations(folder,'native');raw[key+'-native']=nr;native_table[h,n]=nv;native_full[h,n]=nf
            agreement[key]={k:compare(nv[k],v[k]) for k in v}
    h_rates={};dt_rates={};q_differences={};mixed={};full_dt={};full_q={}
    for n,q in itertools.product((4,8,16),(3,5,7)):
        h_rates[f'n{n}-q{q}']={k:rates([table[h,n,q][k] for h in (2,4,8)]) for k in table[2,n,q]}
    for h,q in itertools.product((2,4,8),(3,5,7)):
        full_dt[f'h{h}-q{q}']={k:rates([full[h,n,q][k] for n in (4,8,16)]) for k in full[h,4,q]}
        dt_rates[f'h{h}-q{q}']={k:rates([table[h,n,q][k] for n in (4,8,16)]) for k in table[h,4,q]}
    for h,n in itertools.product((2,4,8),(4,8,16)):
        full_q[f'h{h}-n{n}']={k:dict(q3_q5=compare(full[h,n,3][k],full[h,n,5][k]),q5_q7=compare(full[h,n,5][k],full[h,n,7][k])) for k in full[h,n,3]}
        q_differences[f'h{h}-n{n}']={k:dict(q3_q5=compare(table[h,n,3][k],table[h,n,5][k]),q5_q7=compare(table[h,n,5][k],table[h,n,7][k]),not_geometric_quadrature_rate=True) for k in table[h,n,3]}
    for h,n,q in itertools.product((2,4),(4,8),(3,5)):
        # Signed eight-corner finite difference. Norms are not additive percentages.
        values={k:np.zeros_like(np.asarray(v,float)) for k,v in table[h,n,q].items()}
        for ih,inn,iq in itertools.product((0,1),repeat=3):
            for k in values:values[k]+=(-1)**(ih+inn+iq)*np.asarray(table[h*2**ih,n*2**inn,q+2*iq][k])
        mixed[f'h{h}-n{n}-q{q}']={k:dict(absolute_l2=norm(v),signed_values=v) for k,v in values.items()}
    pairwise={}
    for h,n,q in itertools.product((2,4),(4,8),(3,5)):
        for axis,left,right in (('h_dt',(2*h,2*n,q),(h,n,q)),('h_q',(2*h,n,q+2),(h,n,q)),('dt_q',(h,2*n,q+2),(h,n,q))):
            if axis=='h_dt':edges=((2*h,n,q),(h,2*n,q))
            elif axis=='h_q':edges=((2*h,n,q),(h,n,q+2))
            else:edges=((h,2*n,q),(h,n,q+2))
            pairwise[f'{axis}-h{h}-n{n}-q{q}']={}
            for k in table[right]:
                v=np.asarray(table[left][k])-np.asarray(table[edges[0]][k])-np.asarray(table[edges[1]][k])+np.asarray(table[right][k])
                pairwise[f'{axis}-h{h}-n{n}-q{q}'][k]=dict(absolute_l2=norm(v),signed_values=v)
    native_h={f'n{n}':{k:rates([native_table[h,n][k] for h in (2,4,8)]) for k in native_table[2,n]} for n in (4,8,16)}
    native_dt={f'h{h}':{k:rates([native_full[h,n][k] for n in (4,8,16)]) for k in native_full[h,4]} for h in (2,4,8)}
    return dict(raw=raw,native_h_directional_rates=native_h,native_full_fields_same_mesh_dt=native_dt,pairwise_interactions=pairwise,whole_path=whole,full_fields_same_mesh_dt=full_dt,full_fields_same_mesh_quadrature=full_q,full_fields_native_reference=full_agreement,h_directional_rates=h_rates,dt_rates=dt_rates,quadrature=q_differences,three_factor_interactions=mixed,native_reference=agreement,
                moment_definition='integral z^k times H, chi*p, or water equivalent density; k=0..3; fixed fsum',
                full_spatial_field_norm=False,finest_is_exact=False)

MATERIAL=set('stress strain void total_stress pore_pressure head plastic_strain equivalent_plastic_strain plastic_work material_return'.split())
CONTACT=set('jump traction reference_jump initial_traction plastic_slip accumulated plastic_dissipation numerical_dissipation endpoint_work elastic_offset initial_energy contact_status contact_open opening gap normal_pressure shear_traction elastic_energy frozen_jump frozen_traction frozen_plastic frozen_initial frozen_reference frozen_rest frozen_accumulated frozen_pd frozen_nd frozen_work frozen_initial_energy frozen_status joint_water_pressure'.split())

def canonical_joint(f,layout):
    ek=layout['element_keys'];eo=sorted(range(len(ek)),key=lambda i:tuple(ek[i]));ck=layout['contact_keys']
    co=sorted(range(len(ck)),key=lambda i:json.dumps(ck[i],separators=(',',':')))
    out={}
    for k,v in f.items():
        a=np.asarray(v)
        if k in MATERIAL and len(a):
            gp=len(a)//len(ek);a=a[[e*gp+j for e in eo for j in range(gp)]]
        elif k in CONTACT:a=a[co]
        out[k]=a
    return out

def order_vectors(folder,engine,kind):
    layout=read(folder/(('' if kind=='soil' else 'a04follow-')+engine+'-layout.json'));vectors={};previous=None;loads=[];reactions=[];coordinates=[]
    for item in rows(folder/(engine+'-fields.jsonl')):
        row=item['row'];f=item['fields'];coordinates.append(coordinate(row))
        if kind=='soil':
            ek=layout['element_keys'];cells=f['cell_ids'];order=sorted(range(len(cells)),key=lambda i:tuple(ek[cells[i]]))
            f={**f,'cell_head_m':np.asarray(f['cell_head_m'])[order]}
        else:f=canonical_joint(f,layout)
        u=np.asarray(f['u']);du=u-(np.zeros_like(u) if previous is None else previous);previous=u.copy()
        loads.append(dot(f['applied_load'],du));reactions.append(dot(f['reaction'],du))
        for k in ('u','traction','cell_head_m','reaction','fluid_internal','water_body','endpoint_work','plastic_dissipation','numerical_dissipation'):
            if k in f:vectors.setdefault(k,[]).append(f[k])
        vectors.setdefault('fixed_nodal_load_work',[]).append(math.fsum(loads))
        vectors.setdefault('fixed_nodal_reaction_work',[]).append(math.fsum(reactions))
    return {k:np.concatenate([np.asarray(v,float).ravel(order="C") for v in values]) for k,values in vectors.items()},coordinates,layout

def normalized_layout(layout,kind):
    if kind=='soil':
        ek=layout['element_keys'];value={}
        for name in ('full','cut'):
            h=layout[name];indices=sorted(range(len(h['elements'])),key=lambda i:tuple(ek[h['elements'][i]]));points=[]
            for i in indices:
                points.append(sorted(zip(h['positions'][i],h['weights'][i]),key=lambda v:tuple(v[0])))
            value[name]=dict(elements=[ek[h['elements'][i]] for i in indices],points=points,face_keys=h['face_keys'],
                              owners=sorted((tuple(key),sorted((tuple(ek[e]),j) for e,j in owners)) for key,owners in zip(h['face_keys'],h['owners'])))
        return value
    ck=layout['contact_keys'];co=sorted(range(len(ck)),key=lambda i:json.dumps(ck[i],separators=(',',':')))
    return dict(elements=sorted(layout['element_keys']),contacts=[ck[i] for i in co],owner_elements=[layout['owner_elements'][i] for i in co] if 'owner_elements' in layout else None,positions=[layout['positions'][i] for i in co] if 'positions' in layout else None,normal=np.asarray(layout['normal'])[co].tolist(),weight=np.asarray(layout['weight'])[co].tolist(),c=np.asarray(layout['c'])[co].tolist())

def order_report(root,surface=None):
    report={};audits={}
    for kind,eng in itertools.product(('soil','dry','hydro'),('native','independent')):
        vectors={};layouts={};coords=None
        for e,f,q in itertools.product((False,True),repeat=3):
            key=f'e{int(e)}f{int(f)}q{int(q)}'
            folder=surface/key if surface is not None and kind=='hydro' and eng=='independent' and q else root/'order'/kind/key
            v,c,lay=order_vectors(folder,eng,kind)
            if coords is None:coords=c
            elif coords!=c:raise ValueError('ordering coordinates changed')
            vectors[(e,f,q)]=v;layouts[key]=normalized_layout(lay,kind)
        baseline=layouts['e0f0q0'];same={k:v==baseline for k,v in layouts.items()}
        if not all(same.values()):raise ValueError(('point/ownership mapping changed',kind,eng,same))
        effects={}
        for axis,name in enumerate(('E','F','Q')):
            pairs=[]
            for flags in itertools.product((False,True),repeat=3):
                if flags[axis]:continue
                other=list(flags);other[axis]=True
                pairs.append(dict(left=flags,right=other,fields={k:compare(vectors[flags][k],vectors[tuple(other)][k]) for k in vectors[flags]}))
            effects[name]=pairs
        report[kind+'-'+eng]=dict(states=len(coords),effects=effects,head_available=kind=='soil',joint_head=None if kind!='soil' else 'not_a_joint',roundoff_attribution_requires_identical_point_sets=True)
        audits[kind+'-'+eng]=dict(all_point_sets_and_owners_exact=True,permutations=same,normalized_layout=baseline)
    return report,audits

def exchange_report(root):
    vectors={};raw={};partitions={}
    for eng,p,r in itertools.product(('native','independent'),('N','I'),('N','I')):
        key=eng+'-p'+p+'r'+r;folder=root/'exchange'/key
        data=list(rows(folder/'fields.jsonl'));raw[key]=[v['row'] for v in data]
        vectors[key]={k:[v['fields'][k] for v in data] for k in data[0]['fields']}
        if eng=='native':
            parts=list(rows(folder/'raw/constitutive.jsonl'));partitions[key]=[v['partitions'] for v in parts]
        else:partitions[key]=[v['contact_partitions'] for v in read(folder/'result.json')['rows']]
    def pair(a,b):
        aa=vectors[a];bb=vectors[b]
        ca=np.asarray(aa['contact_status']);cb=np.asarray(bb['contact_status']);mask=ca!=cb
        def margin(v):
            tr=np.asarray(v['traction']);press=np.asarray(v['effective_contact_pressure']);limit=np.asarray(v['friction_limit'])
            # This fixture's horizontal seam has normal [0,0,1] up to side sign.
            shear=np.sqrt(np.sum(tr[...,:2]**2,axis=-1))
            return np.abs(limit-shear)/(1+limit+shear)
        ma=margin(aa);mb=margin(bb)
        pa=np.asarray(aa['plastic_slip']);pb=np.asarray(bb['plastic_slip'])
        da=np.diff(pa,axis=0,prepend=np.zeros_like(pa[:1]));db=np.diff(pb,axis=0,prepend=np.zeros_like(pb[:1]))
        zeros=np.all(da==0,axis=-1)&np.all(db==0,axis=-1)
        return dict(fields={k:compare(aa[k],bb[k]) for k in ('u','traction','plastic_slip','plastic_dissipation','numerical_dissipation','endpoint_work','reaction')},
                    contact_category_mismatches=int(np.count_nonzero(mask)),
                    category_observation=dict(both_return_margins_le_1e14=int(np.count_nonzero(mask&(ma<=1e-14)&(mb<=1e-14))),
                        both_saved_plastic_increments_exact_zero=int(np.count_nonzero(mask&zeros)),
                        maximum_scaled_return_margin=max(float(np.max(ma[mask],initial=0)),float(np.max(mb[mask],initial=0))),
                        maximum_plastic_increment_m=max(float(np.max(np.abs(da[mask]),initial=0)),float(np.max(np.abs(db[mask]),initial=0))),
                        classification_band_is_not_acceptance=True),
                    native=str(a),reference=str(b))
    matches={p+r:pair('native-p'+p+'r'+r,'independent-p'+p+'r'+r) for p,r in itertools.product(('N','I'),repeat=2)}
    effects={eng:{axis:[pair(eng+'-p'+p+'r'+r,eng+'-p'+pp+'r'+rr) for p,r,pp,rr in
                    (('N','N','I','N'),('N','I','I','I')) if axis=='P'] if axis=='P' else
                    [pair(eng+'-p'+p+'r'+r,eng+'-p'+pp+'r'+rr) for p,r,pp,rr in (('N','N','N','I'),('I','N','I','I'))] for axis in ('P','R')} for eng in ('native','independent')}
    return dict(original_fresh=pair('native-pNrN','independent-pIrI'),fully_exchanged=pair('native-pIrI','independent-pNrN'),same_policy=matches,effects=effects,rows=raw,
                partition_counts={k:dict(Counter(str(x) for row in vals for x in (row[-6:] if k.startswith('native') else row))) for k,vals in partitions.items()},
                old_saved_values_not_rewritten=dict(dissipation_relative_l2=0.0212002682,stick_slip_mismatches=133))

def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--surface',type=Path)
    a=p.parse_args();order,audits=order_report(a.run,a.surface)
    write(a.output,dict(passed=False,a04_complete=False,diagnostic_only=True,soil=soil_report(a.run),ordering=order,exchange=exchange_report(a.run),not_new_acceptance=True))
    write(a.output.with_name('order-point-set-audit.json'),audits)

if __name__=='__main__':main()
