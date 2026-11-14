"""Read actual input geometry; compare each shared ID/soil binding to PR25 algebra."""
import argparse,importlib.util,itertools,json
from pathlib import Path
import numpy as np
from b03_cases import HERE,read,write
from study import sha
spec=importlib.util.spec_from_file_location('pr25_geometry',HERE/'reference-pr25/geometry_reference.py')
geo=importlib.util.module_from_spec(spec);spec.loader.exec_module(geo)


def audit(folder,out):
    rows=[]
    for case in sorted(folder.iterdir()):
        if not (case/'stages.json').exists():continue
        m,c=read(case/'mesh.json'),read(case/'stages.json');cfg=c['beam_system'];node={n['id']:n for n in cfg['nodes']}
        sections={s['id']:s for s in cfg['plate_sections']};faces=[]
        for p in cfg['plates']:
            xyz=np.array([node[k]['position_m'] for k in p['nodes']]);n=np.cross(xyz[1]-xyz[0],xyz[2]-xyz[0]);n/=np.linalg.norm(n)
            e=p.get('attachment_offset_m',0.);t=sections[p['section_id']]['thickness_m']
            faces.append(dict(id=p['id'],nodes=p['nodes'],n=n,e=e,t=t,xyz=xyz))
        pairs=[]
        for a,b in itertools.combinations(faces,2):
            for i,ida in enumerate(a['nodes']):
                for j,idb in enumerate(b['nodes']):
                    na,nb=node[ida],node[idb]
                    shared=ida==idb
                    same_soil=na.get('soil_node_set') is not None and nb.get('soil_node_set') is not None and m['node_sets'][na['soil_node_set']]==m['node_sets'][nb['soil_node_set']]
                    if not (shared or same_soil):continue
                    arm_difference=b['e']*b['n']-a['e']*a['n']
                    # Common attachment position means opposite neutral displacement.
                    neutral1=a['xyz'][i]-a['e']*a['n'];neutral2=b['xyz'][j]-b['e']*b['n']
                    from_neutral=geo.mismatch(neutral1,a['n'],a['e'],neutral2,b['n'],b['e'])
                    common_neutral=geo.mismatch([0,0,0],a['n'],a['e'],[0,0,0],b['n'],b['e'])
                    theta=np.arccos(np.clip(a['n']@b['n'],-1,1))
                    stable=np.sqrt((b['e']-a['e'])**2+4*a['e']*b['e']*np.sin(theta/2)**2) if a['e']*b['e']>=0 else np.sqrt(a['e']**2+b['e']**2-2*a['e']*b['e']*np.cos(theta))
                    pairs.append(dict(faces=[a['id'],b['id']],node_ids=[ida,idb],shared_structural_id=shared,same_soil_attachment=same_soil,
                                      normals=[a['n'].tolist(),b['n'].tolist()],e_m=[a['e'],b['e']],theta_deg=float(np.rad2deg(theta)),
                                      arm_difference_m=arm_difference.tolist(),arm_difference_norm_m=float(np.linalg.norm(arm_difference)),
                                      production_shared_id_bound_m=64*np.finfo(float).eps*b['t'],
                                      neutral_difference_m=(neutral2-neutral1).tolist(),
                                      common_neutral_PR25=common_neutral,actual_neutral_PR25=from_neutral,
                                      closed_length_m=float(stable),closed_length_error_m=abs(float(stable)-np.linalg.norm(arm_difference))))
        rows.append(dict(case=case.name,pairs=pairs,accepted=None,scope='geometry only; actual CLI result is separate'))
    series=geo.series()
    result=dict(passed=False,b03_complete=False,accepted=None,reference_commit='11c3d1bc5b913b311526842a3ab5c330d5ff087c',
                reference_sha256=sha(HERE/'reference-pr25/geometry_reference.py'),actual_cases=rows,
                PR25_fresh_geometry_series_count=series['case_count'],
                PR25_series_maximum_vector_error_m=max(r['closed_vector_error_m'] for r in series['cases']),
                PR25_series_maximum_norm_error_m=max(r['closed_norm_error_m'] for r in series['cases']),
                actual_maximum_length_error_m=max((p['closed_length_error_m'] for r in rows for p in r['pairs']),default=0.))
    write(out,result);write(out.with_name('pr25-geometry-series.json'),series)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('folder',type=Path);p.add_argument('output',type=Path);a=p.parse_args();audit(a.folder,a.output)
