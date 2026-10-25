"""Free/bonded tendon segments and physical soil transfer, without smoothing."""
import csv
import json
import xml.etree.ElementTree as ET
import meshio
import numpy as np
from beam_io import empty_frame


def export(run,output):
    root=ET.Element('VTKFile',type='Collection',version='0.1',byte_order='LittleEndian');collection=ET.SubElement(root,'Collection');rows=[]
    definitions=run['case']['anchor_system']['anchors'];ids={d['id']:i for i,d in enumerate(definitions)}
    for fi,frame in enumerate(run['frames']):
        points=[];displacement=[];transfer=[];lines=[];anchor_ids=[];segments=[];force=[];length=[];taut=[]
        for a in frame['anchor_results']['anchors']:
            first=len(points);points+=a['point_positions_m'];displacement+=a['point_displacement_m'];transfer+=a['point_soil_force_kn']
            state=frame['state']['anchors']['anchors'][ids[a['id']]]
            for si,(l,n,t) in enumerate(zip(a['segment_lengths_m'],a['segment_forces_kn'],a['segment_taut'])):
                lines.append([first+si,first+si+1]);anchor_ids.append(ids[a['id']]);segments.append(si);length.append(l);force.append(n);taut.append(t)
                rows.append({'stage_index':fi,'stage':frame['name'],'anchor_id':a['id'],'segment_index':si,'free_span':si==0,
                    'length_m':l,'tension_kn':n,'taut':t,'extension_from_cad_m':a['segment_extension_m'][si],
                    'reference_extension_m':state['reference_extensions_m'][si],'force_control_active':a['force_control_active']})
        name=f'anchors-{fi:03}.vtu'
        if not lines:empty_frame(output/name)
        else:
            pd={'total_displacement_from_cad_m':np.asarray(displacement),'force_from_anchor_to_soil_kn':np.asarray(transfer)}
            cd={'anchor_index':[np.asarray(anchor_ids)],'segment_index':[np.asarray(segments)],'free_span':[np.asarray(segments)==0],
                'length_m':[np.asarray(length)],'tension_kn':[np.asarray(force)],'taut':[np.asarray(taut,dtype=int)]}
            cd['free_span']=[cd['free_span'][0].astype(int)]
            meshio.write_points_cells(output/name,points,[('line',lines)],point_data=pd,cell_data=cd)
            check=meshio.read(output/name);assert np.array_equal(check.cells_dict['line'],lines)
            for k,v in pd.items():assert np.allclose(check.point_data[k],v)
            for k,v in cd.items():assert np.allclose(check.cell_data[k][0],v[0])
        ET.SubElement(collection,'DataSet',timestep=str(fi),group=frame['name'],part='0',file=name)
    ET.ElementTree(root).write(output/'anchors.pvd',encoding='utf-8',xml_declaration=True)
    with (output/'anchor-segments.csv').open('x',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else ['stage_index','anchor_id']);writer.writeheader();writer.writerows(rows)
    (output/'anchor-fields.json').write_text(json.dumps({'schema':'soil-fem-anchor-export/1','anchors':ids,
        'model':'Straight tension-only tendon, unbonded first segment, perfectly bonded subsequent segments; no pullout/slip capacity.',
        'kinematics':'CAD chord geometry and physical total displacement, independent of soil display RESET.',
        'point_force':'Force exerted by the anchor on its embedded host; self-equilibrated. No additional equivalent prestress is applied.'},indent=2))
