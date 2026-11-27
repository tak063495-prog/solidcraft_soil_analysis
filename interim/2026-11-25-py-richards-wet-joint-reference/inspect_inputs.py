"""Read cited immutable inputs/documents; import no production solver."""
import argparse
import hashlib
import json
from pathlib import Path
from reference import write,PARAMETERS
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def inventory():
    sources=json.loads((HERE/'sources.json').read_text(encoding='utf-8'))
    for name,digest in sources['files'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
            raise ValueError('cited source differs: '+name)
    folder=ROOT/'interim/2026-11-22-wet-joint-transient-water/inputs'
    hydraulics=json.loads((folder/'Richards-VG-4.json').read_text(encoding='utf-8'))
    mechanics=json.loads((folder/'mechanics-water-8.json').read_text(encoding='utf-8'))
    return dict(base_commit=sources['base_commit'],source_hashes_verified=True,
      C04_water_materials=hydraulics['materials'],C04_contact=mechanics['joint_system'],
      C04_initial_head=hydraulics['initial_head'],reference_inputs=PARAMETERS,
      identified_comparison_quantities=['both soil-side total head traces','interface total head',
       'instantaneous oriented flow on both faces','actual hydraulic storage change','cumulative exterior/reservoir inflow',
       'physical seconds and dt','hydraulic opening and mechanical opening separately','fixed generation and IDs',
       'physical pressure vs effective normal pressure','Coulomb shear/plastic slip/friction dissipation'],
      source_policy_meanings=dict(CONTINUOUS='both sides SOIL; share only hydraulic head',
         SOIL='explicit soil_side a/b; opposite impermeable structure excluded from hydraulics',
         PRESCRIBED='fixed total head on adjacent active SOIL faces; null level follows exterior water level'),
      unconfirmed=['Richards-to-JOINT production trace mapping','SOIL transient closed/open connection',
        'contact-triggered PRESCRIBED reservoir switching','two-way production controller/schema/CLI',
        'generation-change head/water/history transfer and initial contact force',
        'A04-specific fully coupled material/geometry/boundary dataset','coupled numerical acceptance tolerance',
        'real fracture aperture-dependent hydraulic calibration','general 3D bulk mechanical effective-stress solution'],
      tolerance=dict(coupled_acceptance='未定義',accepted=None,
        existing_C04_hydraulics=dict(relative_flow=1e-8,absolute_flow_m3_s=1e-16,head_m=1e-8),
        source='interim/2026-11-22-wet-joint-transient-water/REPRODUCE.ja.md: 出力と確認基準'),
      A04='未達',C04='未達')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();write(a.output,inventory())
