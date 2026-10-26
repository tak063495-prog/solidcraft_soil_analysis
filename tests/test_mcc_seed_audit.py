"""Original-coefficient checks for affine seeds and global audit entry points."""
import copy,json,math,os
from pathlib import Path
import subprocess,sys,tempfile,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from mcc_seed_audit import verify_affine_seed


class AffineSeedAudit(unittest.TestCase):
    def test_reconstruction_rejects_modified_state_coefficients_and_nonzero_plastic_history(self):
        mesh=dict(nodes=[[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]],tetrahedra=[[0,1,2,3]],
            region_ids=[1],material_ids=[1],regions=[dict(id=1,body='clay')])
        field=dict(body='clay',reference_position_m=[0.,0.,0.],
            effective_stress_at_reference_compression_positive_kpa=[100.,100.,100.,0.,0.,0.],
            effective_stress_gradient_kpa_per_m=[[0.,0.,-12.],[0.,0.,-12.],[0.,0.,-20.],[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]],
            pc_over_mean_effective=1.8,void_ratio=dict(mode='from_mcc_reference'))
        m=dict(critical_state_m=1.4,**{'lambda':.12,'kappa':.025},poisson=.2,
            reference_pressure_kpa=1.,specific_volume_ncl_at_reference=3.)
        ground=dict(materials=[dict(id=1,mcc=m)],bodies=[dict(body='clay',material_id=1)],
            initialization=dict(geostatic=dict(seed=dict(kind='mcc_affine_fields',fields=[field]))))
        p=(97.+97.+95.)/3;pc=p*1.8
        states=[dict(effective_stress_compression_positive_kpa=[97.,97.,95.,0.,0.,0.],
            preconsolidation_pressure_kpa=pc,void_ratio=3.-.12*math.log(pc)+.025*math.log(1.8)-1.,
            plastic_strain_compression_positive=[0.]*6)]
        verify_affine_seed(mesh,ground,states)
        for failure in ('stress','pc','void','plastic','extra','inventory','coefficients'):
            trial=copy.deepcopy(states);definition=copy.deepcopy(ground)
            if failure=='stress':trial[0]['effective_stress_compression_positive_kpa'][0]+=.001
            elif failure=='pc':trial[0]['preconsolidation_pressure_kpa']+=.001
            elif failure=='void':trial[0]['void_ratio']+=1e-7
            elif failure=='plastic':trial[0]['plastic_strain_compression_positive'][0]=1e-15
            elif failure=='extra':trial[0]['unknown']=0.
            elif failure=='inventory':trial=[]
            else:definition['initialization']['geostatic']['seed']['fields'][0]['pc_over_mean_effective']+=.01
            with self.subTest(failure=failure),self.assertRaisesRegex(ValueError,'affine MCC'):
                verify_affine_seed(mesh,definition,trial)

    def test_affine_seed_global_audits_cover_both_mixed_forms_and_displacement_form(self):
        from prepare_mcc_structures import inputs
        import fem_reference as ref
        from verify_projected_stage_package import audit
        from verify_geostatic import verify_initial
        from verify_high_order import prepare
        root=Path(__file__).resolve().parents[1]
        binary=root/'target/release'/('soil-fem-3d.exe' if os.name=='nt' else 'soil-fem-3d')
        if not binary.is_file():self.skipTest('build the release CLI before global audit tests')
        for schema in ('soil-fem-projected-mcc-stages/1','soil-fem-mixed-mcc-stages/1','soil-fem-stages/1'):
            mesh,_,case=inputs(root);case['schema']=schema
            mesh['face_sets']={'top':[list(face) for face,owners in ref.topology(mesh).items()
                if len(owners)==1 and all(mesh['nodes'][n][2]==1. for n in face)]}
            ground=json.loads((root/'examples/ground-geostatic-mcc-column.json').read_text())
            ground['materials'][1]['mcc']['specific_volume_ncl_at_reference']=3.1
            ground['initialization']['geostatic']['seed']=dict(kind='mcc_affine_fields',fields=[
                dict(body=b['body'],reference_position_m=[0.,0.,0.],
                    effective_stress_at_reference_compression_positive_kpa=[100.,100.,100.,0.,0.,0.],
                    effective_stress_gradient_kpa_per_m=[[0.,0.,0.]]*6,
                    pc_over_mean_effective=1.4,void_ratio=dict(mode='from_mcc_reference')) for b in ground['bodies']])
            case['stages']=case['stages'][1:3];case['verification_reject_once_stages']=['install']
            with tempfile.TemporaryDirectory() as tmp:
                paths=[Path(tmp)/name for name in ('mesh.json','ground.json','stages.json')]
                for path,value in zip(paths,(mesh,ground,case),strict=True):path.write_text(json.dumps(value))
                result=Path(tmp)/'result.json'
                subprocess.run([str(binary),'solve-stages',*map(str,paths),str(result)],check=True,capture_output=True,text=True)
                run=json.loads(result.read_text());self.assertTrue(run['complete'])
                def check(definition,data):
                    if schema=='soil-fem-stages/1':
                        q,offset=prepare(mesh,data['integration_points']);verify_initial(mesh,definition,data,q,offset)
                    else:self.assertEqual(len(audit(mesh,definition,case,data)),3)
                check(ground,run)
                for failure in ('history','recipe'):
                    changed=copy.deepcopy(run);definition=copy.deepcopy(ground)
                    if failure=='history':changed['geostatic']['seed']['material_states'][0]['plastic_strain_compression_positive'][0]=1e-15
                    else:
                        definition['initialization']['geostatic']['seed']['fields'][0]['pc_over_mean_effective']+=.01
                        changed['ground']=ref.normalize_ground(copy.deepcopy(definition))
                    with self.subTest(schema=schema,failure=failure),self.assertRaisesRegex(ValueError,'affine MCC'):
                        check(definition,changed)


if __name__=='__main__':unittest.main()
