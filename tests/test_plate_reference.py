"""Independent polynomial plate physics and raw sample export contracts."""
import os
for _name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[_name]='1'
import sys
from pathlib import Path
import tempfile
import unittest
import numpy as np
import meshio
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from plate_reference import Element
import plate_io


class PlateReference(unittest.TestCase):
    def specimen(self):
        points=np.array([[0.,0.,0.],[2.,0.,0.],[.4,1.3,0.]])
        section=dict(young_kpa=2e8,poisson=.25,thickness_m=.02,drilling_factor=1e-4)
        element=Element(points,section);d=np.zeros(18)
        for i,(x,y,_) in enumerate(points):
            d[6*i:6*i+6]=[.002*x+.001*y,-.003*y+.004*x,
                -(.005*x*x+.007*y*y)/2-.002*x*y,-.007*y-.002*x,.005*x+.002*y,.0015]
        return points,section,element,d

    def test_polynomial_constraints_reproduce_constant_strain_curvature_work_and_rigid_modes(self):
        points,section,element,d=self.specimen();report=element.recover(d)
        strain=np.array([.002,-.003,.005]);curvature=np.array([.005,.007,.004])
        for sample in report['samples']:
            np.testing.assert_allclose(sample['membrane_strain'],strain,atol=1e-14,rtol=0.)
            np.testing.assert_allclose(sample['curvature_per_m'],curvature,atol=1e-14,rtol=0.)
            self.assertLess(abs(sample['drilling_mismatch_rad']),1e-15)
        expected=element.area/2*(strain@element.membrane@strain+curvature@element.bending@curvature)
        self.assertAlmostEqual(float(d@element.stiffness@d)/2,expected,places=10)
        for axis in range(6):
            rigid=np.zeros(18)
            for i,p in enumerate(points):
                if axis<3:rigid[6*i+axis]=1.
                else:
                    w=np.eye(3)[axis-3];rigid[6*i:6*i+3]=np.cross(w,p);rigid[6*i+3:6*i+6]=w
            self.assertLess(np.linalg.norm(element.stiffness@rigid),1e-12*np.linalg.norm(element.stiffness)*np.linalg.norm(rigid))
        values=np.linalg.eigvalsh(element.stiffness)
        self.assertEqual(np.count_nonzero(values>1e-10*values.max()),12)

    def test_plate_export_preserves_local_raw_samples_and_inactive_frames(self):
        points,_,element,d=self.specimen();plate=element.recover(d);plate['id']='wall'
        cfg=dict(nodes=[dict(id=str(i)) for i in range(3)],plates=[dict(id='wall',nodes=['0','1','2'])])
        empty=dict(beam_results=dict(plates=[]),state=dict(beams=dict(active_plates=[False])))
        active=dict(beam_results=dict(plates=[plate],node_positions_m=points),
            state=dict(beams=dict(active_plates=[True],node_displacement_m=d.reshape(3,6)[:,:3],node_rotation_rad=d.reshape(3,6)[:,3:])))
        run=dict(case=dict(beam_system=cfg),frames=[dict(name='before',**empty),dict(name='loaded',**active),dict(name='removed',**empty)])
        with tempfile.TemporaryDirectory() as tmp:
            output=Path(tmp);plate_io.export(run,output)
            raw=meshio.read(output/'plate-samples-001.vtu')
            np.testing.assert_allclose(raw.point_data['curvature_per_m'],np.tile([.005,.007,.004],(3,1)),rtol=0.,atol=1e-14)
            self.assertEqual(len(raw.points),3)
            self.assertTrue((output/'plates.pvd').exists())
            self.assertEqual(len((output/'plate-samples.csv').read_text(encoding='utf-8-sig').splitlines()),4)


if __name__=='__main__':unittest.main()
