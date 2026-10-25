import unittest
import numpy as np
from surface_projection_reference import upward_area

class ProjectionTests(unittest.TestCase):
    def test_curved_signed_projection_exact(self):
        r=np.array([0.,1.,0.,.5,.5,0.]);s=np.array([0.,0.,1.,0.,.5,.5]);p=np.column_stack((r+.5*r*r,(r-.5)*s,s))
        area,error=upward_area(p);self.assertLess(abs(area-7/192),error+1e-14)
        area,error=upward_area(p[[0,2,1,5,4,3]]);self.assertLess(abs(area-23/192),error+1e-14)
    def test_planar_shifted_projection(self):
        area,error=upward_area(np.array([[1e8,1e8,1e8],[1e8+3,1e8,1e8],[1e8,1e8+3,1e8]]));self.assertAlmostEqual(area,4.5);self.assertEqual(error,0.)

if __name__=='__main__':unittest.main()
