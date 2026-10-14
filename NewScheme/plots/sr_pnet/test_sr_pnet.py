import unittest
import numpy as np
from check_sr_pnet import calculate,poisson_errors


class Compatibility(unittest.TestCase):
    def test_known_chi2(self):
        z=np.zeros((5,1))
        parts=dict(DY=np.array([[10.],[20.],[1.],[4.],[2.]]),TTplusWW=z,VZ=z,
                   Data=np.full((5,1),10.))
        mc,d,mask,rows=calculate(parts)
        self.assertTrue(mask[0]); self.assertEqual(rows[0]['chi2'],0)
        self.assertAlmostEqual(rows[1]['chi2'],100./14)
        self.assertEqual(rows[0]['ndof'],rows[1]['ndof'])

    def test_shared_mask(self):
        z=np.zeros((5,2)); dy=np.array([[10.,1.],[20.,1.],[1.,1.],[4.,1.],[2.,1.]])
        parts=dict(DY=dy,TTplusWW=z,VZ=z,Data=np.tile([10.,2.],(5,1)))
        mc,d,mask,rows=calculate(parts)
        np.testing.assert_array_equal(mask,[True,False])
        self.assertEqual(rows[0]['data'],12)

    def test_zero_data_upper_error(self):
        e=poisson_errors(np.array([0.,10.]))
        self.assertEqual(e[0,0],0); self.assertGreater(e[1,0],1)


if __name__=='__main__': unittest.main()
