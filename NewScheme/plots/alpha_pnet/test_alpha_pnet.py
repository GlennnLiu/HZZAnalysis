import unittest
import numpy as np
from check_alpha_pnet import paired_alpha


class PairedAlpha(unittest.TestCase):
    def test_equal_constant_sf_cancels(self):
        t=np.array([[3.],[6.],[5.],[20.],[10.]])
        c=np.array([[5.],[10.],[13.],[52.],[26.]])
        v,mask=paired_alpha(t,c)
        np.testing.assert_allclose(v[0],v[1]); np.testing.assert_allclose(v[5],1)
        np.testing.assert_allclose(v[6:],0,atol=1e-14)

    def test_region_specific_constant_sf(self):
        t=np.array([[3.],[6.],[5.],[20.],[10.]])
        c=np.array([[5.],[15.],[13.],[117.],[39.]])
        v,mask=paired_alpha(t,c)
        np.testing.assert_allclose(v[5],2./3)
        np.testing.assert_allclose(v[6],0,atol=1e-14)

    def test_invalid_denominator_mask(self):
        h=np.zeros((5,1)); v,mask=paired_alpha(h,h)
        self.assertFalse(mask[0]); self.assertTrue(np.isnan(v).all())


if __name__=='__main__': unittest.main()
