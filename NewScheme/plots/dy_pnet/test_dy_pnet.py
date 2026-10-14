import unittest
import numpy as np
from plot_dy_pnet import ratio


class PairedRatio(unittest.TestCase):
    def test_constant_scale_zero_ratio_variance(self):
        r,e=ratio(np.array([[3.],[3.6],[5.],[7.2],[6.]]))
        np.testing.assert_allclose(r,[1.2])
        np.testing.assert_allclose(e,[0.],atol=1e-8)

    def test_varying_scale_covariance(self):
        # Weights [1,2], SF [1,2]: same-event ratio error, not independent samples.
        r,e=ratio(np.array([[3.],[5.],[5.],[17.],[9.]]))
        np.testing.assert_allclose(r,[5./3])
        np.testing.assert_allclose(e,[np.sqrt(8./81)])

    def test_nonpositive_denominator_mask(self):
        r,e=ratio(np.array([[0.,-1.],[0.,-2.],[0.,1.],[0.,4.],[0.,2.]]))
        self.assertTrue(np.isnan(r).all() and np.isnan(e).all())


if __name__=='__main__': unittest.main()
