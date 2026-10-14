import unittest
import numpy as np
from check_alpha_closure import prediction,statistics


class Closure(unittest.TestCase):
    def test_variance_components(self):
        a=np.array([[.5],[.6],[.01],[.02],[.005],[1.2],[.001],[.02]])
        p=prediction(a,np.array([[100.],[120.]]),np.array([[20.],[4.]]))
        np.testing.assert_allclose(p[:,0],[70,80,134,247.2,90])

    def test_shared_prediction_cancellation(self):
        a=np.array([[.5],[.5],[.01],[.01],[.01],[1.],[0.],[0.]])
        p=prediction(a,np.array([[100.],[120.]]),np.array([[20.],[4.]]))
        np.testing.assert_allclose(p[2]+p[3]-2*p[4],0)

    def test_statistic(self):
        p=np.array([[10.],[12.],[1.],[2.],[1.]])
        r=statistics(p,np.array([[10.],[10.]]),np.array([True]))
        self.assertEqual(r[0]['chi2'],0); self.assertAlmostEqual(r[1]['chi2'],4./12)


if __name__=='__main__': unittest.main()
