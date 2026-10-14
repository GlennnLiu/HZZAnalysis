"""Small deterministic tests for the independent-region closure calculation."""
import unittest
import os
for variable in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[variable] = '1'
import numpy as np
from auxiliary_alpha_closure import MET, closure, region


class ClosureTests(unittest.TestCase):
    def test_unity_and_independent_variance(self):
        c = [np.array([v]) for v in [10., 10., 100., 100.]]
        t = [np.array([v]) for v in [20., 20., 200., 200.]]
        a, p, vp, k, vk = closure(c, t)
        np.testing.assert_allclose([a[0], p[0], vp[0], k[0], vk[0]],
                                   [2., 20., 46., 1., .165])

    def test_negative_target_is_not_clipped(self):
        c = [np.array([v]) for v in [10., 10., 100., 100.]]
        t = [np.array([v]) for v in [-5., 20., 200., 200.]]
        self.assertEqual(closure(c, t)[3][0], -.25)

    def test_nonpositive_source_is_flagged(self):
        c = [np.array([v]) for v in [-1., 10., 100., 100.]]
        t = [np.array([v]) for v in [20., 20., 200., 200.]]
        self.assertTrue(np.isnan(closure(c, t)[3][0]))

    def test_subtraction_and_sumw2(self):
        cells = {p: np.array([[[y, v]]]) for p, y, v in
                 [('Data', 20., 20.), ('TTplusWW', 4., 2.),
                  ('VZ', 1., .5), ('DY', 10., 3.)]}
        z, vz, d, vd, _ = region(cells, [0])
        np.testing.assert_allclose([z[0], vz[0], d[0], vd[0]], [15., 22.5, 10., 3.])

    def test_met_boundary(self):
        np.testing.assert_array_equal(np.searchsorted(MET, [150., 150.01], side='left'), [2, 3])


if __name__ == '__main__':
    unittest.main()
