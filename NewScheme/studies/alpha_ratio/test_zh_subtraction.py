"""Analytic regression checks; no real data or production writes."""
import unittest
import numpy as np
from compare_zh_subtraction import derive, combine, hist


def fixture():
    h = {}
    for reg in ['CR', 'target']:
        for proc, values in [('Data', [100., 80.]), ('DY', [50., 20.]),
                             ('TTplusWW', [10., 5.]), ('VZ', [5., 3.]),
                             ('ZHbb', [4., 2.]), ('ZHbb_unit_tag', [4., 2.])]:
            h[reg, proc] = np.stack([values, np.array(values)*.1])
    h['target', 'DY'][0] *= 2
    return h


class TestStudy(unittest.TestCase):
    def test_subtraction_identity(self):
        r = derive(fixture())
        np.testing.assert_allclose(r['pred_nozh'][0]-r['pred_withzh'][0], r['transferred_zh'][0])
        np.testing.assert_allclose(r['obs_nozh'][0]-r['obs_withzh'][0], r['target', 'ZHbb'][0])
        np.testing.assert_allclose(r['pred_total_withzh'][0]-r['pred_nozh'][0], r['delta_total_prediction'][0])
        np.testing.assert_allclose(r['residual_withzh'][0]-r['residual_nozh'][0], r['delta_closure_residual'][0])

    def test_correlated_difference_variance(self):
        r = derive(fixture())
        v = r['pred_nozh'][1]+r['pred_withzh'][1]-2*r['cov_before_after_prediction']
        np.testing.assert_allclose(v, r['delta_dy_prediction'][1])

    def test_combined_is_channel_sum(self):
        r = derive(fixture()); second = derive(fixture())
        c = combine([r, second])
        np.testing.assert_allclose(c['pred_withzh'], 2*r['pred_withzh'])
        np.testing.assert_allclose(c['alpha'], np.stack([r['alpha'][0], r['alpha'][1]/2]))

    def test_signed_weights_and_entire_tail(self):
        a = dict(mass=np.array([299., 300., 3500., 10000.]), weight=np.array([20., 2., -1., 3.]))
        h = hist(a, [300., 3500., np.inf])
        np.testing.assert_allclose(h, [[2., 2.], [4., 10.]])

    def test_nonpositive_dy_rejected(self):
        h = fixture(); h['CR', 'DY'][0, 1] = 0
        with self.assertRaises(ValueError): derive(h)


if __name__ == '__main__': unittest.main()
