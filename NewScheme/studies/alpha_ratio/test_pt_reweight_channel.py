import unittest
import numpy as np
from pt_reweight_core import fit_pt, predictions, Estimate
from test_pt_reweight import fixture
from pt_reweight_channel import channel_predictions, pt_predictions, common_pt_groups, channel_mass_groups


class TestChannelPTReweight(unittest.TestCase):
    def test_own_channel_closure_and_isolation(self):
        space = fixture()
        for ch in (0, 1):
            fit = fit_pt(space, [1], [0, 2], channels=[ch])
            result = pt_predictions(space, 1, [ch], [0, 2], fit['groups'], {ch: fit})
            np.testing.assert_allclose(result['dy_reweighted'].value, result['observed'].value)
            np.testing.assert_allclose(result['dy_reweighted'].jac, result['observed'].jac, atol=1e-12)
            for proc in (0, 1):
                self.assertTrue(np.all(fit['factor'].jac[:, space.indices(proc, [1], [1-ch], [0, 1, 2])] == 0))
                self.assertTrue(np.all(fit['factor'].jac[:, space.indices(proc, [0, 2], [0, 1], [0, 1, 2])] == 0))

    def test_combined_sum_covariance(self):
        space = fixture()
        # Force genuinely different lepton-channel corrections.
        space.flat[space.indices(0, [1], [1], [0, 2])] *= 1.3
        fits = {ch: fit_pt(space, [1], [0, 2], channels=[ch]) for ch in (0, 1)}
        groups = [np.array([0]), np.array([1])]
        a = channel_predictions(space, 1, [0], 1, groups, fits, (0, 2))
        b = channel_predictions(space, 1, [1], 1, groups, fits, (0, 2))
        combined = channel_predictions(space, 1, [0, 1], 1, groups, fits, (0, 2))
        for key in ('dy_reweighted', 'alpha_reweighted', 'dy_normonly'):
            np.testing.assert_allclose(combined[key].value, a[key].value+b[key].value)
            np.testing.assert_allclose(combined[key].covariance(space.variance),
                a[key].covariance(space.variance)+b[key].covariance(space.variance), atol=1e-10)
        self.assertNotIn('transfer_reweighted', combined)

    def test_no_target_or_other_channel_in_prediction(self):
        space = fixture()
        fits = {ch: fit_pt(space, [1], [0], channels=[ch]) for ch in (0, 1)}
        groups = [np.array([0]), np.array([1])]
        result = channel_predictions(space, 1, [0], 2, groups, fits, (0,))
        for key in ('dy_reweighted', 'alpha_reweighted'):
            self.assertTrue(np.all(result[key].jac[:, space.indices(0, [1], [0], [1, 2])] == 0))
            for proc in (0, 1):
                self.assertTrue(np.all(result[key].jac[:, space.indices(proc, [1], [1], [0, 1, 2])] == 0))

    def test_common_coarsening(self):
        fits = {0: {'groups': [np.array([0, 1]), np.array([2]), np.array([3])]},
                1: {'groups': [np.array([0]), np.array([1]), np.array([2, 3])]}}
        groups = common_pt_groups(fits, [0, 1], 4)
        self.assertEqual([g.tolist() for g in groups], [[0, 1], [2, 3]])

    def test_finite_difference_combined_channel_predictions(self):
        space = fixture()
        pg = [np.array([i]) for i in range(3)]
        def evaluate():
            fits = {}
            for ch in (0, 1):
                z = space.linear(0, [1], [ch], [0, 2], ptgroups=pg)
                dy = space.linear(1, [1], [ch], [0, 2], ptgroups=pg)
                fits[ch] = dict(groups=pg, factor=z/dy, observation=z, dy=dy)
            return channel_predictions(space, 1, [0, 1], 1,
                                       [np.array([0]), np.array([1])], fits, (0, 2))
        result = evaluate()
        for index in np.random.default_rng(9).choice(space.flat.size, 24, replace=False):
            value, epsilon = space.flat[index], 1e-3
            space.flat[index] = value+epsilon
            up = evaluate()
            space.flat[index] = value-epsilon
            down = evaluate()
            space.flat[index] = value
            for key in ('dy_reweighted', 'alpha_reweighted', 'dy_normonly'):
                np.testing.assert_allclose((up[key].value-down[key].value)/(2*epsilon),
                                           result[key].jac[:, index], rtol=2e-6, atol=2e-8)

    def test_constant_channel_factors_cancel_alpha(self):
        space = fixture()
        fits = {ch: fit_pt(space, [1], [0, 2], channels=[ch]) for ch in (0, 1)}
        for ch, fit in fits.items():
            fit['factor'] = Estimate(np.full(len(fit['groups']), 1.2+.4*ch), np.zeros_like(fit['factor'].jac))
        result = channel_predictions(space, 1, [0, 1], 1,
                                     [np.array([0]), np.array([1])], fits, (0, 2))
        np.testing.assert_allclose(result['alpha_nominal'].value, result['alpha_reweighted'].value)
        np.testing.assert_allclose(result['alpha_nominal'].jac, result['alpha_reweighted'].jac, atol=1e-12)


if __name__ == '__main__':
    unittest.main()
